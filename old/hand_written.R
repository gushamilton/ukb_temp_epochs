library(tidyverse)
library(arrow)
library(data.table)
library(lubridate)
library(RcppRoll)
library(patchwork)
library(lme4)
library(ICC)
library(corrplot)

#bring in some subset of data

COVARIATE_FILE_PATH <- "combined_covariates.csv.gz"
TEMPERATURE_DATA_PATH <- "combined_temperature_data.parquet"
OUTPUT_DIR <- "results_test_plots"



# --- TEST SECTION: Select top 50 IDs directly from the Parquet file ---
print("--- Running in TEST MODE ---")
temp_dataset <- open_dataset(TEMPERATURE_DATA_PATH)



all_ids_to_process <- temp_dataset %>%
  distinct(id) %>%
  head(50) %>%
  collect() %>%
  pull(id)

raw_temp_data <- temp_dataset %>%
  filter(id %in% all_ids_to_process) %>%
  select(id, time, temp, enmoTrunc) %>%
  collect() %>%
  mutate(time = as_datetime(time))


# Do some filtering


covars <- fread(COVARIATE_FILE_PATH)
participant_counts <- covars %>% count(participant_id)
qc_passed  <- covars %>%
  filter(goodWear == 1, goodCal == 1, wearDays >= 4)


processed_data <- raw_temp_data %>%
filter(id %in% qc_passed$id)  %>%
  # Ensure time is a datetime object
  mutate(time = as_datetime(time)) %>%
  # Group by each participant to apply rolling functions correctly
  group_by(id) %>%
  # Step B: Z-score Filtering to remove local spikes
  mutate(
    rolling_mean = RcppRoll::roll_mean(temp, n = 21, fill = NA, align = "center"),
    rolling_sd = RcppRoll::roll_sd(temp, n = 21, fill = NA, align = "center"),
    z_score = (temp - rolling_mean) / rolling_sd,
    temp_z_filtered = if_else(abs(z_score) > 3, NA_real_, temp)
  ) %>%
  # Step C: Smoothing the Z-score filtered data
  mutate(
    temp_smoothed = RcppRoll::roll_median(
      temp_z_filtered,
      n = 7,
      fill = NA,
      align = "center"
    )
  ) %>%
  ungroup() # Ungroup when finished


# --- 4. Downsample Data to 5-Minute Epochs ---
# THIS IS THE KEY PERFORMANCE IMPROVEMENT STEP
cat("Downsampling data to 5-minute epochs...\n")

downsampled_data <- processed_data %>%
  # Remove rows where smoothed temp is NA before aggregation
  filter(!is.na(temp_smoothed)) %>%
  group_by(id) %>%
  # Create a new time column representing the start of each 5-minute interval
  mutate(time_5min = floor_date(time, unit = "5 minutes")) %>%
  group_by(id, time_5min) %>%
  # Calculate the mean of the smoothed temperature for each 5-minute window
  summarise(
    temp_final = mean(temp_smoothed, na.rm = TRUE),
    .groups = 'drop' # Drop grouping structure after summarising
  )

cat("Downsampling complete. Data size reduced significantly.\n")




# Prepare data for the MLM
model_data_all <- downsampled_data %>%
  mutate(
    hour_of_day = hour(time_5min) + minute(time_5min) / 60,
    time_radian = 2 * pi * hour_of_day / 24
  )





# --- 6. Loop Through Four Random Participants for Visualization ---
set.seed(42) # for reproducibility
random_ids <- sample(unique(processed_data$id), 4)

for (i in 1:length(random_ids)) {
  
  current_id <- random_ids[i]
  cat("\n--- Generating plot", i, "for participant:", current_id, "---\n")
  
  # Use the original processed data for the first three plots to show the steps
  participant_data_vis <- processed_data %>% filter(id == current_id)
  
  # Use the downsampled data for the modeling plot
  model_data_single <- model_data_all %>% filter(id == current_id)
  
  # --- Fit Simple Circadian Model (LM) for this one participant on their downsampled data ---
  lm_model <- lm(temp_final ~ sin(time_radian) + cos(time_radian), data = model_data_single)
  model_data_single$temp_predicted_lm <- predict(lm_model)
  
  # --- Get Predictions from the Global MLM for this participant ---
  # CORRECTED: Removed allow.new.levels = TRUE to ensure random effects are used
  
  # --- Create the Plots ---
  plot_a <- ggplot(participant_data_vis, aes(x = time, y = temp)) +
    geom_line(color = "gray50", alpha = 0.8) +
    labs(title = "A: Raw Data", y = "Temp (°C)", x = "") +
    theme_bw()
  
  plot_b <- ggplot(participant_data_vis, aes(x = time)) +
    geom_line(aes(y = temp_z_filtered), color = "dodgerblue") +
    geom_point(data = . %>% filter(is.na(temp_z_filtered)), aes(y = temp), color = "red", size = 0.5) +
    labs(title = "B: Z-Score Filtered", y = "Temp (°C)", x = "") +
    theme_bw()
  
  plot_c <- ggplot(participant_data_vis, aes(x = time, y = temp_smoothed)) +
    geom_line(color = "navy") +
    labs(title = "C: Smoothed Signal", y = "Temp (°C)", x = "Time") +
    theme_bw()
  
  plot_d <- ggplot(model_data_single, aes(x = time_5min)) +
    geom_point(aes(y = temp_final), color = "gray70", alpha = 0.5, size = 0.8) +
    geom_line(aes(y = temp_predicted_lm, color = "Simple Model (LM)"), linewidth = 1.2) +
    scale_color_manual(values = c("Simple Model (LM)" = "darkorange")) +
    labs(
      title = "D: Fitted Circadian Models on 5-min Data",
      subtitle = "Comparing individual (LM)",
      y = "Temp (°C)",
      x = "Time",
      color = "Model Type"
    ) +
    theme_bw() +
    theme(legend.position = "bottom")
  
  # --- Combine and Save the Plot ---
  final_plot <- (plot_a | plot_b) / (plot_c | plot_d) +
    plot_annotation(title = paste("Preprocessing & Modeling Steps for Participant:", current_id))
  
  file_name <- paste0("participant_", i, "_visualization.pdf")
  ggsave(file_name, final_plot, width = 12, height = 9)
  
  cat("Saved plot to", file_name, "\n")
}


# Function to compute first peak amplitude and frequency for a participant
get_first_peak <- function(df) {
  df <- df %>% arrange(time)
  temp_adjust <- df$temp - mean(df$temp, na.rm = TRUE)
  N <- length(temp_adjust)
  if (N < 10) return(data.frame(id = unique(df$id), first_peak_freq = NA, first_peak_amp = NA)) # skip too-short series
  ft <- fft(temp_adjust)
  freq <- (0:(N-1)) / N
  amplitude <- Mod(ft)
  half <- 2:floor(N/2) # skip DC (freq=0) at index 1

  # Find first local maximum (peak) in amplitude spectrum (excluding DC)
  amp_half <- amplitude[half]
  freq_half <- freq[half]
  # Find local maxima
  peak_idx <- which(diff(sign(diff(amp_half))) == -2) + 1
  if (length(peak_idx) == 0) {
    # fallback: take global max (excluding DC)
    max_idx <- which.max(amp_half)
    first_peak <- max_idx
  } else {
    first_peak <- peak_idx[1]
  }
  first_peak_freq <- freq_half[first_peak]
  first_peak_amp <- amp_half[first_peak]
  data.frame(id = unique(df$id), first_peak_freq = first_peak_freq, first_peak_amp = first_peak_amp)
}

# Apply to all participants
fourier_peaks <- processed_data %>%
  group_by(id) %>%
  group_split() %>%
  map_df(get_first_peak)

# Plot distribution of first peak frequency
ggplot(fourier_peaks, aes(x = first_peak_freq)) +
  geom_histogram(bins = 50, fill = "skyblue", color = "black") +
  labs(
    title = "Distribution of First Peak Frequency (Fourier) Across Participants",
    x = "Frequency (cycles per sample)",
    y = "Count"
  ) +
  theme_bw()

# Plot distribution of first peak amplitude
ggplot(fourier_peaks, aes(x = first_peak_amp)) +
  geom_histogram(bins = 50, fill = "orange", color = "black") +
  labs(
    title = "Distribution of First Peak Amplitude (Fourier) Across Participants",
    x = "Amplitude",
    y = "Count"
  ) +
  theme_bw()

# --- 5. Align All Participants to 24-Hour Clock and Plot Mean + SE ---

# Add hour of day (0-23) for each sample
processed_data <- processed_data %>%
  mutate(hour = hour(time) + minute(time)/60)

# Compute mean and standard error for each hour across all participants
hourly_summary <- processed_data %>%
  group_by(hour) %>%
  summarise(
    mean_temp = mean(temp_smoothed, na.rm = TRUE),
    se_temp = sd(temp_smoothed, na.rm = TRUE) / sqrt(sum(!is.na(temp_smoothed)))
  )

# Plot mean temperature over 24 hours with standard error ribbon
ggplot(hourly_summary, aes(x = hour, y = mean_temp)) +
  geom_line(color = "blue", size = 1) +
  geom_ribbon(aes(ymin = mean_temp - se_temp, ymax = mean_temp + se_temp), alpha = 0.2, fill = "blue") +
  scale_x_continuous(breaks = seq(0, 24, by = 2), limits = c(0, 24)) +
  labs(
    title = "Average Smoothed Temperature Across 24-Hour Clock",
    x = "Hour of Day",
    y = "Mean Temperature (smoothed, °C)"
  ) +
  theme_bw()

# --- 5b. Align All Participants to 24-Hour Clock and Plot Mean + SE for Activity ---

# Compute mean and standard error for each hour across all participants for activity
hourly_activity_summary <- processed_data %>%
  group_by(hour) %>%
  summarise(
    mean_enmo = mean(enmoTrunc, na.rm = TRUE),
    se_enmo = sd(enmoTrunc, na.rm = TRUE) / sqrt(sum(!is.na(enmoTrunc)))
  )

# Plot mean activity over 24 hours with standard error ribbon
ggplot(hourly_activity_summary, aes(x = hour, y = mean_enmo)) +
  geom_line(color = "darkgreen", size = 1) +
  geom_ribbon(aes(ymin = mean_enmo - se_enmo, ymax = mean_enmo + se_enmo), alpha = 0.2, fill = "darkgreen") +
  scale_x_continuous(breaks = seq(0, 24, by = 2), limits = c(0, 24)) +
  labs(
    title = "Average Activity (ENMO) Across 24-Hour Clock",
    x = "Hour of Day",
    y = "Mean ENMO (smoothed, g)"
  ) +
  theme_bw()
# --- 5c. Plot Inverse Activity on Temperature Plot ---

# First, rescale mean_enmo to match the range of mean_temp for overlay
# We'll use a simple min-max scaling to the range of mean_temp
activity_min <- min(hourly_activity_summary$mean_enmo, na.rm = TRUE)
activity_max <- max(hourly_activity_summary$mean_enmo, na.rm = TRUE)
temp_min <- min(hourly_summary$mean_temp, na.rm = TRUE)
temp_max <- max(hourly_summary$mean_temp, na.rm = TRUE)

# Invert activity (so high activity = low value) and rescale
hourly_activity_summary <- hourly_activity_summary %>%
  mutate(
    mean_enmo_inv = activity_max - mean_enmo, # invert
    mean_enmo_inv_rescaled = (mean_enmo_inv - min(mean_enmo_inv, na.rm = TRUE)) / 
                             (max(mean_enmo_inv, na.rm = TRUE) - min(mean_enmo_inv, na.rm = TRUE)) * 
                             (temp_max - temp_min) + temp_min
  )

# Plot temperature and inverse activity on the same plot
ggplot() +
  geom_line(data = hourly_summary, aes(x = hour, y = mean_temp), color = "blue", size = 1) +
  geom_ribbon(data = hourly_summary, aes(x = hour, ymin = mean_temp - se_temp, ymax = mean_temp + se_temp), 
              alpha = 0.2, fill = "blue") +
  geom_line(data = hourly_activity_summary, aes(x = hour, y = mean_enmo_inv_rescaled), color = "red", size = 1, linetype = "dashed") +
  scale_x_continuous(breaks = seq(0, 24, by = 2), limits = c(0, 24)) +
  labs(
    title = "Temperature and Inverse Activity Across 24-Hour Clock",
    x = "Hour of Day",
    y = "Mean Temperature (smoothed, °C)",
    caption = "Red dashed line: Inverse Activity (rescaled to temp axis)"
  ) +
  theme_bw()

  # --- 5d. Correlation Between Temperature and Activity ---

  # 1. Correlation between mean temp and mean activity (averaged across all participants, per hour)
  avg_temp_activity_corr <- cor(
    hourly_summary$mean_temp,
    hourly_activity_summary$mean_enmo,
    use = "complete.obs"
  )
  cat("Correlation between mean temperature and mean activity (across all participants, per hour):", round(avg_temp_activity_corr, 3), "\n")

  # 2. For a random sample of 100 individuals, calculate correlation between temp and activity per hour
  set.seed(42)
  # Assume 'qc_data' is the main data frame with columns: id, hour, temp, enmoTrunc
  # If your data frame is named differently, adjust accordingly

  # Get 100 random participant IDs
  sampled_ids <- sample(unique(processed_data$id), 10)

  # For each participant, calculate correlation between temp and activity per hour
  indiv_corrs <- sapply(sampled_ids, function(pid) {
    df <- processed_data %>%
      filter(id == pid) %>%
      group_by(hour) %>%
      summarise(
        mean_temp = mean(temp, na.rm = TRUE),
        mean_enmo = mean(enmoTrunc, na.rm = TRUE)
      )
    # Only compute correlation if there are at least 5 hours with data
    if (nrow(df) >= 5) {
      cor(df$mean_temp, df$mean_enmo, use = "complete.obs")
    } else {
      NA
    }
  })

  # Present the mean and distribution of these correlations
  mean_indiv_corr <- mean(indiv_corrs, na.rm = TRUE)
  cat("Mean correlation between temperature and activity per hour (across 100 individuals):", round(mean_indiv_corr, 3), "\n")

  # Plot the distribution
  indiv_corrs_df <- data.frame(correlation = indiv_corrs)
  ggplot(indiv_corrs_df, aes(x = correlation)) +
    geom_histogram(binwidth = 0.05, fill = "skyblue", color = "black") +
    labs(
      title = "Distribution of Per-Individual Temp-Activity Correlations (per hour)",
      x = "Correlation (temp vs activity, per hour)",
      y = "Count"
    ) +
    theme_bw()


    # --- Compare Fourier Transform Amplitude and Frequency Across Individuals ---

    # 1. Compute mean and standard deviation of first peak frequency and amplitude across all participants
    mean_peak_freq <- mean(fourier_peaks$first_peak_freq, na.rm = TRUE)
    sd_peak_freq <- sd(fourier_peaks$first_peak_freq, na.rm = TRUE)
    mean_peak_amp <- mean(fourier_peaks$first_peak_amp, na.rm = TRUE)
    sd_peak_amp <- sd(fourier_peaks$first_peak_amp, na.rm = TRUE)

    cat("Mean first peak frequency (across individuals):", round(mean_peak_freq, 4), 
        "SD:", round(sd_peak_freq, 4), "\n")
    cat("Mean first peak amplitude (across individuals):", round(mean_peak_amp, 4), 
        "SD:", round(sd_peak_amp, 4), "\n")

    # 2. Correlation between first peak frequency and amplitude across individuals
    freq_amp_corr <- cor(fourier_peaks$first_peak_freq, fourier_peaks$first_peak_amp, use = "complete.obs")
    cat("Correlation between first peak frequency and amplitude (across individuals):", 
        round(freq_amp_corr, 3), "\n")

    # 3. Plot joint distribution (scatterplot) of first peak frequency vs amplitude
    ggplot(fourier_peaks, aes(x = first_peak_freq, y = first_peak_amp)) +
      geom_point(alpha = 0.6, color = "purple") +
      labs(
        title = "First Fourier Peak: Frequency vs Amplitude (per individual)",
        x = "First Peak Frequency (cycles per sample)",
        y = "First Peak Amplitude"
      ) +
      theme_bw()

      # Calculate the Fourier transform first peak for activity (enmoTrunc) for each participant
      get_first_peak_activity <- function(df) {
        df <- df %>% arrange(time)
        enmo_adjust <- df$enmoTrunc - mean(df$enmoTrunc, na.rm = TRUE)
        N <- length(enmo_adjust)
        if (N < 10) return(data.frame(id = unique(df$id), first_peak_freq_enmo = NA, first_peak_amp_enmo = NA))
        ft <- fft(enmo_adjust)
        freq <- (0:(N-1)) / N
        amplitude <- Mod(ft)
        half <- 2:floor(N/2)
        amp_half <- amplitude[half]
        freq_half <- freq[half]
        peak_idx <- which(diff(sign(diff(amp_half))) == -2) + 1
        if (length(peak_idx) == 0) {
          max_idx <- which.max(amp_half)
          first_peak <- max_idx
        } else {
          first_peak <- peak_idx[1]
        }
        first_peak_freq_enmo <- freq_half[first_peak]
        first_peak_amp_enmo <- amp_half[first_peak]
        data.frame(id = unique(df$id), first_peak_freq_enmo = first_peak_freq_enmo, first_peak_amp_enmo = first_peak_amp_enmo)
      }

      # Apply to all participants
      fourier_peaks_enmo <- processed_data %>%
        group_by(id) %>%
        group_split() %>%
        map_df(get_first_peak_activity)

      # Merge temperature and activity Fourier peak data by id
      fourier_peaks_merged <- fourier_peaks %>%
        left_join(fourier_peaks_enmo, by = "id")

      # Correlate amplitude of first peak of temperature with amplitude of first peak of activity
      amp_corr <- cor(fourier_peaks_merged$first_peak_amp, fourier_peaks_merged$first_peak_amp_enmo, use = "complete.obs")
      cat("Correlation between amplitude of first Fourier peak (temperature vs activity):", round(amp_corr, 3), "\n")

      # Optionally, plot the relationship
      ggplot(fourier_peaks_merged, aes(x = first_peak_amp_enmo, y = first_peak_amp)) +
        geom_point(alpha = 0.6, color = "darkred") +
        labs(
          title = "First Fourier Peak Amplitude: Temperature vs Activity (per individual)",
          x = "First Peak Amplitude (Activity, ENMO)",
          y = "First Peak Amplitude (Temperature)"
        ) +
        theme_bw()


        # Plot 4 participants: individual inverse temp vs activity

        set.seed(123)
        # Sample 4 random participant IDs
        random_ids_inv <- sample(unique(processed_data$id), 4)

        for (i in seq_along(random_ids_inv)) {
          pid <- random_ids_inv[i]
          df <- processed_data %>% filter(id == pid)

          # Inverse temperature (z-score for better comparability)
          df <- df %>%
            mutate(
              inv_temp = -scale(temp_smoothed),
              enmo_z = scale(enmoTrunc)
            )

          p <- ggplot(df, aes(x = time)) +
            geom_line(aes(y = inv_temp, color = "Inverse Temp (z)"), size = 1) +
            geom_line(aes(y = enmo_z, color = "Activity (ENMO, z)"), size = 1, alpha = 0.7) +
            scale_color_manual(values = c("Inverse Temp (z)" = "blue", "Activity (ENMO, z)" = "darkgreen")) +
            labs(
              title = paste("Participant", pid, "- Inverse Temp vs Activity"),
              x = "Time",
              y = "z-score",
              color = ""
            ) +
            theme_bw() +
            theme(legend.position = "top")

          ggsave(sprintf("participant_%s_invtemp_vs_activity.pdf", pid), p, width = 10, height = 4)
          print(p)
        }

        # Overlay all individuals, aligning their time so that each starts at 0
        # (i.e., time since their own first sample)

        # Compute time since first sample for each participant, remove extreme temp_smoothed (Z > 5), plot on real scale
        processed_data_aligned <- processed_data %>%
          group_by(id) %>%
          arrange(time) %>%
          mutate(
            time_since_start = as.numeric(difftime(time, min(time), units = "hours")),
            temp_z = as.numeric(scale(temp_smoothed))
          ) %>%
          ungroup() %>%
          filter(!is.na(temp_smoothed), abs(temp_z) <= 5)

        # Plot all individuals' temperature (real scale), aligned at time 0
        p_all <- ggplot(processed_data_aligned, aes(x = time_since_start, y = temp_smoothed, group = id)) +
          geom_line(alpha = 0.3, color = "blue") +
          labs(
            title = "All Participants: Temperature (Aligned at Start, Extreme Z > 5 Removed)",
            x = "Time Since Start (hours)",
            y = "Temperature (°C)"
          ) +
          theme_bw()

        ggsave("all_participants_temperature_aligned.pdf", p_all, width = 12, height = 6)
        print(p_all)

        # Calculate the average residual temperature after removing activity effects over 5 days

        # First, fit a linear model of temp_smoothed ~ enmoTrunc for each participant
        # and get the residuals (i.e., temp not explained by activity)
        processed_data_resid <- processed_data %>%
          filter(!is.na(temp_smoothed), !is.na(enmoTrunc)) %>%
          group_by(id) %>%
          mutate(
            # Fit linear model per participant and extract residuals
            temp_resid = {
              fit <- lm(temp_smoothed ~ enmoTrunc)
              resid(fit)
            }
          ) %>%
          ungroup()

        # Align time for each participant (time since their own start, in hours)
        processed_data_resid <- processed_data_resid %>%
          group_by(id) %>%
          arrange(time) %>%
          mutate(
            time_since_start = as.numeric(difftime(time, min(time), units = "hours"))
          ) %>%
          ungroup()

        # Limit to first 5 days (120 hours)
        processed_data_resid_5d <- processed_data_resid %>%
          filter(time_since_start <= 120)

        # Bin time into, e.g., 30-minute intervals for averaging
        processed_data_resid_5d <- processed_data_resid_5d %>%
          mutate(
            time_bin = floor(time_since_start * 2) / 2  # 0.5 hour bins
          )

        # Calculate mean, SEM, and SE of residual temp at each time_bin
        summary_resid <- processed_data_resid_5d %>%
          group_by(time_bin) %>%
          summarise(
            mean_resid = mean(temp_resid, na.rm = TRUE),
            sem_resid = sd(temp_resid, na.rm = TRUE) / sqrt(n()),  # SEM (standard error of the mean)
            se_resid = sd(temp_resid, na.rm = TRUE) / sqrt(n()),   # SE (identical to SEM here)
            n = n()
          ) %>%
          ungroup()

        # Plot the average residual temperature over 5 days, with SE ribbon
        p_resid <- ggplot(summary_resid, aes(x = time_bin, y = mean_resid)) +
          geom_line(color = "purple", size = 1) +
          geom_ribbon(aes(ymin = mean_resid - se_resid, ymax = mean_resid + se_resid), alpha = 0.2, fill = "purple") +
          labs(
            title = "Average Residual Temperature (After Activity Removed) Over 5 Days",
            x = "Time Since Start (hours)",
            y = "Residual Temperature (°C)"
          ) +
          theme_bw()

        ggsave("average_residual_temp_over_5days.pdf", p_resid, width = 12, height = 5)
        print(p_resid)

        # Plot 4 random participants' residual temperature over time

        set.seed(123)  # For reproducibility
        random_ids <- sample(unique(processed_data_resid_5d$id), 4)

        for (i in seq_along(random_ids)) {
          pid <- random_ids[i]
          pdata <- processed_data_resid_5d %>% filter(id == pid)

          # Bin and summarise for this participant
          pdata_summary <- pdata %>%
            group_by(time_bin) %>%
            summarise(
              mean_resid = mean(temp_resid, na.rm = TRUE),
              se_resid = sd(temp_resid, na.rm = TRUE) / sqrt(n()),
              n = n()
            ) %>%
            ungroup()

          p <- ggplot(pdata_summary, aes(x = time_bin, y = mean_resid)) +
            geom_line(color = "darkgreen", size = 1) +
            geom_ribbon(aes(ymin = mean_resid - se_resid, ymax = mean_resid + se_resid), alpha = 0.2, fill = "darkgreen") +
            labs(
              title = paste("Residual Temperature for Participant", pid),
              x = "Time Since Start (hours)",
              y = "Residual Temperature (°C)"
            ) +
            theme_bw()

          fname <- paste0("residual_temp_participant_", i, "_", pid, ".pdf")
          ggsave(fname, p, width = 10, height = 4)
          print(p)
        }

        # Calculate R2 of the linear regression temp_smoothed ~ enmoTrunc (across all data)
        fit_all <- lm(temp ~ enmoTrunc, data = processed_data_resid)
        r2 <- summary(fit_all)$r.squared
        cat(sprintf("R-squared of temp_smoothed ~ enmoTrunc: %.4f\n", r2))

broom::glance(fit_all)
# Calculate R2 of temp_smoothed ~ enmoTrunc for each participant
library(dplyr)
library(ggplot2)
library(broom)

# Remove rows with missing values in temp_smoothed or enmoTrunc
processed_data_resid_clean <- processed_data_resid %>%
  filter(!is.na(temp_smoothed), !is.na(enmoTrunc), !is.na(id))

# For each participant, fit lm and extract R2
r2_by_id <- processed_data_resid_clean %>%
  group_by(id) %>%
  do({
    fit = lm(temp_smoothed ~ enmoTrunc, data = .)
    tibble(r2 = summary(fit)$r.squared)
  }) %>%
  ungroup()

# Print summary statistics
cat(sprintf("Median R2 across participants: %.4f\n", median(r2_by_id$r2, na.rm = TRUE)))
cat(sprintf("Mean R2 across participants: %.4f\n", mean(r2_by_id$r2, na.rm = TRUE)))

# Plot histogram of R2 values
p_r2_hist <- ggplot(r2_by_id, aes(x = r2)) +
  geom_histogram(bins = 30, fill = "skyblue", color = "black") +
  labs(
    title = "Histogram of R-squared: temp_smoothed ~ enmoTrunc (per participant)",
    x = "R-squared",
    y = "Number of Participants"
  ) +
  theme_bw()

ggsave("histogram_r2_temp_smoothed_vs_enmoTrunc.pdf", p_r2_hist, width = 8, height = 5)
print(p_r2_hist)
 