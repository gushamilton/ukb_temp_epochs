# --- 1. Setup: Load Libraries and Define Paths ---

# Use pacman to install and load packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, arrow, data.table, RcppRoll, patchwork, lubridate, lme4, broom)

# Define relative paths
COVARIATE_FILE_PATH <- "combined_covariates.csv.gz"
TEMPERATURE_DATA_PATH <- "combined_temperature_data.parquet"
OUTPUT_DIR <- "results_test_plots"

# Create the output directory if it doesn't exist
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Define the number of participants to process
N_PARTICIPANTS_TO_LOAD <- 100


# --- 2. Load Covariates and Select a Subset of Participants ---

print("--- Starting: Loading covariates and selecting participants ---")

if (!file.exists(COVARIATE_FILE_PATH)) {
  stop("Covariate file not found. Please check the COVARIATE_FILE_PATH.")
}

# Load the covariate data
covars <- fread(COVARIATE_FILE_PATH)

# Get a list of the first N unique participant IDs from the covariate file
# Note: We are assuming the 'id' column in the covariates matches the 'id' in the Parquet file.
# If the column name is different (e.g., 'participant_id'), you should change it here.
if (!"id" %in% names(covars)) {
    stop("The covariate file must contain an 'id' column for filtering.")
}
selected_participant_ids <- unique(covars$id)[1:N_PARTICIPANTS_TO_LOAD]

print(paste("Selected the first", N_PARTICIPANTS_TO_LOAD, "unique participants for analysis."))
print("--- Finished: Covariate loading and participant selection ---")


# --- 3. Load Temperature Data for Selected Participants Only ---

print("--- Starting: Loading temperature data for the selected subset ---")

if (!file.exists(TEMPERATURE_DATA_PATH)) {
  stop("Temperature data file not found. Please check the TEMPERATURE_DATA_PATH.")
}

# Use the arrow package to efficiently filter the Parquet file before loading into memory
temp_dataset <- open_dataset(TEMPERATURE_DATA_PATH)

# Filter the dataset for the selected participant IDs and then collect the results into an R data.frame
d <- temp_dataset %>%
  filter(id %in% selected_participant_ids) %>%
  collect()

print(paste("Successfully loaded data for", length(unique(d$id)), "participants."))
print(paste("Total rows loaded:", nrow(d)))
print("--- Finished: Data loading ---")


# --- 4. Process and Smooth Data for the Subset ---

cat("Processing all loaded participants...\n")

processed_data <- d %>%
  # Ensure time is a datetime object
  mutate(time = as_datetime(time)) %>%
  # Group by each participant to apply rolling functions correctly
  group_by(id) %>%
  # Step B: Z-score Filtering to remove local spikes
  mutate(
    rolling_mean = RcppRoll::roll_mean(temp, n = 121, fill = NA, align = "center"),
    rolling_sd = RcppRoll::roll_sd(temp, n = 121, fill = NA, align = "center"),
    z_score = (temp - rolling_mean) / rolling_sd,
    temp_z_filtered = if_else(abs(z_score) > 3, NA_real_, temp)
  ) %>%
  # Step C: Smoothing the Z-score filtered data
  mutate(
    temp_smoothed = RcppRoll::roll_median(
      temp_z_filtered,
      n = 21,
      fill = NA,
      align = "center"
    )
  ) %>%
  ungroup() # Ungroup when finished

cat("Processing complete.\n")


# --- 5. Downsample Data to 5-Minute Epochs ---
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

cat("Downsampling complete.\n")


# --- 6. Fit Global Mixed-Effects Model (MLM) on Downsampled Data ---
cat("Fitting global mixed-effects model (lmer)...\n")

# Prepare data for the MLM
model_data_all <- downsampled_data %>%
  mutate(
    hour_of_day = hour(time_5min) + minute(time_5min) / 60,
    time_radian = 2 * pi * hour_of_day / 24
  )

# Fit the MLM with random intercepts and slopes for each participant
mlm_model <- lmer(temp_final ~ sin(time_radian) + cos(time_radian) + (1 + sin(time_radian) + cos(time_radian) | id), data = model_data_all)

cat("MLM fitting complete.\n")


# --- 7. Loop Through Four Random Participants for Visualization ---
set.seed(42) # for reproducibility
# Ensure we don't try to sample more IDs than we have
num_to_sample <- min(4, length(unique(processed_data$id)))
random_ids <- sample(unique(processed_data$id), num_to_sample)

for (i in 1:length(random_ids)) {
  
  current_id <- random_ids[i]
  cat("\n--- Generating plot", i, "for participant:", current_id, "---\n")
  
  participant_data_vis <- processed_data %>% filter(id == current_id)
  model_data_single <- model_data_all %>% filter(id == current_id)
  
  # Fit Simple Circadian Model (LM) for this one participant
  lm_model <- lm(temp_final ~ sin(time_radian) + cos(time_radian), data = model_data_single)
  model_data_single$temp_predicted_lm <- predict(lm_model)
  
  # Get Predictions from the Global MLM for this participant
  model_data_single$temp_predicted_mlm <- predict(mlm_model, newdata = model_data_single)
  
  # Create the Plots
  plot_a <- ggplot(participant_data_vis, aes(x = time, y = temp)) +
    geom_line(color = "gray50", alpha = 0.8) + labs(title = "A: Raw Data", y = "Temp (°C)", x = "") + theme_bw()
  
  plot_b <- ggplot(participant_data_vis, aes(x = time)) +
    geom_line(aes(y = temp_z_filtered), color = "dodgerblue") +
    geom_point(data = . %>% filter(is.na(temp_z_filtered)), aes(y = temp), color = "red", size = 0.5) +
    labs(title = "B: Z-Score Filtered", y = "Temp (°C)", x = "") + theme_bw()
  
  plot_c <- ggplot(participant_data_vis, aes(x = time, y = temp_smoothed)) +
    geom_line(color = "navy") + labs(title = "C: Smoothed Signal", y = "Temp (°C)", x = "Time") + theme_bw()
  
  plot_d <- ggplot(model_data_single, aes(x = time_5min)) +
    geom_point(aes(y = temp_final), color = "gray70", alpha = 0.5, size = 0.8) +
    geom_line(aes(y = temp_predicted_lm, color = "Simple Model (LM)"), linewidth = 1.2) +
    geom_line(aes(y = temp_predicted_mlm, color = "Mixed Model (MLM)"), linewidth = 1.2) +
    scale_color_manual(values = c("Simple Model (LM)" = "darkorange", "Mixed Model (MLM)" = "purple")) +
    labs(
      title = "D: Fitted Circadian Models on 5-min Data",
      subtitle = "Comparing individual (LM) vs. global (MLM) fit",
      y = "Temp (°C)", x = "Time", color = "Model Type"
    ) +
    theme_bw() + theme(legend.position = "bottom")
  
  # Combine and Save the Plot
  final_plot <- (plot_a | plot_b) / (plot_c | plot_d) +
    plot_annotation(title = paste("Preprocessing & Modeling Steps for Participant:", current_id))
  
  file_name <- file.path(OUTPUT_DIR, paste0("participant_", current_id, "_visualization.pdf"))
  ggsave(file_name, final_plot, width = 12, height = 9)
  
  cat("Saved plot to", file_name, "\n")
}


# --- 8. Additional Analyses on the Subset ---

# Note: All the following analyses are performed on the `processed_data` dataframe,
# which contains the subset of 100 participants.

# --- 8a. Fourier Analysis ---
cat("\n--- Starting Fourier Analysis ---\n")

get_first_peak <- function(df) {
  df <- df %>% arrange(time) %>% filter(!is.na(temp_smoothed))
  temp_adjust <- df$temp_smoothed - mean(df$temp_smoothed, na.rm = TRUE)
  N <- length(temp_adjust)
  if (N < 10) return(data.frame(id = unique(df$id), first_peak_freq = NA, first_peak_amp = NA))
  
  ft <- fft(temp_adjust)
  freq <- (0:(N-1)) / N
  amplitude <- Mod(ft)
  half <- 2:floor(N/2)
  
  amp_half <- amplitude[half]
  freq_half <- freq[half]
  peak_idx <- which(diff(sign(diff(amp_half))) == -2) + 1
  
  first_peak <- if (length(peak_idx) == 0) which.max(amp_half) else peak_idx[1]
  
  data.frame(id = unique(df$id), first_peak_freq = freq_half[first_peak], first_peak_amp = amp_half[first_peak])
}

fourier_peaks <- processed_data %>%
  group_by(id) %>%
  group_split() %>%
  map_df(get_first_peak)

# Plot distribution of first peak frequency
p_freq <- ggplot(fourier_peaks, aes(x = first_peak_freq)) +
  geom_histogram(bins = 50, fill = "skyblue", color = "black") +
  labs(title = "Distribution of First Peak Frequency (Fourier)", x = "Frequency", y = "Count") + theme_bw()
ggsave(file.path(OUTPUT_DIR, "fourier_freq_dist.pdf"), p_freq)

# Plot distribution of first peak amplitude
p_amp <- ggplot(fourier_peaks, aes(x = first_peak_amp)) +
  geom_histogram(bins = 50, fill = "orange", color = "black") +
  labs(title = "Distribution of First Peak Amplitude (Fourier)", x = "Amplitude", y = "Count") + theme_bw()
ggsave(file.path(OUTPUT_DIR, "fourier_amp_dist.pdf"), p_amp)


# --- 8b. 24-Hour Averages ---
cat("--- Starting 24-Hour Average Analysis ---\n")

# Add hour of day (0-23.99) for each sample
processed_data_hourly <- processed_data %>%
  mutate(hour = hour(time) + minute(time)/60)

# Compute mean and standard error for each hour across all participants
hourly_summary <- processed_data_hourly %>%
  group_by(hour = round(hour * 2)/2) %>% # Binning to 30 mins for smoother plot
  summarise(
    mean_temp = mean(temp_smoothed, na.rm = TRUE),
    se_temp = sd(temp_smoothed, na.rm = TRUE) / sqrt(sum(!is.na(temp_smoothed))),
    .groups = 'drop'
  )

# Plot mean temperature over 24 hours
p_hourly <- ggplot(hourly_summary, aes(x = hour, y = mean_temp)) +
  geom_line(color = "blue", size = 1) +
  geom_ribbon(aes(ymin = mean_temp - se_temp, ymax = mean_temp + se_temp), alpha = 0.2, fill = "blue") +
  scale_x_continuous(breaks = seq(0, 24, by = 2), limits = c(0, 24)) +
  labs(title = "Average Smoothed Temperature Across 24-Hour Clock", x = "Hour of Day", y = "Mean Temperature (°C)") +
  theme_bw()
ggsave(file.path(OUTPUT_DIR, "hourly_average_temp.pdf"), p_hourly)

# --- 8c. R-squared of Activity vs. Temperature per Participant ---
cat("--- Starting R-squared Analysis (Temp vs Activity) ---\n")

# Check if 'enmoTrunc' column exists
if ("enmoTrunc" %in% names(processed_data)) {
    
    processed_data_resid_clean <- processed_data %>%
      filter(!is.na(temp_smoothed), !is.na(enmoTrunc), !is.na(id))

    # For each participant, fit lm and extract R2
    r2_by_id <- processed_data_resid_clean %>%
      group_by(id) %>%
      do({
        # Ensure there's enough data to fit a model
        if(nrow(.) > 10 && length(unique(.$enmoTrunc)) > 1) {
            fit = lm(temp_smoothed ~ enmoTrunc, data = .)
            tibble(r2 = summary(fit)$r.squared)
        } else {
            tibble(r2 = NA_real_)
        }
      }) %>%
      ungroup()

    # Print summary statistics
    cat(sprintf("Median R2 across participants: %.4f\n", median(r2_by_id$r2, na.rm = TRUE)))
    cat(sprintf("Mean R2 across participants: %.4f\n", mean(r2_by_id$r2, na.rm = TRUE)))

    # Plot histogram of R2 values
    p_r2_hist <- ggplot(r2_by_id, aes(x = r2)) +
      geom_histogram(bins = 30, fill = "skyblue", color = "black", na.rm = TRUE) +
      labs(
        title = "Histogram of R-squared: temp_smoothed ~ enmoTrunc (per participant)",
        x = "R-squared",
        y = "Number of Participants"
      ) +
      theme_bw()

    ggsave(file.path(OUTPUT_DIR, "histogram_r2_temp_vs_activity.pdf"), p_r2_hist, width = 8, height = 5)
    print(p_r2_hist)

} else {
    cat("Skipping R-squared analysis: 'enmoTrunc' column not found in the data.\n")
}

cat("\n--- Analysis script finished successfully! ---\n")
