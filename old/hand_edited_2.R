# ==============================================================================
# R SCRIPT FOR WEARABLE DEVICE TEMPERATURE & ACTIVITY ANALYSIS
#
# Description:
# This script processes, analyzes, and visualizes temperature and activity data
# from wearable devices. It performs the following key steps:
#   1. Loads temperature and covariate data.
#   2. Applies quality control filters to select valid participants.
#   3. Preprocesses raw temperature data by filtering outliers and smoothing.
#   4. Downsamples data for modeling and visualization efficiency.
#   5. Visualizes the preprocessing and modeling steps for a sample of participants.
#   6. Analyzes and plots aggregate 24-hour patterns for temperature and activity.
#   7. Investigates the correlation between temperature and activity.
#   8. Uses Fourier analysis to examine periodic signals in the data.
#   9. Models the relationship between temperature and activity, calculating R-squared.
#   10. Calculates standard and activity-adjusted circadian metrics.
#   11. Generates a correlation plot of all derived metrics.
#   12. Saves all plots and tables to a dedicated directory and archives it.
#
# ==============================================================================


# --- 1. SETUP: LOAD LIBRARIES AND DEFINE CONSTANTS ---

# Ensure all necessary packages are installed and loaded.
packages <- c("tidyverse", "arrow", "data.table", "lubridate", "RcppRoll", "patchwork", "broom", "corrplot")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

# --- Configuration ---
# Define file paths and output directories for easy management.
COVARIATE_FILE_PATH   <- "combined_covariates.csv.gz"
TEMPERATURE_DATA_PATH <- "combined_temperature_data.parquet"
OUTPUT_DIR            <- "results_plots" # All plots will be saved here

# Create the output directory if it doesn't already exist.
cat("Creating output directory:", OUTPUT_DIR, "\n")
dir.create(OUTPUT_DIR, showWarnings = FALSE)


# --- 2. DATA LOADING AND INITIAL QC ---

# --- Load Covariate Data and Apply Quality Control (QC) ---
cat("Loading covariate data and applying QC filters...\n")
covars <- fread(COVARIATE_FILE_PATH)

# Filter for participants who meet the quality criteria:
# - goodWear: Sufficient wear time.
# - goodCal: Successful calibration.
# - wearDays: At least 4 days of valid data.
qc_passed_ids <- covars %>%
  filter(goodWear == 1, goodCal == 1, wearDays >= 4) %>%
  pull(id)

cat(sprintf("Found %d participants passing initial QC.\n", length(qc_passed_ids)))


# --- Load Temperature Data for QC-passed Participants ---
# For efficiency, we only load data for the participants who passed QC.
cat("Loading temperature and activity data from Parquet file...\n")
temp_dataset <- open_dataset(TEMPERATURE_DATA_PATH)

# Using a subset for demonstration purposes. Remove `[1:50]` to run on all QC-passed IDs.
raw_data <- temp_dataset %>%
  filter(id %in% qc_passed_ids[1:50]) %>%
  select(id, time, temp, enmoTrunc) %>%
  collect() %>%
  mutate(time = as_datetime(time)) # Ensure 'time' is in datetime format

cat(sprintf("Loaded data for %d participants.\n", length(unique(raw_data$id))))


# --- 3. DATA PREPROCESSING ---

cat("Starting data preprocessing: Z-score filtering and smoothing...\n")

processed_data <- raw_data %>%
  # Group by participant to apply rolling calculations independently.
  group_by(id) %>%
  # Arrange by time to ensure correct rolling window application.
  arrange(time) %>%
  # Step A: Z-score Filtering to remove local spikes/outliers.
  mutate(
    rolling_mean = RcppRoll::roll_mean(temp, n = 21, fill = NA, align = "center"),
    rolling_sd = RcppRoll::roll_sd(temp, n = 21, fill = NA, align = "center"),
    z_score = (temp - rolling_mean) / rolling_sd,
    temp_z_filtered = if_else(abs(z_score) > 3, NA_real_, temp)
  ) %>%
  # Step B: Smoothing the Z-score filtered data using a rolling median.
  mutate(
    temp_smoothed = RcppRoll::roll_median(
      temp_z_filtered,
      n = 7,
      fill = NA,
      align = "center"
    )
  ) %>%
  # Add a column for hour of the day, which will be used for aggregation.
  mutate(hour = hour(time) + minute(time)/60) %>%
  ungroup() # Ungroup for subsequent global operations.

cat("Preprocessing complete.\n")


# --- 4. DOWNSAMPLE DATA FOR EFFICIENCY ---
cat("Downsampling data to 5-minute epochs for modeling...\n")

downsampled_data <- processed_data %>%
  # Remove rows where smoothed temp is NA before aggregation
  filter(!is.na(temp_smoothed)) %>%
  # Create a new time column representing the start of each 5-minute interval
  mutate(time_5min = floor_date(time, unit = "5 minutes")) %>%
  group_by(id, time_5min) %>%
  # Calculate the mean of key variables for each 5-minute window
  summarise(
    temp_final = mean(temp_smoothed, na.rm = TRUE),
    enmo_final = mean(enmoTrunc, na.rm = TRUE),
    hour = first(hour), # Keep the hour for cosinor models
    .groups = 'drop' # Drop grouping structure after summarising
  )

cat("Downsampling complete.\n")


# --- 5. VISUALIZE PREPROCESSING & MODELING STEPS ---

cat("Generating plots to visualize preprocessing and model fit for 4 random participants...\n")
set.seed(42) # for reproducibility
# Ensure we don't try to sample more IDs than available.
n_vis_sample <- min(4, length(unique(processed_data$id)))
random_ids_vis <- sample(unique(processed_data$id), n_vis_sample)

for (i in seq_along(random_ids_vis)) {
  current_id <- random_ids_vis[i]
  cat("  - Plotting for participant:", current_id, "\n")
  
  # Use original processed data for high-resolution preprocessing plots
  participant_data_vis <- processed_data %>% filter(id == current_id)
  # Use downsampled data for the model fit plot
  participant_model_data <- downsampled_data %>% filter(id == current_id)
  
  # --- Create Plots A, B, C (Preprocessing) ---
  plot_a <- ggplot(participant_data_vis, aes(x = time, y = temp)) +
    geom_line(color = "gray50", alpha = 0.8) +
    labs(title = "A: Raw Temperature Data", y = "Temp (°C)", x = "") + theme_bw()
  
  plot_b <- ggplot(participant_data_vis, aes(x = time)) +
    geom_line(aes(y = temp_z_filtered), color = "dodgerblue") +
    geom_point(data = . %>% filter(is.na(temp_z_filtered)), aes(y = temp), color = "red", size = 0.5) +
    labs(title = "B: Z-Score Filtered Data", y = "Temp (°C)", x = "") + theme_bw()
  
  plot_c <- ggplot(participant_data_vis, aes(x = time, y = temp_smoothed)) +
    geom_line(color = "navy") +
    labs(title = "C: Smoothed Signal", y = "Temp (°C)", x = "Time") + theme_bw()

  # --- Create Plot D (Model Fit) ---
  # Fit a simple cosinor model just for this participant's visualization
  lm_model <- lm(temp_final ~ cos(2*pi*hour/24) + sin(2*pi*hour/24), data = participant_model_data)
  participant_model_data$temp_predicted <- predict(lm_model)
  
  plot_d <- ggplot(participant_model_data, aes(x = time_5min)) +
    geom_point(aes(y = temp_final), color = "gray70", alpha = 0.5, size = 0.8) +
    geom_line(aes(y = temp_predicted), color = "darkorange", linewidth = 1.2) +
    labs(title = "D: Cosinor Model Fit on 5-min Data", y = "Temp (°C)", x = "Time") +
    theme_bw()
  
  # --- Combine and Save the Plot ---
  final_vis_plot <- (plot_a | plot_b) / (plot_c | plot_d) +
    plot_annotation(title = paste("Preprocessing & Modeling Steps for Participant:", current_id))
  
  # Save the combined plot to the output directory.
  file_name <- file.path(OUTPUT_DIR, paste0("01_preprocessing_steps_participant_", i, ".pdf"))
  ggsave(file_name, final_vis_plot, width = 12, height = 9)
}


# --- 6. AGGREGATE 24-HOUR PATTERN ANALYSIS ---

cat("Analyzing and plotting aggregate 24-hour patterns...\n")

# --- Temperature Pattern ---
hourly_temp_summary <- processed_data %>%
  filter(!is.na(temp_smoothed)) %>% group_by(hour) %>%
  summarise(mean_val = mean(temp_smoothed, na.rm = TRUE), se_val = sd(temp_smoothed, na.rm = TRUE) / sqrt(n()), .groups = 'drop')
plot_24hr_temp <- ggplot(hourly_temp_summary, aes(x = hour, y = mean_val)) +
  geom_ribbon(aes(ymin = mean_val - se_val, ymax = mean_val + se_val), alpha = 0.2, fill = "blue") +
  geom_line(color = "blue", linewidth = 1) + scale_x_continuous(breaks = seq(0, 24, by = 2), limits = c(0, 24)) +
  labs(title = "Average Smoothed Temperature Across 24-Hour Cycle", x = "Hour of Day", y = "Mean Temperature (°C)") + theme_bw()
ggsave(file.path(OUTPUT_DIR, "02_aggregate_24hr_temperature.pdf"), plot_24hr_temp, width = 10, height = 6)

# --- Activity (ENMO) Pattern ---
hourly_activity_summary <- processed_data %>%
  filter(!is.na(enmoTrunc)) %>% group_by(hour) %>%
  summarise(mean_val = mean(enmoTrunc, na.rm = TRUE), se_val = sd(enmoTrunc, na.rm = TRUE) / sqrt(n()), .groups = 'drop')
plot_24hr_activity <- ggplot(hourly_activity_summary, aes(x = hour, y = mean_val)) +
  geom_ribbon(aes(ymin = mean_val - se_val, ymax = mean_val + se_val), alpha = 0.2, fill = "darkgreen") +
  geom_line(color = "darkgreen", linewidth = 1) + scale_x_continuous(breaks = seq(0, 24, by = 2), limits = c(0, 24)) +
  labs(title = "Average Activity (ENMO) Across 24-Hour Cycle", x = "Hour of Day", y = "Mean ENMO (g)") + theme_bw()
ggsave(file.path(OUTPUT_DIR, "03_aggregate_24hr_activity.pdf"), plot_24hr_activity, width = 10, height = 6)


# --- 7. TEMPERATURE AND ACTIVITY RELATIONSHIP ANALYSIS ---

cat("Analyzing relationship between temperature and activity...\n")

# --- Overlay Plot ---
activity_min <- min(hourly_activity_summary$mean_val, na.rm = TRUE); activity_max <- max(hourly_activity_summary$mean_val, na.rm = TRUE)
temp_min <- min(hourly_temp_summary$mean_val, na.rm = TRUE); temp_max <- max(hourly_temp_summary$mean_val, na.rm = TRUE)
hourly_activity_summary_rescaled <- hourly_activity_summary %>%
  mutate(inv_mean_val = activity_max - mean_val, rescaled_val = (inv_mean_val - min(inv_mean_val, na.rm=T)) / (max(inv_mean_val, na.rm=T) - min(inv_mean_val, na.rm=T)) * (temp_max - temp_min) + temp_min)
plot_overlay <- ggplot() +
  geom_line(data = hourly_temp_summary, aes(x = hour, y = mean_val, color = "Temperature"), linewidth = 1) +
  geom_ribbon(data = hourly_temp_summary, aes(x = hour, ymin = mean_val - se_val, ymax = mean_val + se_val), fill = "blue", alpha = 0.2) +
  geom_line(data = hourly_activity_summary_rescaled, aes(x = hour, y = rescaled_val, color = "Inverse Activity"), linewidth = 1, linetype = "dashed") +
  scale_color_manual(name = "Metric", values = c("Temperature" = "blue", "Inverse Activity" = "red")) +
  scale_x_continuous(breaks = seq(0, 24, by = 2), limits = c(0, 24)) +
  labs(title = "Temperature and Inverse Activity Across 24-Hour Cycle", subtitle = "Inverse activity is rescaled to match the temperature axis range.", x = "Hour of Day", y = "Mean Temperature (°C)") +
  theme_bw() + theme(legend.position = "bottom")
ggsave(file.path(OUTPUT_DIR, "04_overlay_temp_vs_inv_activity.pdf"), plot_overlay, width = 10, height = 6)

# --- Correlation Analysis ---
avg_temp_activity_corr <- cor(hourly_temp_summary$mean_val, hourly_activity_summary$mean_val, use = "complete.obs")
cat(sprintf("Correlation between aggregate hourly mean temperature and activity: %.3f\n", avg_temp_activity_corr))


# --- 8. FOURIER & R-SQUARED ANALYSIS ---

cat("Performing Fourier and R-squared analysis...\n")

# --- Fourier Analysis ---
get_first_peak <- function(data_vec, id_val) {
  if (length(na.omit(data_vec)) < 20) return(data.frame(id = id_val, peak_freq = NA, peak_amp = NA))
  vec_adjusted <- data_vec - mean(data_vec, na.rm = TRUE); vec_adjusted[is.na(vec_adjusted)] <- 0
  N <- length(vec_adjusted); ft <- fft(vec_adjusted); freq <- (0:(N-1)) / N; amplitude <- Mod(ft) / (N/2)
  half_indices <- 2:floor(N/2); amp_half <- amplitude[half_indices]; freq_half <- freq[half_indices]
  peak_indices <- which(diff(sign(diff(amp_half))) == -2) + 1
  first_peak_idx <- if (length(peak_indices) == 0) which.max(amp_half) else peak_indices[1]
  data.frame(id = id_val, peak_freq = freq_half[first_peak_idx], peak_amp = amp_half[first_peak_idx])
}
fourier_peaks <- processed_data %>% group_by(id) %>%
  summarise(temp_peak = list(get_first_peak(temp_smoothed, unique(id))), enmo_peak = list(get_first_peak(enmoTrunc, unique(id))), .groups = 'drop') %>%
  unnest(temp_peak, names_sep = "_temp") %>% unnest(enmo_peak, names_sep = "_enmo") %>% select(id = id_temp, starts_with("peak"))
ggsave(file.path(OUTPUT_DIR, "05_fft_temperature_frequency_dist.pdf"), ggplot(fourier_peaks, aes(x = peak_freq_temp)) + geom_histogram(bins = 50, fill = "skyblue", color = "black") + labs(title = "Distribution of Dominant Temperature Frequencies", x = "Frequency (cycles/sample)", y = "Count") + theme_bw(), width = 8, height = 5)
ggsave(file.path(OUTPUT_DIR, "06_fft_amplitude_correlation.pdf"), ggplot(fourier_peaks, aes(x = peak_amp_enmo, y = peak_amp_temp)) + geom_point(alpha = 0.6, color = "darkred") + labs(title = "Dominant Peak Amplitude: Temperature vs. Activity", x = "Activity Peak Amplitude", y = "Temperature Peak Amplitude") + theme_bw(), width = 8, height = 6)
amp_corr <- cor(fourier_peaks$peak_amp_temp, fourier_peaks$peak_amp_enmo, use = "complete.obs")
cat(sprintf("Correlation between temperature and activity peak amplitudes: %.3f\n", amp_corr))

# --- R-Squared Analysis ---
r2_by_id <- processed_data %>% filter(!is.na(temp_smoothed) & !is.na(enmoTrunc)) %>% group_by(id) %>%
  do(broom::glance(lm(temp_smoothed ~ enmoTrunc, data = .))) %>% ungroup() %>% select(id, r.squared)
cat(sprintf("Median R-squared (temp ~ activity) across participants: %.4f\n", median(r2_by_id$r.squared, na.rm = TRUE)))
ggsave(file.path(OUTPUT_DIR, "07_r_squared_distribution.pdf"), ggplot(r2_by_id, aes(x = r.squared)) + geom_histogram(bins = 30, fill = "purple", color = "black") + labs(title = "Distribution of R-squared (temp ~ activity)", x = "R-squared", y = "Count") + theme_bw(), width = 8, height = 5)


# --- 9. PER-PARTICIPANT METRIC CALCULATION & CORRELATION ---

cat("Calculating summary metrics (Standard & Activity-Adjusted) for each participant...\n")

# Function to fit cosinor models and extract circadian parameters.
calculate_circadian_metrics <- function(df) {
  
  # --- Model 1: Unadjusted Cosinor ---
  unadj_model <- lm(temp_final ~ cos(2*pi*hour/24) + sin(2*pi*hour/24), data = df)
  coeffs_unadj <- coef(unadj_model)
  
  if (any(is.na(coeffs_unadj)) || length(coeffs_unadj) < 3) {
    mesor <- NA; amplitude <- NA; acrophase_h <- NA
  } else {
    mesor <- coeffs_unadj[1]
    amplitude <- sqrt(coeffs_unadj[2]^2 + coeffs_unadj[3]^2)
    acrophase_rad <- atan2(coeffs_unadj[3], coeffs_unadj[2])
    acrophase_h <- (-acrophase_rad * (24 / (2 * pi))) %% 24
  }
  
  # --- Model 2: Activity-Adjusted Cosinor ---
  # Center enmo to make the intercept (mesor) interpretable as the rhythm at the mean activity level.
  df <- df %>% mutate(enmo_centered = enmo_final - mean(enmo_final, na.rm = TRUE))
  
  adj_model <- lm(temp_final ~ enmo_centered + cos(2*pi*hour/24) + sin(2*pi*hour/24), data = df)
  coeffs_adj <- coef(adj_model)
  
  if (any(is.na(coeffs_adj)) || length(coeffs_adj) < 4) {
    mesor_adj <- NA; amplitude_adj <- NA; acrophase_h_adj <- NA
  } else {
    mesor_adj <- coeffs_adj[1] # Intercept is mesor at mean activity
    amplitude_adj <- sqrt(coeffs_adj[3]^2 + coeffs_adj[4]^2)
    acrophase_rad_adj <- atan2(coeffs_adj[4], coeffs_adj[3])
    acrophase_h_adj <- (-acrophase_rad_adj * (24 / (2 * pi))) %% 24
  }

  return(tibble(mesor, amplitude, acrophase_h, mesor_adj, amplitude_adj, acrophase_h_adj))
}

# Apply the function to each participant's downsampled data.
circadian_metrics <- downsampled_data %>%
  group_by(id) %>%
  nest() %>%
  mutate(metrics = map(data, calculate_circadian_metrics)) %>%
  select(id, metrics) %>%
  unnest(metrics)

# --- Combine all metrics into a single summary tibble ---
participant_summary_metrics <- circadian_metrics %>%
  left_join(r2_by_id, by = "id") %>%
  left_join(fourier_peaks, by = "id") %>%
  rename(
    r_squared_activity = r.squared,
    fft_freq_temp = peak_freq_temp,
    fft_amp_temp = peak_amp_temp,
    fft_freq_enmo = peak_freq_enmo,
    fft_amp_enmo = peak_amp_enmo
  )

cat("Summary metrics tibble created:\n")
print(head(participant_summary_metrics))

# Save the summary metrics to a CSV file.
summary_filename <- file.path(OUTPUT_DIR, "participant_summary_metrics.csv")
write_csv(participant_summary_metrics, summary_filename)
cat(sprintf("Metrics summary saved to: %s\n", summary_filename))


# --- Create and Save Correlation Plot ---
cat("Generating correlation plot of all metrics...\n")

metrics_for_corr <- participant_summary_metrics %>%
  select_if(is.numeric) %>%
  select_if(~!all(is.na(.))) # Exclude any all-NA columns

corr_matrix <- cor(metrics_for_corr, use = "pairwise.complete.obs")

# Save the plot directly to a PDF for better control over dimensions.
pdf(file.path(OUTPUT_DIR, "08_metrics_correlation_plot.pdf"), width = 10, height = 10)
corrplot(corr_matrix,
         method = "color",
         type = "upper",
         order = "hclust",
         addCoef.col = "black",
         tl.col = "black", tl.srt = 45,
         diag = FALSE,
         title = "Correlation Matrix of Participant Metrics",
         mar = c(0,0,2,0))
dev.off() # Close the PDF device.

cat("Correlation plot saved.\n")


# --- 10. CREATE FINAL ARCHIVE ---

cat("Archiving all plots and data into a single file...\n")
tar(
  tarfile = paste0(OUTPUT_DIR, ".tar.gz"),
  files = OUTPUT_DIR,
  compression = "gzip"
)
cat(sprintf("Successfully created archive: %s.tar.gz\n", OUTPUT_DIR))

cat("\n--- SCRIPT FINISHED SUCCESSFULLY ---\n")


#add in local data

location <- as.tibble(fread("data_participant.tsv")) 
temps_at_location <- fread("https://github.com/gushamilton/ukb_temp_epochs/raw/refs/heads/main/data/ukb_temps_per_date.tsv.gz") %>% distinct()

linked_id <- location %>%
 select(participant_id = eid, ukb_ac = p54_i0) %>%
 right_join(covars)
 select(participant_id = eid, ukb_ac = p54_i0) %>%
 right_join(qc_passed_id)

temps_at_location %>% head(30) %>%
    mutate(p54_i0 = word(centre_name)) %>%
    
