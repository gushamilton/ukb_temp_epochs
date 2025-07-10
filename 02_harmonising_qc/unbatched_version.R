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
#   5. Links external daily ambient temperature data based on location and date.
#   6. Visualizes the preprocessing and modeling steps for a sample of participants.
#   7. Analyzes and plots aggregate 24-hour patterns for temperature and activity.
#   8. Uses Fourier analysis to examine periodic signals in the data.
#   9. Models the relationship between temperature and activity, calculating R-squared.
#   10. Calculates standard, activity-adjusted, and ambient-adjusted circadian metrics,
#       as well as several novel phenotypes based on the user's research plan.
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
LOCATION_FILE_PATH    <- "data_participant.tsv" # Path to your participant location file
OUTPUT_DIR            <- "results_plots" # All plots will be saved here

# Create the output directory if it doesn't already exist.
cat("Creating output directory:", OUTPUT_DIR, "\n")
dir.create(OUTPUT_DIR, showWarnings = FALSE)


# --- 2. DATA LOADING AND INITIAL QC ---

# --- Load Covariate Data and Apply Quality Control (QC) ---
cat("Loading covariate data and applying QC filters...\n")
covars <- fread(COVARIATE_FILE_PATH)

# Filter for participants who meet the quality criteria.
qc_passed_ids <- covars %>%
  filter(goodWear == 1, goodCal == 1, wearDays >= 4) %>%
  pull(id) # Use the complex ID from covariates file

cat(sprintf("Found %d participant records passing initial QC.\n", length(qc_passed_ids)))


# --- Load Temperature Data for QC-passed Participants ---
cat("Loading temperature and activity data from Parquet file...\n")
temp_dataset <- open_dataset(TEMPERATURE_DATA_PATH)

# Using a subset for demonstration purposes. Remove `[1:50]` to run on all QC-passed IDs.
raw_data <- temp_dataset %>%
  filter(id %in% qc_passed_ids[1:50]) %>%
  select(id, time, temp, enmoTrunc) %>%
  collect() %>%
  mutate(time = as_datetime(time))

cat(sprintf("Loaded data for %d participants.\n", length(unique(raw_data$id))))


# --- 3. DATA PREPROCESSING ---

cat("Starting data preprocessing: Z-score filtering and smoothing...\n")
processed_data <- raw_data %>%
  group_by(id) %>%
  arrange(time) %>%
  mutate(
    rolling_mean = RcppRoll::roll_mean(temp, n = 21, fill = NA, align = "center"),
    rolling_sd = RcppRoll::roll_sd(temp, n = 21, fill = NA, align = "center"),
    z_score = (temp - rolling_mean) / rolling_sd,
    temp_z_filtered = if_else(abs(z_score) > 3, NA_real_, temp)
  ) %>%
  mutate(
    temp_smoothed = RcppRoll::roll_median(temp_z_filtered, n = 7, fill = NA, align = "center")
  ) %>%
  mutate(hour = hour(time) + minute(time)/60) %>%
  ungroup()
cat("Preprocessing complete.\n")


# --- 4. DOWNSAMPLE DATA FOR EFFICIENCY ---
cat("Downsampling data to 5-minute epochs for modeling...\n")
downsampled_data <- processed_data %>%
  filter(!is.na(temp_smoothed)) %>%
  mutate(time_5min = floor_date(time, unit = "5 minutes")) %>%
  group_by(id, time_5min) %>%
  summarise(
    temp_final = mean(temp_smoothed, na.rm = TRUE),
    enmo_final = mean(enmoTrunc, na.rm = TRUE),
    hour = first(hour),
    .groups = 'drop'
  )
cat("Downsampling complete.\n")


# --- 5. LINK EXTERNAL AMBIENT TEMPERATURE ---

cat("Loading and linking external ambient temperature data...\n")

# Step 5a: Load external data files
location <- as.tibble(fread(LOCATION_FILE_PATH))
temps_at_location <- fread("https://github.com/gushamilton/ukb_temp_epochs/raw/refs/heads/main/data/ukb_temps_per_date.tsv.gz")

# Step 5b: Prepare data for joining
processed_data_with_eid <- downsampled_data %>%
  mutate(eid = as.integer(str_extract(id, "^[0-9]+")))

location_clean <- location %>%
  select(eid, assessment_centre = p54_i0)

temps_at_location_clean <- temps_at_location %>%
  mutate(
    assessment_centre = str_to_title(word(centre_name, 1)),
    date = as_date(date)
  ) %>%
  group_by(assessment_centre, date) %>%
  summarise(
    ambient_temp = mean(temp, na.rm = TRUE),
    .groups = 'drop'
  )

# Step 5c: Link datasets efficiently
data_with_location <- processed_data_with_eid %>%
  left_join(location_clean, by = "eid") %>%
  mutate(date = as_date(time_5min))

required_temps <- data_with_location %>%
  distinct(assessment_centre, date) %>%
  filter(!is.na(assessment_centre), !is.na(date))

temps_at_location_filtered <- temps_at_location_clean %>%
  semi_join(required_temps, by = c("assessment_centre", "date"))

cat(sprintf("Filtered external temperature data to %d unique centre-date rows.\n", nrow(temps_at_location_filtered)))

# Step 5d: Join ambient temp into the main downsampled data frame
# This becomes the primary data frame for all subsequent modeling.
downsampled_data <- data_with_location %>%
  left_join(temps_at_location_filtered, by = c("assessment_centre", "date"))

cat("Ambient temperature linking complete.\n")


# --- 6. VISUALIZE PREPROCESSING & MODELING STEPS ---

cat("Generating plots to visualize preprocessing and model fit for 4 random participants...\n")
set.seed(42)
n_vis_sample <- min(4, length(unique(processed_data$id)))
random_ids_vis <- sample(unique(processed_data$id), n_vis_sample)

for (i in seq_along(random_ids_vis)) {
  current_id <- random_ids_vis[i]
  cat("  - Plotting for participant:", current_id, "\n")
  
  participant_data_vis <- processed_data %>% filter(id == current_id)
  participant_model_data <- downsampled_data %>% filter(id == current_id)
  
  plot_a <- ggplot(participant_data_vis, aes(x = time, y = temp)) +
    geom_line(color = "gray50", alpha = 0.8) + labs(title = "A: Raw Temperature Data", y = "Temp (°C)", x = "") + theme_bw()
  plot_b <- ggplot(participant_data_vis, aes(x = time)) +
    geom_line(aes(y = temp_z_filtered), color = "dodgerblue") +
    geom_point(data = . %>% filter(is.na(temp_z_filtered)), aes(y = temp), color = "red", size = 0.5) +
    labs(title = "B: Z-Score Filtered Data", y = "Temp (°C)", x = "") + theme_bw()
  plot_c <- ggplot(participant_data_vis, aes(x = time, y = temp_smoothed)) +
    geom_line(color = "navy") + labs(title = "C: Smoothed Signal", y = "Temp (°C)", x = "Time") + theme_bw()
  
  lm_model <- lm(temp_final ~ cos(2*pi*hour/24) + sin(2*pi*hour/24), data = participant_model_data)
  participant_model_data$temp_predicted <- predict(lm_model)
  
  plot_d <- ggplot(participant_model_data, aes(x = time_5min)) +
    geom_point(aes(y = temp_final), color = "gray70", alpha = 0.5, size = 0.8) +
    geom_line(aes(y = temp_predicted), color = "darkorange", linewidth = 1.2) +
    labs(title = "D: Cosinor Model Fit on 5-min Data", y = "Temp (°C)", x = "Time") + theme_bw()
  
  final_vis_plot <- (plot_a | plot_b) / (plot_c | plot_d) +
    plot_annotation(title = paste("Preprocessing & Modeling Steps for Participant:", current_id))
  
  file_name <- file.path(OUTPUT_DIR, paste0("01_preprocessing_steps_participant_", i, ".pdf"))
  ggsave(file_name, final_vis_plot, width = 12, height = 9)
}


# --- 7. AGGREGATE 24-HOUR PATTERN ANALYSIS ---
# Note: This section uses 'processed_data' and is not affected by the ambient temp link.
cat("Analyzing and plotting aggregate 24-hour patterns...\n")

hourly_temp_summary <- processed_data %>%
  filter(!is.na(temp_smoothed)) %>% group_by(hour) %>%
  summarise(mean_val = mean(temp_smoothed, na.rm = TRUE), se_val = sd(temp_smoothed, na.rm = TRUE) / sqrt(n()), .groups = 'drop')
plot_24hr_temp <- ggplot(hourly_temp_summary, aes(x = hour, y = mean_val)) +
  geom_ribbon(aes(ymin = mean_val - se_val, ymax = mean_val + se_val), alpha = 0.2, fill = "blue") +
  geom_line(color = "blue", linewidth = 1) + scale_x_continuous(breaks = seq(0, 24, by = 2), limits = c(0, 24)) +
  labs(title = "Average Smoothed Temperature Across 24-Hour Cycle", x = "Hour of Day", y = "Mean Temperature (°C)") + theme_bw()
ggsave(file.path(OUTPUT_DIR, "02_aggregate_24hr_temperature.pdf"), plot_24hr_temp, width = 10, height = 6)


# --- 8. FOURIER & R-SQUARED ANALYSIS ---
cat("Performing Fourier and R-squared analysis...\n")
# Fourier Analysis
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
  unnest(temp_peak, names_sep = "_") %>% unnest(enmo_peak, names_sep = "_") %>%
  select(id, peak_freq_temp = temp_peak_peak_freq, peak_amp_temp = temp_peak_peak_amp, peak_freq_enmo = enmo_peak_peak_freq, peak_amp_enmo = enmo_peak_peak_amp)

# R-Squared Analysis
r2_by_id <- processed_data %>% filter(!is.na(temp_smoothed) & !is.na(enmoTrunc)) %>% group_by(id) %>%
  do(broom::glance(lm(temp_smoothed ~ enmoTrunc, data = .))) %>% ungroup() %>% select(id, r.squared)


# --- 9. PER-PARTICIPANT METRIC CALCULATION & CORRELATION ---

cat("Calculating summary metrics (Standard, Activity- & Ambient-Adjusted) for each participant...\n")

# Function to fit cosinor models and extract circadian parameters.
calculate_circadian_metrics <- function(df) {
  # --- Phenotypes from Temperature Distribution (as per research plan) ---
  p5 <- quantile(df$temp_final, 0.05, na.rm = TRUE)
  p95 <- quantile(df$temp_final, 0.95, na.rm = TRUE)
  
  # A simple, non-parametric measure of rhythm strength
  temp_range <- p95 - p5
  # Relative Amplitude, robust to mean shifts
  relative_amplitude <- (p95 - p5) / (p95 + p5)
  
  # --- Model 1: Unadjusted Cosinor ---
  unadj_model <- lm(temp_final ~ cos(2*pi*hour/24) + sin(2*pi*hour/24), data = df)
  coeffs_unadj <- coef(unadj_model)
  if (any(is.na(coeffs_unadj)) || length(coeffs_unadj) < 3) { mesor <- NA; amplitude <- NA; acrophase_h <- NA
  } else {
    mesor <- coeffs_unadj[1]; amplitude <- sqrt(coeffs_unadj[2]^2 + coeffs_unadj[3]^2); acrophase_rad <- atan2(coeffs_unadj[3], coeffs_unadj[2]); acrophase_h <- (-acrophase_rad * (24 / (2 * pi))) %% 24
  }
  
  # --- Model 2: Activity-Adjusted Cosinor ---
  df_act <- df %>% mutate(enmo_centered = enmo_final - mean(enmo_final, na.rm = TRUE))
  adj_model <- lm(temp_final ~ enmo_centered + cos(2*pi*hour/24) + sin(2*pi*hour/24), data = df_act)
  coeffs_adj <- coef(adj_model)
  if (any(is.na(coeffs_adj)) || length(coeffs_adj) < 4) { mesor_adj <- NA; amplitude_adj <- NA; acrophase_h_adj <- NA; activity_effect_slope <- NA
  } else {
    mesor_adj <- coeffs_adj[1]; amplitude_adj <- sqrt(coeffs_adj[3]^2 + coeffs_adj[4]^2); acrophase_rad_adj <- atan2(coeffs_adj[4], coeffs_adj[3]); acrophase_h_adj <- (-acrophase_rad_adj * (24 / (2 * pi))) %% 24
    # NEW PHENOTYPE: Quantifies the masking effect of activity
    activity_effect_slope <- coeffs_adj[2]
  }
  
  # --- Model 3: Ambient Temperature-Adjusted Cosinor ---
  df_amb <- df %>% mutate(ambient_temp_centered = ambient_temp - mean(ambient_temp, na.rm = TRUE))
  if (all(is.na(df_amb$ambient_temp_centered))) { # Check if ambient temp is available
    mesor_ambient <- NA; amplitude_ambient <- NA; acrophase_h_ambient <- NA; environmental_sensitivity <- NA; residual_variability <- NA
  } else {
    ambient_model <- lm(temp_final ~ ambient_temp_centered + cos(2*pi*hour/24) + sin(2*pi*hour/24), data = df_amb)
    coeffs_ambient <- coef(ambient_model)
    if (any(is.na(coeffs_ambient)) || length(coeffs_ambient) < 4) { mesor_ambient <- NA; amplitude_ambient <- NA; acrophase_h_ambient <- NA; environmental_sensitivity <- NA; residual_variability <- NA
    } else {
      mesor_ambient <- coeffs_ambient[1]; amplitude_ambient <- sqrt(coeffs_ambient[3]^2 + coeffs_ambient[4]^2); acrophase_rad_ambient <- atan2(coeffs_ambient[4], coeffs_ambient[3]); acrophase_h_ambient <- (-acrophase_rad_ambient * (24 / (2 * pi))) %% 24
      # NEW PHENOTYPE: Quantifies sensitivity to external temperature
      environmental_sensitivity <- coeffs_ambient[2]
      # NEW PHENOTYPE: Quantifies temperature "lability" or instability
      residual_variability <- sd(residuals(ambient_model), na.rm = TRUE)
    }
  }
  
  return(tibble(mesor, amplitude, acrophase_h, mesor_adj, amplitude_adj, acrophase_h_adj, mesor_ambient, amplitude_ambient, acrophase_h_ambient, temp_range, relative_amplitude, activity_effect_slope, environmental_sensitivity, residual_variability))
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
metrics_for_corr <- participant_summary_metrics %>% ungroup() %>%
  select_if(is.numeric) %>% select_if(~!all(is.na(.)))
corr_matrix <- cor(metrics_for_corr, use = "pairwise.complete.obs")

pdf(file.path(OUTPUT_DIR, "08_metrics_correlation_plot.pdf"), width = 10, height = 10)
corrplot(corr_matrix, method = "color", type = "upper", order = "hclust", addCoef.col = "black", tl.col = "black", tl.srt = 45, diag = FALSE, title = "Correlation Matrix of Participant Metrics", mar = c(0,0,2,0))
dev.off()
cat("Correlation plot saved.\n")


# --- 10. CREATE FINAL ARCHIVE ---
cat("Archiving all plots and data into a single file...\n")
tar(tarfile = paste0(OUTPUT_DIR, ".tar.gz"), files = OUTPUT_DIR, compression = "gzip")
cat(sprintf("Successfully created archive: %s.tar.gz\n", OUTPUT_DIR))
cat("\n--- SCRIPT FINISHED SUCCESSFULLY ---\n")
