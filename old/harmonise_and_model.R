
# COMBINED QC, VISUALIZATION, AND PHENOTYPE DERIVATION SCRIPT
#
# This script provides a comprehensive pipeline for processing wrist-worn temperature
# and activity data from the UK Biobank, focusing on the derivation of robust
# thermoregulatory phenotypes for Genome-Wide Association Studies (GWAS).
#
# Key functionalities include:
# 1.  **Data Loading & Preprocessing:** Reads raw epoch-level data, applies
#     robust outlier filtering (Z-score based), and smoothing (rolling median).
# 2.  **Data Downsampling:** Aggregates high-frequency 30-second data into
#     5-minute epochs for computational efficiency and noise reduction.
# 3.  **Circadian Modeling (Mixed-Effects Models - MLM):**
#     *   Fits a standard MLM to extract circadian parameters (Mesor, Amplitude, Acrophase)
#         from raw temperature data.
#     *   Calculates temperature residuals after accounting for linear activity effects,
#         and then fits a second MLM to these "activity-adjusted" residuals to derive
#         circadian parameters that are less confounded by physical activity.
# 4.  **Phenotype Extraction:** Derives key quantitative traits for each participant:
#     *   Mesor, Amplitude, Acrophase from both "circadian-only" and "activity-adjusted" MLMs.
#     *   Individual-specific "Activity Effect Slope" (how much temperature changes per
#         unit activity after circadian rhythm is removed).
#     *   FFT-derived amplitude and frequency from both raw and activity-adjusted temperature.
# 5.  **Comprehensive Visualization:** Generates a suite of plots to:
#     *   Illustrate the data cleaning and modeling steps for individual participants.
#     *   Show overall population-level patterns (e.g., average daily temperature/activity).
#     *   Compare and correlate the various derived phenotypes (FFT vs. MLM, circadian-only vs. activity-adjusted).
# 6.  **Correlation Matrix:** Provides a quantitative summary of relationships between all
#     derived phenotypes, aiding in the selection of independent traits for GWAS.
#
# NOTE: This script operates on a **SUBSET** of the full UK Biobank data for
# demonstration and and development purposes. The `file_path` variable must be updated
# to point to your specific data file.
#
# Output plots are saved to the 'plots/v2' subdirectory within the '02_harmonising_qc' folder.

# --- 1. Load Libraries ---
# Ensure all necessary R packages are installed and loaded.
# pacman is used for convenient loading/installation of multiple packages.
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, arrow, RcppRoll, patchwork, lubridate, lme4, broom, corrplot)

# --- 2. Load and Prepare Data ---
# Define the path to the input Parquet file containing temperature epoch data.
# This file is assumed to be generated from Part 2 of the overall research plan.
file_path <- "/Users/fh6520/Downloads/temperature_epochs_first100.parquet" # UPDATE THIS PATH TO YOUR DATA
if (!file.exists(file_path)) {
  stop("Data file not found. Please update the 'file_path' variable to your Parquet data.")
}
d <- read_parquet(file_path)

# --- 3. Full Pre-processing Pipeline ---
# This section implements the core data cleaning and smoothing steps.
# It transforms raw, noisy 30-second temperature data into a cleaner signal.
cat("Processing all loaded participants (outlier filtering and smoothing)...
")

processed_data <- d %>%
  # Convert 'time' column to a proper datetime object for time-based operations.
  mutate(time = as_datetime(time)) %>%
  # Group data by participant ID to apply rolling functions independently for each individual.
  group_by(id) %>%
  # Step A: Z-score Filtering to remove local, non-physiological spikes.
  # A rolling Z-score is robust to local signal characteristics and avoids bias from overall trends.
  mutate(
    rolling_mean = RcppRoll::roll_mean(temp, n = 121, fill = NA, align = "center"), # 121 points = 60.5 minutes window
    rolling_sd = RcppRoll::roll_sd(temp, n = 121, fill = NA, align = "center"),
    z_score = (temp - rolling_mean) / rolling_sd,
    temp_z_filtered = if_else(abs(z_score) > 3, NA_real_, temp) # Replace outliers (Z > 3) with NA
  ) %>%
  # Step B: Smoothing the Z-score filtered data to remove high-frequency noise.
  # A rolling median is used for its robustness to any remaining outliers.
  mutate(
    temp_smoothed = RcppRoll::roll_median(
      temp_z_filtered,
      n = 21,          # 21 points = 10.5 minutes window (21 * 30s/point)
      fill = NA,
      align = "center"
    )
  ) %>%
  # Calculate residuals of RAW temp ~ enmoTrunc for each participant at 30-second level.
  # This represents temperature variation NOT linearly explained by activity at the unsmoothed level.
  mutate(
    temp_residuals_from_activity_30s = {
      if (sum(!is.na(temp) & !is.na(enmoTrunc)) >= 2 && var(enmoTrunc, na.rm = TRUE) > 1e-9) {
        lm_fit <- lm(temp ~ enmoTrunc, data = cur_data()) # Use cur_data() for current group's data
        residuals(lm_fit)
      } else {
        NA_real_ # Assign NA if model cannot be fitted
      }
    }
  ) %>%
  ungroup() # Ungroup the data after per-individual calculations are complete.

cat("Processing complete.
")

# --- 4. Downsample Data to 5-Minute Epochs ---
# This step reduces data volume and computational load for subsequent modeling,
# while retaining the essential circadian signal.
cat("Downsampling data to 5-minute epochs...
")

downsampled_data <- processed_data %>%
  # Remove rows where smoothed temperature is NA before aggregation.
  filter(!is.na(temp_smoothed)) %>%
  group_by(id) %>%
  # Create a new time column representing the start of each 5-minute interval.
  mutate(time_5min = floor_date(time, unit = "5 minutes")) %>%
  # Group by the new 5-minute epoch to aggregate points within each interval.
  group_by(id, time_5min) %>%
  # Calculate the mean of the smoothed temperature and ENMO for each 5-minute window.
  summarise(
    temp_final = mean(temp_smoothed, na.rm = TRUE),
    enmo_final = mean(enmoTrunc, na.rm = TRUE),
    # Calculate the mean of the 30-second activity-adjusted residuals for 5-minute epoch.
    temp_adjusted_by_activity_5min = mean(temp_residuals_from_activity_30s, na.rm = TRUE),
    .groups = 'drop' # Drop grouping structure after summarising.
  )

cat("Downsampling complete.
")

# --- 5. Modeling ---
# This section fits mixed-effects models to derive circadian parameters.
cat("Fitting models...
")

# Prepare data for MLM by creating harmonic time variables (sine and cosine)
# which capture the 24-hour periodicity of circadian rhythms.
model_data_all <- downsampled_data %>%
  mutate(
    hour_of_day = hour(time_5min) + minute(time_5min) / 60,
    time_radian = 2 * pi * hour_of_day / 24 # Convert hour of day to radians (0 to 2*pi)
  )

# Model 1: Standard MLM (Circadian Only)
# This model estimates the population-level circadian rhythm and individual deviations.
# Fixed effects: Sine and cosine terms for 24-hour rhythm.
# Random effects: Individual-specific intercepts, sine, and cosine slopes.
mlm_model <- lmer(temp_final ~ sin(time_radian) + cos(time_radian) +
                    (1 + sin(time_radian) + cos(time_radian) | id),
                  data = model_data_all)

# Model 2: MLM on Activity-Adjusted Temperature (temp_adjusted_by_activity_5min)
# This model aims to capture the "purer" circadian rhythm after removing the linear
# confounding effect of physical activity at the 30-second level, then aggregated.
mlm_model_active <- lmer(temp_adjusted_by_activity_5min ~ sin(time_radian) + cos(time_radian) +
                           (1 + sin(time_radian) + cos(time_radian) | id),
                         data = model_data_all)

cat("Modeling complete.
")

# --- 6. Visualization ---
# Generates various plots to visualize the data processing, modeling, and derived phenotypes.
cat("Generating visualizations...
")

# Define the output directory for all plots.
output_dir <- "/Users/fh6520/R/temp_ukb/02_harmonising_qc/plots/v2"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE) # Create directory if it doesn't exist

set.seed(42) # Set seed for reproducibility of random participant selection
random_ids <- sample(unique(processed_data$id), 4) # Select 4 random participants for detailed plots

for (i in 1:length(random_ids)) {
  current_id <- random_ids[i]
  cat(paste0("
--- Generating detailed plots for participant: ", current_id, " ---
"))

  # Filter data for the current participant for plotting.
  participant_data_vis <- processed_data %>% filter(id == current_id)
  model_data_single <- model_data_all %>% filter(id == current_id)

  # --- Predictions for Individual Participant Plots ---
  # Predict temperature from the MLM (circadian only).
  model_data_single$temp_predicted_mlm_circadian <- predict(mlm_model, newdata = model_data_single)

  # Calculate residuals from the MLM (circadian only) for the current participant.
  # These are the deviations from the circadian rhythm.
  model_data_single$mlm_residuals_circadian <- model_data_single$temp_final - model_data_single$temp_predicted_mlm_circadian

  # Fit LM of residuals vs activity for the current participant to show activity effect.
  # This models how activity explains the *non-circadian* variation.
  if (nrow(model_data_single) > 1 && !all(is.na(model_data_single$mlm_residuals_circadian)) && !all(is.na(model_data_single$enmo_final))) {
    lm_residual_activity <- lm(mlm_residuals_circadian ~ enmo_final, data = model_data_single)
    # Predict the activity effect on residuals and add back to MLM (circadian only) prediction.
    # This shows the circadian rhythm PLUS the activity-explained deviation.
    model_data_single$temp_predicted_mlm_plus_activity_effect <- model_data_single$temp_predicted_mlm_circadian + predict(lm_residual_activity, newdata = model_data_single)
  } else {
    model_data_single$temp_predicted_mlm_plus_activity_effect <- NA_real_ # Assign NA if model cannot be fitted
  }

  # Predict temperature from the MLM on activity-adjusted temperature.
  model_data_single$temp_predicted_mlm_active <- predict(mlm_model_active, newdata = model_data_single)

  # --- Generate 5-Panel Visualization for Each Participant ---
  # Panel A: Raw Data
  plot_a <- ggplot(participant_data_vis, aes(x = time, y = temp)) +
    geom_line(color = "gray50", alpha = 0.8) +
    labs(title = "A: Raw Data", y = "Temperature (°C)", x = "Time") +
    theme_bw()

  # Panel B: Z-Score Filtered Data
  plot_b <- ggplot(participant_data_vis, aes(x = time)) +
    geom_line(aes(y = temp_z_filtered), color = "dodgerblue") +
    geom_point(data = . %>% filter(is.na(temp_z_filtered)), aes(y = temp), color = "red", size = 0.5) +
    labs(title = "B: Z-Score Filtered Data", y = "Temperature (°C)", x = "Time") +
    theme_bw()

  # Panel C: Smoothed Temperature Signal
  plot_c <- ggplot(participant_data_vis, aes(x = time, y = temp_smoothed)) +
    geom_line(color = "navy") +
    labs(title = "C: Smoothed Temperature Signal", y = "Temperature (°C)", x = "Time") +
    theme_bw()

  # Panel D: Model Comparison & Activity Effect
  # Compares the circadian-only MLM with the circadian MLM plus the activity effect on residuals.
  plot_d <- ggplot(model_data_single, aes(x = time_5min)) +
    geom_point(aes(y = temp_final), color = "gray70", alpha = 0.5, size = 0.8) +
    geom_line(aes(y = temp_predicted_mlm_circadian, color = "MLM (Circadian Only)"), linewidth = 1) +
    geom_line(aes(y = temp_predicted_mlm_plus_activity_effect, color = "MLM (Circadian) + Activity Effect"), linewidth = 1) +
    scale_color_manual(values = c(
      "MLM (Circadian Only)" = "purple",
      "MLM (Circadian) + Activity Effect" = "darkred"
    )) +
    labs(
      title = "D: Model Comparison & Activity Effect",
      subtitle = "Circadian-Only vs. Circadian + Activity Effect on Residuals",
      y = "Temperature (°C)", x = "Time", color = "Model Type"
    ) +
    theme_bw() + theme(legend.position = "bottom")

  # Panel E: MLM on Activity-Adjusted Temperature
  # Shows the circadian rhythm after the linear effect of activity has been removed.
  plot_e <- ggplot(model_data_single, aes(x = time_5min)) +
    geom_point(aes(y = temp_adjusted_by_activity_5min), color = "gray70", alpha = 0.5, size = 0.8) +
    geom_line(aes(y = temp_predicted_mlm_active, color = "MLM on Activity-Adjusted Temp"), linewidth = 1) +
    scale_color_manual(values = c("MLM on Activity-Adjusted Temp" = "blue")) +
    labs(
      title = "E: MLM on Activity-Adjusted Temperature",
      subtitle = "Circadian Rhythm After Linear Activity Effect Removed (Unsmoothed Adjustment)",
      y = "Activity-Adjusted Temp (°C)", x = "Time", color = "Model Type"
    ) +
    theme_bw() + theme(legend.position = "bottom")

  # Combine all 5 plots into a single multi-panel PDF.
  final_plot <- (plot_a | plot_b) / (plot_c | plot_d) / plot_e +
    plot_annotation(title = paste("Preprocessing & Modeling Steps for Participant:", current_id))

  # Construct the full file path for saving the plot.
  full_file_path <- paste0(output_dir, "/", paste0("participant_", current_id, "_visualization.pdf"))
  ggsave(filename = full_file_path, plot = final_plot, width = 12, height = 12)

  cat("Saved plot to", full_file_path, "
")
}

# --- 7. Exploratory Data Analysis ---
# This section performs various exploratory analyses and generates plots
# to understand the distributions and relationships of key parameters.
cat("Performing exploratory data analysis...
")

# FFT Analysis Function
# This function calculates the amplitude and frequency of the first (dominant)
# peak in the Fourier spectrum of a given time series.
get_first_peak <- function(df, value_col_name, time_col_name) {
  # Arrange data by the specified time column.
  df <- df %>% arrange(!!sym(time_col_name))
  # Extract the time series data.
  value_data <- df[[value_col_name]]
  # Adjust data by subtracting its mean for FFT.
  value_adjust <- value_data - mean(value_data, na.rm = TRUE)
  N <- length(value_adjust)

  # Return NA if there are not enough data points for a meaningful FFT.
  if (N < 10) return(data.frame(id = unique(df$id), first_peak_freq = NA, first_peak_amp = NA)) # Return NA if not enough data

  # Perform Fast Fourier Transform.
  ft <- fft(value_adjust)
  # Calculate frequencies corresponding to each FFT component.
  freq <- (0:(N-1)) / N
  # Calculate amplitude of each FFT component, normalized to be comparable to original signal amplitude.
  amplitude <- Mod(ft) / N * 2
  # Consider only the first half of the spectrum (excluding DC component) as it's symmetric.
  half <- 2:floor(N/2)
  amp_half <- amplitude[half]
  freq_half <- freq[half]

  # Find the first local maximum (peak) in the amplitude spectrum.
  peak_idx <- which(diff(sign(diff(amp_half))) == -2) + 1
  if (length(peak_idx) == 0) {
    # Fallback: if no local maximum, take the global maximum (excluding DC).
    max_idx <- which.max(amp_half)
    first_peak <- max_idx
  } else {
    first_peak <- peak_idx[1] # Take the first identified peak.
  }
  first_peak_freq <- freq_half[first_peak]
  first_peak_amp <- amp_half[first_peak]

  # Return a data frame with the participant ID, peak frequency, and peak amplitude.
  data.frame(id = unique(df$id), first_peak_freq = first_peak_freq, first_peak_amp = first_peak_amp)
}

# Apply FFT to raw temperature data (from processed_data, which is 30-second data).
fourier_peaks_raw <- processed_data %>%
  group_by(id) %>%
  group_split() %>%
  map_df(~get_first_peak(., "temp", "time")) # Use "temp" and "time" columns
names(fourier_peaks_raw) <- c("id", "first_peak_freq_raw", "first_peak_amp_raw")

# Apply FFT to activity-adjusted temperature data (from processed_data, which is 30-second data).
fourier_peaks_active <- processed_data %>%
  group_by(id) %>%
  group_split() %>%
  map_df(~get_first_peak(., "temp_residuals_from_activity_30s", "time")) # Use "temp_residuals_from_activity_30s" and "time"
names(fourier_peaks_active) <- c("id", "first_peak_freq_active", "first_peak_amp_active")

# Plot distributions of FFT frequencies and amplitudes for both raw and activity-adjusted temperature.
p_fft_freq_raw <- ggplot(fourier_peaks_raw, aes(x = first_peak_freq_raw)) +
  geom_histogram(bins = 50, fill = "skyblue", color = "black") +
  labs(title = "FFT Frequency (Raw Temp)", x = "Frequency (cycles/sample)", y = "Count") +
  theme_bw()

p_fft_amp_raw <- ggplot(fourier_peaks_raw, aes(x = first_peak_amp_raw)) +
  geom_histogram(bins = 50, fill = "orange", color = "black") +
  labs(title = "FFT Amplitude (Raw Temp)", x = "Amplitude", y = "Count") +
  theme_bw()

p_fft_freq_active <- ggplot(fourier_peaks_active, aes(x = first_peak_freq_active)) +
  geom_histogram(bins = 50, fill = "skyblue", color = "black") +
  labs(title = "FFT Frequency (Activity-Adjusted Temp)", x = "Frequency (cycles/sample)", y = "Count") +
  theme_bw()

p_fft_amp_active <- ggplot(fourier_peaks_active, aes(x = first_peak_amp_active)) +
  geom_histogram(bins = 50, fill = "orange", color = "black") +
  labs(title = "FFT Amplitude (Activity-Adjusted Temp)", x = "Amplitude", y = "Count") +
  theme_bw()

# Save the combined FFT distribution plots.
ggsave(file.path(output_dir, "fft_distributions_raw_vs_active.pdf"),
       (p_fft_freq_raw | p_fft_amp_raw) / (p_fft_freq_active | p_fft_amp_active),
       width = 12, height = 10)

# --- Overall and Individual-level Temperature-Activity Correlations ---

# Add hour of day (0-23) for each sample for hourly summaries.
processed_data_hourly <- processed_data %>%
  mutate(hour = hour(time) + minute(time)/60)

# Compute mean and standard error for each hour across all participants.
hourly_summary <- processed_data_hourly %>%
  group_by(hour) %>%
  summarise(
    mean_temp = mean(temp_smoothed, na.rm = TRUE),
    se_temp = sd(temp_smoothed, na.rm = TRUE) / sqrt(sum(!is.na(temp_smoothed))),
    mean_enmo = mean(enmoTrunc, na.rm = TRUE),
    .groups = 'drop'
  )

# Plot mean temperature and inverse activity over 24 hours.
# Rescale inverse activity to match the temperature axis for visual comparison.
activity_min <- min(hourly_summary$mean_enmo, na.rm = TRUE)
activity_max <- max(hourly_summary$mean_enmo, na.rm = TRUE)
temp_min <- min(hourly_summary$mean_temp, na.rm = TRUE)
temp_max <- max(hourly_summary$mean_temp, na.rm = TRUE)

hourly_summary <- hourly_summary %>%
  mutate(
    mean_enmo_inv = activity_max - mean_enmo, # Invert activity (high activity = low value)
    mean_enmo_inv_rescaled = (mean_enmo_inv - min(mean_enmo_inv, na.rm = TRUE)) /
                             (max(mean_enmo_inv, na.rm = TRUE) - min(mean_enmo_inv, na.rm = TRUE)) *
                             (temp_max - temp_min) + temp_min # Rescale to temperature range
  )

p_overall_temp_activity <- ggplot() +
  geom_line(data = hourly_summary, aes(x = hour, y = mean_temp, color = "Mean Temperature"), linewidth = 1) +
  geom_ribbon(data = hourly_summary, aes(x = hour, ymin = mean_temp - se_temp, ymax = mean_temp + se_temp),
              alpha = 0.2, fill = "blue") +
  geom_line(data = hourly_summary, aes(x = hour, y = mean_enmo_inv_rescaled, color = "Inverse Activity (Rescaled)"), linewidth = 1, linetype = "dashed") +
  scale_x_continuous(breaks = seq(0, 24, by = 2), limits = c(0, 24)) +
  scale_color_manual(values = c("Mean Temperature" = "blue", "Inverse Activity (Rescaled)" = "red")) +
  labs(
    title = "Average Temperature and Inverse Activity Across 24-Hour Clock",
    x = "Hour of Day",
    y = "Temperature / Rescaled Activity",
    color = "" # Remove legend title
  ) +
  theme_bw() + theme(legend.position = "bottom")

ggsave(file.path(output_dir, "overall_temp_activity_aligned.pdf"), p_overall_temp_activity, width = 12, height = 6)

# Individual-level correlations between hourly temperature and activity.
set.seed(123) # Set seed for reproducibility
# Sample up to 100 unique participant IDs for correlation calculation.
sampled_ids_corr <- sample(unique(processed_data_hourly$id), min(100, length(unique(processed_data_hourly$id))))

indiv_corrs <- sapply(sampled_ids_corr, function(pid) {
  df <- processed_data_hourly %>%
    filter(id == pid) %>%
    group_by(hour) %>%
    summarise(
      mean_temp = mean(temp_smoothed, na.rm = TRUE),
      mean_enmo = mean(enmoTrunc, na.rm = TRUE),
      .groups = 'drop'
    )
  # Compute correlation only if there are at least 5 hourly data points.
  if (nrow(df) >= 5) {
    cor(df$mean_temp, df$mean_enmo, use = "complete.obs")
  } else {
    NA_real_ # Assign NA if not enough data
  }
})

indiv_corrs_df <- data.frame(correlation = indiv_corrs)
p_indiv_corr_hist <- ggplot(indiv_corrs_df, aes(x = correlation)) +
  geom_histogram(binwidth = 0.05, fill = "skyblue", color = "black") +
  labs(
    title = "Distribution of Per-Individual Hourly Temp-Activity Correlations",
    x = "Correlation (Temperature vs. Activity)",
    y = "Count"
  ) +
  theme_bw()

ggsave(file.path(output_dir, "individual_temp_activity_correlations.pdf"), p_indiv_corr_hist, width = 8, height = 5)


# --- Extract MLM Parameters (Mesor, Amplitude, Acrophase) ---

# Function to extract Mesor, Amplitude, and Acrophase from lmer random effects.
# These parameters are derived from the individual-specific random effects
# combined with the population-level fixed effects.
extract_mlm_params <- function(mlm_object, model_type) {
  # Extract random effects for each participant.
  ranef_df <- as.data.frame(ranef(mlm_object)$id)

  # Assign column names for the random effects (intercept, sine, cosine).
  names(ranef_df) <- c("intercept", "sin_coef", "cos_coef")

  # Add fixed effects (population-level averages) to random effects for each participant.
  fixed_effects <- fixef(mlm_object)

  ranef_df <- ranef_df %>%
    rownames_to_column(var = "id") %>%
    mutate(
      intercept_total = intercept + fixed_effects["(Intercept)"],
      sin_coef_total = sin_coef + fixed_effects["sin(time_radian)"],
      cos_coef_total = cos_coef + fixed_effects["cos(time_radian)"]
    ) %>%
    rowwise() %>% # Perform calculations row by row for each participant.
    mutate(
      Mesor = intercept_total, # Mesor is the rhythm-adjusted mean.
      Amplitude = sqrt(sin_coef_total^2 + cos_coef_total^2), # Amplitude of the rhythm.
      Acrophase_rad = atan2(-sin_coef_total, cos_coef_total), # Acrophase in radians.
      Acrophase_hr = (Acrophase_rad %% (2 * pi)) / (2 * pi) * 24 # Convert Acrophase to 0-24h scale.
    ) %>%
    ungroup() %>%
    select(id, Mesor, Amplitude, Acrophase_hr) # Select only the derived parameters.

  # Rename columns to indicate the model type (e.g., _circadian, _active).
  names(ranef_df) <- paste0(names(ranef_df), "_", model_type)
  names(ranef_df)[1] <- "id" # Keep 'id' column name consistent.
  return(ranef_df)
}

# Extract parameters from the "circadian only" MLM.
mlm_params_circadian <- extract_mlm_params(mlm_model, "circadian")
# Extract parameters from the "activity-adjusted" MLM.
mlm_params_active <- extract_mlm_params(mlm_model_active, "active")

# Extract individual-specific activity slopes from the residual model.
# This quantifies how much an individual's temperature changes per unit activity
# after accounting for their circadian rhythm.
individual_activity_slopes <- model_data_all %>%
  group_by(id) %>%
  do({
    current_participant_data <- .
    # Filter out NA values for robust regression.
    df_filtered <- current_participant_data %>%
      filter(!is.na(temp_final), !is.na(enmo_final))

    activity_slope_val <- NA_real_ # Initialize slope to NA

    # Check if there's enough data for this participant after filtering (at least 2 points for lm).
    if (nrow(df_filtered) >= 2) {
      # Calculate residuals from the MLM (circadian only) for this participant.
      temp_predicted_mlm_single <- predict(mlm_model, newdata = df_filtered)
      mlm_residuals_single <- df_filtered$temp_final - temp_predicted_mlm_single

      # Check for sufficient variance in both residuals and ENMO for lm.
      if (var(mlm_residuals_single, na.rm = TRUE) > 1e-9 && var(df_filtered$enmo_final, na.rm = TRUE) > 1e-9) {
        lm_res_act <- lm(mlm_residuals_single ~ df_filtered$enmo_final)
        # Check if the coefficient for 'df_filtered$enmo_final' exists before extracting.
        if ("df_filtered$enmo_final" %in% names(coef(lm_res_act))) {
          activity_slope_val <- coef(lm_res_act)[["df_filtered$enmo_final"]]
        }
      }
    }
    data.frame(activity_slope = activity_slope_val)
  }) %>%
  ungroup()

names(individual_activity_slopes)[2] <- "Activity_Effect_Slope" # Rename for clarity

# Merge all derived parameters into a single comprehensive dataframe.
all_params <- fourier_peaks_raw %>%
  left_join(fourier_peaks_active, by = "id") %>%
  left_join(mlm_params_circadian, by = "id") %>%
  left_join(mlm_params_active, by = "id") %>%
  left_join(individual_activity_slopes, by = "id")

# --- Comparisons between FFT and MLM Parameters ---
# These plots visually compare different methods of deriving similar phenotypes.

# FFT Amplitude (Raw) vs. MLM Amplitude (Circadian Only)
p_fft_amp_vs_mlm_amp_circadian <- ggplot(all_params, aes(x = first_peak_amp_raw, y = Amplitude_circadian)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "FFT Amplitude (Raw) vs. MLM Amplitude (Circadian Only)",
       x = "FFT Amplitude (Raw)", y = "MLM Amplitude (Circadian Only)") +
  theme_bw()

# FFT Frequency (Raw) vs. MLM Acrophase (Circadian Only)
p_fft_freq_vs_mlm_acro_circadian <- ggplot(all_params, aes(x = first_peak_freq_raw, y = Acrophase_hr_circadian)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "FFT Frequency (Raw) vs. MLM Acrophase (Circadian Only)",
       x = "FFT Frequency", y = "MLM Acrophase (Hours)") +
  theme_bw()

# MLM Amplitude: Circadian Only vs. Activity-Adjusted
p_mlm_amp_comp_active <- ggplot(all_params, aes(x = Amplitude_circadian, y = Amplitude_active)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "MLM Amplitude: Circadian Only vs. Activity-Adjusted",
       x = "MLM Amplitude (Circadian Only)", y = "MLM Amplitude (Activity-Adjusted)") +
  theme_bw()

# MLM Mesor: Circadian Only vs. Activity-Adjusted
p_mlm_mesor_comp_active <- ggplot(all_params, aes(x = Mesor_circadian, y = Mesor_active)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "MLM Mesor: Circadian Only vs. Activity-Adjusted",
       x = "MLM Mesor (Circadian Only)", y = "MLM Mesor (Activity-Adjusted)") +
  theme_bw()

# MLM Acrophase: Circadian Only vs. Activity-Adjusted
p_mlm_acro_comp_active <- ggplot(all_params, aes(x = Acrophase_hr_circadian, y = Acrophase_hr_active)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "MLM Acrophase: Circadian Only vs. Activity-Adjusted",
       x = "MLM Acrophase (Circadian Only)", y = "MLM Acrophase (Activity-Adjusted)") +
  theme_bw()

# Save the combined parameter comparison plots.
ggsave(file.path(output_dir, "parameter_comparisons.pdf"),
       (p_fft_amp_vs_mlm_amp_circadian | p_fft_freq_vs_mlm_acro_circadian) /
       (p_mlm_amp_comp_active | p_mlm_mesor_comp_active | p_mlm_acro_comp_active),
       width = 15, height = 10)

# --- Correlation Matrix of All Key Parameters ---
# Visualizes the inter-correlations between all derived phenotypes.
# This helps identify independent traits for GWAS.

# Select only the numeric columns for the correlation matrix and remove any NAs.
corr_data <- all_params %>%
  select(-id) %>% # Exclude 'id' column
  drop_na() # Remove rows with any NA values (participants with incomplete data)

# Compute the Pearson correlation matrix.
correlation_matrix <- cor(corr_data)

# Plot the correlation matrix using corrplot.
# Increased size for better readability of labels and coefficients.
png(file.path(output_dir, "parameter_correlation_matrix.png"), width = 1000, height = 1000)
corrplot::corrplot(correlation_matrix, method = "circle", type = "upper",
                   tl.col = "black", tl.srt = 45, # Text label color and string rotation
                   addCoef.col = "black", number.cex = 0.7, # Add correlation coefficients
                   col = colorRampPalette(c("blue", "white", "red"))(200)) # Color palette
dev.off() # Close the PNG device

# --- Additional Plots (Appendix) ---
# These plots provide further exploratory insights and are included as an appendix.

# Plot 4 random participants: individual inverse temp vs activity.
set.seed(123) # For reproducibility
random_ids_inv <- sample(unique(processed_data$id), 4)

for (pid in random_ids_inv) {
  df <- processed_data %>% filter(id == pid)
  df <- df %>%
    mutate(
      inv_temp = -scale(temp_smoothed), # Inverse temperature (z-scored for comparability)
      enmo_z = scale(enmoTrunc) # Z-scored ENMO for comparability
    )
  p <- ggplot(df, aes(x = time)) +
    geom_line(aes(y = inv_temp, color = "Inverse Temp (z)"), linewidth = 1) +
    geom_line(aes(y = enmo_z, color = "Activity (ENMO, z)"), linewidth = 1, alpha = 0.7) +
    scale_color_manual(values = c("Inverse Temp (z)" = "blue", "Activity (ENMO, z)" = "darkgreen")) +
    labs(
      title = paste("Participant", pid, "- Inverse Temp vs Activity"),
      x = "Time",
      y = "z-score",
      color = ""
    ) +
    theme_bw() + theme(legend.position = "top")
  ggsave(file.path(output_dir, paste0("participant_", pid, "_invtemp_vs_activity.pdf")), p, width = 10, height = 4)
}

# Overlay all individuals, aligning their time so that each starts at 0.
# This visualizes the overall population trend of temperature over time.
processed_data_aligned <- processed_data %>%
  group_by(id) %>%
  arrange(time) %>%
  mutate(
    time_since_start = as.numeric(difftime(time, min(time), units = "hours")),
    temp_z = as.numeric(scale(temp_smoothed))
  ) %>%
  ungroup() %>%
  filter(!is.na(temp_smoothed), abs(temp_z) <= 5) # Filter out extreme Z-scores for cleaner visualization

p_all_temp_aligned <- ggplot(processed_data_aligned, aes(x = time_since_start, y = temp_smoothed, group = id)) +
  geom_line(alpha = 0.3, color = "blue") +
  labs(
    title = "All Participants: Temperature (Aligned at Start, Extreme Z > 5 Removed)",
    x = "Time Since Start (hours)",
    y = "Temperature (°C)"
  ) +
  theme_bw()

ggsave(file.path(output_dir, "all_participants_temperature_aligned.pdf"), p_all_temp_aligned, width = 12, height = 6)

# Calculate the average residual temperature after removing activity effects over 5 days.
# This shows the remaining temperature rhythm after activity is accounted for.
processed_data_resid_for_plot <- processed_data %>%
  filter(!is.na(temp_smoothed), !is.na(enmoTrunc)) %>%
  group_by(id) %>%
  mutate(
    temp_resid = {
      if (sum(!is.na(temp_smoothed) & !is.na(enmoTrunc)) >= 2 && var(enmoTrunc, na.rm = TRUE) > 1e-9) {
        lm_fit <- lm(temp_smoothed ~ enmoTrunc, data = cur_data())
        residuals(lm_fit)
      } else {
        NA_real_
      }
    }
  ) %>%
  ungroup()

processed_data_resid_for_plot <- processed_data_resid_for_plot %>%
  group_by(id) %>%
  arrange(time) %>%
  mutate(
    time_since_start = as.numeric(difftime(time, min(time), units = "hours"))
  ) %>%
  ungroup()

processed_data_resid_5d <- processed_data_resid_for_plot %>%
  filter(time_since_start <= 120) # Limit to first 5 days (120 hours)

processed_data_resid_5d <- processed_data_resid_5d %>%
  mutate(
    time_bin = floor(time_since_start * 2) / 2  # Bin time into 0.5 hour intervals
  )

summary_resid <- processed_data_resid_5d %>%
  group_by(time_bin) %>%
  summarise(
    mean_resid = mean(temp_resid, na.rm = TRUE),
    sem_resid = sd(temp_resid, na.rm = TRUE) / sqrt(n()), # Standard Error of the Mean
    .groups = 'drop'
  )

p_avg_residual_temp <- ggplot(summary_resid, aes(x = time_bin, y = mean_resid)) +
  geom_line(color = "purple", linewidth = 1) +
  geom_ribbon(aes(ymin = mean_resid - sem_resid, ymax = mean_resid + sem_resid), alpha = 0.2, fill = "purple") +
  labs(
    title = "Average Residual Temperature (After Activity Removed) Over 5 Days",
    x = "Time Since Start (hours)",
    y = "Residual Temperature (°C)"
  ) +
  theme_bw()

ggsave(file.path(output_dir, "average_residual_temp_over_5days.pdf"), p_avg_residual_temp, width = 12, height = 5)

# Plot 4 random participants' residual temperature over time.
set.seed(123) # For reproducibility
random_ids_resid <- sample(unique(processed_data_resid_5d$id), 4)

for (pid in random_ids_resid) {
  pdata <- processed_data_resid_5d %>% filter(id == pid)
  pdata_summary <- pdata %>%
    group_by(time_bin) %>%
    summarise(
      mean_resid = mean(temp_resid, na.rm = TRUE),
      sem_resid = sd(temp_resid, na.rm = TRUE) / sqrt(n()),
      .groups = 'drop'
    )
  p <- ggplot(pdata_summary, aes(x = time_bin, y = mean_resid)) +
    geom_line(color = "darkgreen", linewidth = 1) +
    geom_ribbon(aes(ymin = mean_resid - sem_resid, ymax = mean_resid + sem_resid), alpha = 0.2, fill = "darkgreen") +
    labs(
      title = paste("Residual Temperature for Participant", pid),
      x = "Time Since Start (hours)",
      y = "Residual Temperature (°C)"
    ) +
    theme_bw()
  ggsave(file.path(output_dir, paste0("residual_temp_participant_", pid, ".pdf")), p, width = 10, height = 4)
}

cat("Analysis complete. Plots saved to:", output_dir, "
")
