# Title: Visualisation Test Script for Batch Harmonisation and Modeling (Complete)
# Author: Gemini
# Date: 2025-07-08
# Description: This script processes ONLY THE FIRST BATCH of participants from a large
#              Parquet file. It then runs a comprehensive suite of visualisations,
#              including Intra-class Correlation (ICC) analysis, to allow for
#              detailed inspection of the data cleaning and modeling pipeline.

# --- 1. Setup & Configuration ---
library(tidyverse)
library(arrow)
library(data.table)
library(tidyr)
library(purrr)
library(lubridate)
library(zoo)
library(RcppRoll)
library(patchwork)
library(broom)
library(lme4)
library(corrplot)
library(ICC)

# --- Define Global Variables & Paths ---
COVARIATE_FILE_PATH <- "combined_covariates.csv.gz"
TEMPERATURE_DATA_PATH <- "combined_temperature_data.parquet"
OUTPUT_DIR <- "results_test_plots"

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
BATCH_SIZE_PARTICIPANTS <- 50

# --- 2. Load Covariates and Select Participants ---
print("--- Starting: Loading and filtering covariate data ---")

covars <- fread(COVARIATE_FILE_PATH)
participant_counts <- covars %>% count(participant_id)
repeat_participant_ids <- participant_counts %>% filter(n > 1) %>% pull(participant_id)
repeat_covars <- covars %>%
  filter(participant_id %in% repeat_participant_ids)

print(paste0("Found ", length(repeat_participant_ids), " participants with repeat measurements."))

qc_passed_repeats <- repeat_covars %>%
  filter(goodWear == 1, goodCal == 1, wearDays >= 4)

get_season <- function(date) {
  month <- month(date)
  case_when(
    month %in% c(12, 1, 2) ~ "Winter",
    month %in% c(3, 4, 5) ~ "Spring",
    month %in% c(6, 7, 8) ~ "Summer",
    month %in% c(9, 10, 11) ~ "Autumn",
    TRUE ~ NA_character_
  )
}

qc_passed_repeats <- qc_passed_repeats %>%
  mutate(
    startTime = as_datetime(startTime),
    season = get_season(startTime)
  )

qc_ids <- unique(qc_passed_repeats$id)
print(paste0("After QC (wearDays >= 4), analyzing ", length(qc_ids), " individual measurements."))

# --- 3. Open Parquet Dataset ---
print("--- Starting: Opening Parquet dataset ---")
temp_dataset <- open_dataset(TEMPERATURE_DATA_PATH)
print("Parquet dataset opened.")


# --- 4. Modular Data Processing & Modeling Functions ---

load_raw_data_for_ids <- function(ids_to_load, temp_dataset_ref) {
  temp_dataset_ref %>%
    filter(id %in% ids_to_load) %>%
    collect() %>%
    rename(event_timestamp = time) %>%
    mutate(event_timestamp = as_datetime(event_timestamp))
}

preprocess_and_downsample <- function(raw_data_df) {
  raw_data_df %>%
    rename(time = event_timestamp) %>%
    group_by(id) %>%
    arrange(time) %>%
    mutate(
      rolling_mean = RcppRoll::roll_mean(temp, n = 121, fill = NA, align = "center"),
      rolling_sd = RcppRoll::roll_sd(temp, n = 121, fill = NA, align = "center"),
      z_score = (temp - rolling_mean) / rolling_sd,
      temp_z_filtered = if_else(abs(z_score) > 3, NA_real_, temp)
    ) %>%
    mutate(
      temp_smoothed_30s = RcppRoll::roll_median(
        temp_z_filtered,
        n = 21,
        fill = NA,
        align = "center"
      )
    ) %>%
    mutate(
      temp_residuals_from_activity_30s = {
        if (sum(!is.na(temp) & !is.na(enmoTrunc)) >= 10 && var(enmoTrunc, na.rm = TRUE) > 1e-9) {
          model <- lm(temp ~ enmoTrunc, data = cur_data(), na.action = na.exclude)
          residuals(model)
        } else {
          NA_real_
        }
      }
    ) %>%
    ungroup() %>%
    filter(!is.na(temp_smoothed_30s)) %>%
    mutate(timestamp_5min = floor_date(time, "5 minutes")) %>%
    group_by(id, timestamp_5min) %>%
    summarise(
      temp_final = mean(temp_smoothed_30s, na.rm = TRUE),
      temp_adjusted_final = mean(temp_residuals_from_activity_30s, na.rm = TRUE),
      enmo_final = mean(enmoTrunc, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    rename(timestamp = timestamp_5min)
}

calculate_cosinor_params <- function(df, temp_col_name) {
  df <- as.data.frame(df) %>% arrange(timestamp)
  df$time_hours <- as.numeric(difftime(df$timestamp, min(df$timestamp), units = "hours"))
  df$temperature <- df[[temp_col_name]]
  df_model <- df %>% filter(!is.na(temperature) & !is.na(time_hours))
  if (nrow(df_model) < 20 || max(df_model$time_hours, na.rm=TRUE) < 24 || var(df_model$temperature, na.rm = TRUE) < 1e-6) {
    return(tibble(mesor = NA_real_, amplitude = NA_real_, acrophase_hr = NA_real_))
  }
  cosinor_model <- tryCatch({ lm(temperature ~ sin(2 * pi * time_hours / 24) + cos(2 * pi * time_hours / 24), data = df_model) }, error = function(e) NULL)
  if (is.null(cosinor_model) || length(coef(cosinor_model)) < 3) { return(tibble(mesor = NA_real_, amplitude = NA_real_, acrophase_hr = NA_real_)) }
  intercept <- coef(cosinor_model)["(Intercept)"]; beta_sin <- coef(cosinor_model)["sin(2 * pi * time_hours / 24)"]; beta_cos <- coef(cosinor_model)["cos(2 * pi * time_hours / 24)"]
  if (any(is.na(c(intercept, beta_sin, beta_cos)))) { return(tibble(mesor = NA_real_, amplitude = NA_real_, acrophase_hr = NA_real_)) }
  mesor <- intercept; amplitude <- sqrt(beta_sin^2 + beta_cos^2); acrophase_rad <- atan2(-beta_sin, beta_cos); acrophase_hr <- (acrophase_rad %% (2 * pi)) / (2 * pi) * 24
  return(tibble(mesor = mesor, amplitude = amplitude, acrophase_hr = acrophase_hr))
}

calculate_fft_params <- function(df, temp_col_name) {
  temp_series <- df[[temp_col_name]]
  temp_series <- temp_series[!is.na(temp_series)]
  if(all(is.na(temp_series))) return(tibble(fft_dominant_power = NA_real_, fft_dominant_period_hr = NA_real_))
  n <- length(temp_series)
  if (n <= 2) return(tibble(fft_dominant_power = NA_real_, fft_dominant_period_hr = NA_real_))
  fft_result <- fft(temp_series - mean(temp_series, na.rm=TRUE)); power <- Mod(fft_result[2:floor(n/2)])^2 / n
  sampling_rate <- 1 / (5/60); frequencies <- (1:(floor(n/2)-1)) * sampling_rate / n; periods_hr <- 1 / frequencies
  if(length(power) == 0 || all(is.na(power))) return(tibble(fft_dominant_power = NA_real_, fft_dominant_period_hr = NA_real_))
  dominant_power <- max(power, na.rm = TRUE); dominant_period_hr <- periods_hr[which.max(power)]
  return(tibble(fft_dominant_power = dominant_power, fft_dominant_period_hr = dominant_period_hr))
}

run_all_models_on_wear <- function(df) {
  cosinor_std <- calculate_cosinor_params(df, "temp_final"); names(cosinor_std) <- paste0(names(cosinor_std), "_std")
  cosinor_adj <- calculate_cosinor_params(df, "temp_adjusted_final"); names(cosinor_adj) <- paste0(names(cosinor_adj), "_adj")
  fft_std <- calculate_fft_params(df, "temp_final"); names(fft_std) <- paste0(names(fft_std), "_std")
  bind_cols(cosinor_std, cosinor_adj, fft_std)
}

get_first_peak <- function(df, col_name) {
  df <- df %>% arrange(time)
  series <- df[[col_name]]
  series <- series[!is.na(series)]
  N <- length(series)
  if (N < 10) {
    return(tibble(id = unique(df$id), freq = NA_real_, amp = NA_real_))
  }
  series_adjust <- series - mean(series)
  ft <- fft(series_adjust)
  freq <- (0:(N-1)) / N
  amplitude <- Mod(ft)
  half <- 2:floor(N/2)
  if(length(half) < 1) {
    return(tibble(id = unique(df$id), freq = NA_real_, amp = NA_real_))
  }
  amp_half <- amplitude[half]
  freq_half <- freq[half]
  peak_idx <- which(diff(sign(diff(amp_half))) == -2) + 1
  final_peak_idx <- NA_integer_
  if (length(peak_idx) > 0) {
    final_peak_idx <- peak_idx[1]
  } else if (length(amp_half) > 0 && !all(is.na(amp_half))) {
    final_peak_idx <- which.max(amp_half)
  }
  if (is.na(final_peak_idx) || length(final_peak_idx) == 0) {
    return(tibble(id = unique(df$id), freq = NA_real_, amp = NA_real_))
  }
  tibble(id = unique(df$id), freq = freq_half[final_peak_idx], amp = amp_half[final_peak_idx])
}


# --- 5. Main Batch Processing Engine (TEST MODE: ONE BATCH ONLY) ---
print("--- Starting: Main batch processing loop (TEST MODE) ---")

repeat_qc_participant_ids <- unique(qc_passed_repeats$participant_id)
participant_id_chunks <- split(repeat_qc_participant_ids, ceiling(seq_along(repeat_qc_participant_ids) / BATCH_SIZE_PARTICIPANTS))

processed_data_batch_full <- NULL
processed_data_batch_downsampled <- NULL
model_results_batch <- NULL

# This loop will only run once for the first batch and then break
for (i in seq_along(participant_id_chunks)) {
  current_participant_ids_batch <- participant_id_chunks[[i]]
  print(paste0(">>> Processing Batch ", i, " of ", length(participant_id_chunks), " <<<"))
  
  current_wear_ids_batch <- qc_passed_repeats %>%
    filter(participant_id %in% current_participant_ids_batch) %>%
    pull(id)
  
  print(paste0("  1. Loading raw data for ", length(current_wear_ids_batch), " wears..."))
  raw_data_batch <- load_raw_data_for_ids(current_wear_ids_batch, temp_dataset)
  
  if (nrow(raw_data_batch) == 0) {
      print("  No data loaded for this batch, skipping.")
      next
  }
  
  print("  2a. Pre-processing data (filtering, smoothing)...")
  processed_data_batch_full <- raw_data_batch %>%
    rename(time = event_timestamp) %>%
    group_by(id) %>%
    arrange(time) %>%
    mutate(
      rolling_mean = RcppRoll::roll_mean(temp, n = 121, fill = NA, align = "center"),
      rolling_sd = RcppRoll::roll_sd(temp, n = 121, fill = NA, align = "center"),
      z_score = (temp - rolling_mean) / rolling_sd,
      temp_z_filtered = if_else(abs(z_score) > 3, NA_real_, temp)
    ) %>%
    mutate(
      temp_smoothed = RcppRoll::roll_median(temp_z_filtered, n = 21, fill = NA, align = "center")
    ) %>%
    ungroup()

  print("  2b. Downsampling data to 5-minute epochs...")
  processed_data_batch_downsampled <- preprocess_and_downsample(raw_data_batch)

  print("  2c. Running models on downsampled data to get phenotypes...")
  nested_data_batch <- processed_data_batch_downsampled %>%
    group_by(id) %>% nest()
  
  model_results_batch <- nested_data_batch %>%
    mutate(model_output = map(data, run_all_models_on_wear)) %>%
    select(-data) %>%
    unnest(model_output)

  print("--- First batch processed. Exiting loop to begin plotting. ---")
  break
}


# --- 6. Plotting and Visualization ---
if (is.null(processed_data_batch_full)) {
  stop("No data was processed. Cannot generate plots.")
}
print("\n--- Starting: Generating all plots for the first batch ---")

# --- 6a. Individual 4-Panel Plots ---
print("  - Generating 4-panel plots for 4 random participants...")
print("    - Fitting batch-level mixed-effects model (lmer)...")
model_data_batch <- processed_data_batch_downsampled %>%
  mutate(
    hour_of_day = hour(timestamp) + minute(timestamp) / 60,
    time_radian = 2 * pi * hour_of_day / 24
  )

mlm_model_batch <- lmer(temp_final ~ sin(time_radian) + cos(time_radian) + (1 + sin(time_radian) + cos(time_radian) | id), data = model_data_batch)
print("    - Batch-level MLM fitting complete.")

set.seed(42)
random_ids <- sample(unique(processed_data_batch_full$id), 4)

for (k in 1:length(random_ids)) {
  current_id <- random_ids[k]
  cat(paste0("    - Generating plot for participant: ", current_id, "\n"))
  participant_data_vis <- processed_data_batch_full %>% filter(id == current_id)
  model_data_single <- model_data_batch %>% filter(id == current_id)
  lm_model <- lm(temp_final ~ sin(time_radian) + cos(time_radian), data = model_data_single)
  model_data_single$temp_predicted_lm <- predict(lm_model)
  model_data_single$temp_predicted_mlm <- predict(mlm_model_batch, newdata = model_data_single)
  plot_a <- ggplot(participant_data_vis, aes(x = time, y = temp)) + geom_line(color = "gray50", alpha = 0.8) + labs(title = "A: Raw Data", y = "Temp (°C)", x = "") + theme_bw()
  plot_b <- ggplot(participant_data_vis, aes(x = time)) + geom_line(aes(y = temp_z_filtered), color = "dodgerblue") + geom_point(data = . %>% filter(is.na(temp_z_filtered)), aes(y = temp), color = "red", size = 0.5) + labs(title = "B: Z-Score Filtered", y = "Temp (°C)", x = "") + theme_bw()
  plot_c <- ggplot(participant_data_vis, aes(x = time, y = temp_smoothed)) + geom_line(color = "navy") + labs(title = "C: Smoothed Signal", y = "Temp (°C)", x = "Time") + theme_bw()
  
  plot_d <- ggplot(model_data_single, aes(x = timestamp)) + 
    geom_point(aes(y = temp_final), color = "gray70", alpha = 0.5, size = 0.8) + 
    geom_line(aes(y = temp_predicted_lm, color = "Simple Model (LM)"), linewidth = 1.2) + 
    geom_line(aes(y = temp_predicted_mlm, color = "Mixed Model (MLM)"), linewidth = 1.2) + 
    scale_color_manual(values = c("Simple Model (LM)" = "darkorange", "Mixed Model (MLM)" = "purple")) + 
    labs(title = "D: Fitted Circadian Models on 5-min Data", subtitle = "Comparing individual (LM) vs. batch-level (MLM) fit", y = "Temp (°C)", x = "Time", color = "Model Type") + 
    theme_bw() + 
    theme(legend.position = "bottom")
  
  final_plot <- (plot_a | plot_b) / (plot_c | plot_d) + plot_annotation(title = paste("Preprocessing & Modeling Steps for Participant:", current_id))
  file_name <- file.path(OUTPUT_DIR, paste0("participant_", current_id, "_4panel_visualization.pdf"))
  ggsave(file_name, final_plot, width = 12, height = 9)
}

# --- 6b. Fourier Analysis Plots ---
print("  - Generating Fourier analysis plots...")
fourier_peaks_temp <- processed_data_batch_full %>% group_by(id) %>% group_split() %>% map_df(~get_first_peak(., "temp_smoothed"))
fourier_peaks_enmo <- processed_data_batch_full %>% group_by(id) %>% group_split() %>% map_df(~get_first_peak(., "enmoTrunc"))
names(fourier_peaks_enmo) <- c("id", "freq_enmo", "amp_enmo")

p_fft_freq <- ggplot(fourier_peaks_temp, aes(x = freq)) + geom_histogram(bins = 50, fill = "skyblue", color = "black") + labs(title = "Distribution of First Peak Frequency (Temperature)", x = "Frequency", y = "Count") + theme_bw()
p_fft_amp <- ggplot(fourier_peaks_temp, aes(x = amp)) + geom_histogram(bins = 50, fill = "orange", color = "black") + labs(title = "Distribution of First Peak Amplitude (Temperature)", x = "Amplitude", y = "Count") + theme_bw()
ggsave(file.path(OUTPUT_DIR, "fft_distributions.pdf"), p_fft_freq / p_fft_amp, width = 8, height = 8)

# --- 6c. 24-Hour Profile Plots ---
print("  - Generating 24-hour profile plots...")
processed_data_batch_full <- processed_data_batch_full %>% mutate(hour = hour(time) + minute(time)/60)
hourly_summary_temp <- processed_data_batch_full %>% group_by(hour) %>% summarise(mean_val = mean(temp_smoothed, na.rm = TRUE), se_val = sd(temp_smoothed, na.rm = TRUE) / sqrt(n()), .groups = 'drop')
hourly_summary_enmo <- processed_data_batch_full %>% group_by(hour) %>% summarise(mean_val = mean(enmoTrunc, na.rm = TRUE), se_val = sd(enmoTrunc, na.rm = TRUE) / sqrt(n()), .groups = 'drop')

p_24h_temp <- ggplot(hourly_summary_temp, aes(x = hour, y = mean_val)) + geom_line(color = "blue", linewidth = 1) + geom_ribbon(aes(ymin = mean_val - se_val, ymax = mean_val + se_val), alpha = 0.2, fill = "blue") + scale_x_continuous(breaks = seq(0, 24, by = 2)) + labs(title = "Average Temperature Across 24-Hour Clock", x = "Hour of Day", y = "Mean Temperature (°C)") + theme_bw()
p_24h_enmo <- ggplot(hourly_summary_enmo, aes(x = hour, y = mean_val)) + geom_line(color = "darkgreen", linewidth = 1) + geom_ribbon(aes(ymin = mean_val - se_val, ymax = mean_val + se_val), alpha = 0.2, fill = "darkgreen") + scale_x_continuous(breaks = seq(0, 24, by = 2)) + labs(title = "Average Activity (ENMO) Across 24-Hour Clock", x = "Hour of Day", y = "Mean ENMO") + theme_bw()
ggsave(file.path(OUTPUT_DIR, "24hr_profiles_temp_enmo.pdf"), p_24h_temp / p_24h_enmo, width = 10, height = 8)

# --- 6d. Correlation and Alignment Plots ---
print("  - Generating correlation and alignment plots...")
activity_min <- min(hourly_summary_enmo$mean_val, na.rm = TRUE); activity_max <- max(hourly_summary_enmo$mean_val, na.rm = TRUE)
temp_min <- min(hourly_summary_temp$mean_val, na.rm = TRUE); temp_max <- max(hourly_summary_temp$mean_val, na.rm = TRUE)
hourly_summary_enmo <- hourly_summary_enmo %>% mutate(mean_enmo_inv = activity_max - mean_val, mean_enmo_inv_rescaled = (mean_enmo_inv - min(mean_enmo_inv, na.rm=T)) / (max(mean_enmo_inv, na.rm=T) - min(mean_enmo_inv, na.rm=T)) * (temp_max - temp_min) + temp_min)

p_inv_activity <- ggplot() + geom_line(data = hourly_summary_temp, aes(x = hour, y = mean_val, color = "Temperature"), linewidth = 1) + geom_ribbon(data = hourly_summary_temp, aes(x=hour, ymin = mean_val - se_val, ymax = mean_val + se_val), alpha = 0.2, fill = "blue") + geom_line(data = hourly_summary_enmo, aes(x = hour, y = mean_enmo_inv_rescaled, color = "Inverse Activity"), linewidth = 1, linetype = "dashed") + scale_color_manual(values=c("Temperature"="blue", "Inverse Activity"="red")) + scale_x_continuous(breaks = seq(0, 24, by = 2)) + labs(title = "Temperature and Inverse Activity Across 24-Hour Clock", x = "Hour of Day", y = "Mean Temperature (°C)", color="Metric") + theme_bw()
ggsave(file.path(OUTPUT_DIR, "temp_vs_inverse_activity.pdf"), p_inv_activity, width = 10, height = 6)

sampled_ids <- sample(unique(processed_data_batch_full$id), min(100, length(unique(processed_data_batch_full$id))))
indiv_corrs <- map_dbl(sampled_ids, function(pid) {
    df <- processed_data_batch_full %>% filter(id == pid) %>% group_by(hour) %>% summarise(mean_temp = mean(temp_smoothed, na.rm = TRUE), mean_enmo = mean(enmoTrunc, na.rm = TRUE))
    if (nrow(df) >= 5) { cor(df$mean_temp, df$mean_enmo, use = "complete.obs") } else { NA_real_ }
})
p_indiv_corr_hist <- ggplot(data.frame(correlation=indiv_corrs), aes(x = correlation)) + geom_histogram(binwidth = 0.05, fill = "skyblue", color = "black") + labs(title = "Distribution of Per-Individual Temp-Activity Correlations", x = "Correlation", y = "Count") + theme_bw()
ggsave(file.path(OUTPUT_DIR, "individual_correlations_histogram.pdf"), p_indiv_corr_hist, width = 8, height = 5)

processed_data_aligned <- processed_data_batch_full %>%
  group_by(id) %>% arrange(time) %>%
  mutate(time_since_start = as.numeric(difftime(time, min(time), units = "hours")), temp_z = as.numeric(scale(temp_smoothed))) %>%
  ungroup() %>% filter(!is.na(temp_smoothed), abs(temp_z) <= 5)
p_all_aligned <- ggplot(processed_data_aligned, aes(x = time_since_start, y = temp_smoothed, group = id)) + geom_line(alpha = 0.3, color = "blue") + labs(title = "All Participants (Batch 1): Temperature Aligned at Start", x = "Time Since Start (hours)", y = "Temperature (°C)") + theme_bw()
ggsave(file.path(OUTPUT_DIR, "all_participants_aligned.pdf"), p_all_aligned, width = 12, height = 6)

# --- 6e. R-squared analysis ---
print("  - Generating R-squared analysis plots...")
r2_by_id <- processed_data_batch_full %>%
  filter(!is.na(temp_smoothed), !is.na(enmoTrunc)) %>%
  group_by(id) %>%
  do(broom::glance(lm(temp_smoothed ~ enmoTrunc, data = .))) %>%
  ungroup()
p_r2_hist <- ggplot(r2_by_id, aes(x = r.squared)) + geom_histogram(bins = 30, fill = "skyblue", color = "black") + labs(title = "Histogram of R-squared: temp_smoothed ~ enmoTrunc (per participant)", x = "R-squared", y = "Count") + theme_bw()
ggsave(file.path(OUTPUT_DIR, "r_squared_histogram.pdf"), p_r2_hist, width = 8, height = 6)


# --- 7. Correlation Matrix of Derived Phenotypes ---
print("  - Generating correlation matrix plot of derived phenotypes...")

if (!is.null(model_results_batch) && nrow(model_results_batch) > 0) {
  corr_data <- model_results_batch %>%
    select(ends_with("_std"), ends_with("_adj")) %>%
    select_if(is.numeric) %>%
    drop_na()

  if (nrow(corr_data) > 2) {
    correlation_matrix <- cor(corr_data)
    pdf(file.path(OUTPUT_DIR, "phenotype_correlation_matrix.pdf"), width = 10, height = 10)
    corrplot(correlation_matrix,
             method = "circle", type = "upper", tl.col = "black", tl.srt = 45,
             addCoef.col = "black", number.cex = 0.7,
             col = colorRampPalette(c("blue", "white", "red"))(200))
    dev.off()
    print("    - Correlation matrix saved.")
  } else {
    print("    - Not enough complete data to generate correlation matrix.")
  }
} else {
  print("    - No model results generated, skipping correlation matrix.")
}

# --- 8. Intra-class Correlation (ICC) Analysis ---
print("  - Generating Intra-class Correlation (ICC) plots for phenotype stability...")

icc_data_full <- model_results_batch %>%
  left_join(select(qc_passed_repeats, id, participant_id), by = "id")

phenotype_cols <- icc_data_full %>%
  select(ends_with("_std"), ends_with("_adj")) %>%
  names()

icc_results_real <- map_df(phenotype_cols, function(pheno_col) {
  df <- icc_data_full %>%
    select(participant_id, value = !!sym(pheno_col)) %>%
    na.omit()
  valid_subjects <- df %>% count(participant_id) %>% filter(n > 1) %>% pull(participant_id)
  if (length(valid_subjects) < 2) return(tibble(phenotype = pheno_col, icc = NA_real_))
  df_filtered <- df %>% filter(participant_id %in% valid_subjects)
  icc_val <- icc(df_filtered, model = "oneway", type = "consistency", unit = "single")$value
  tibble(phenotype = pheno_col, icc = icc_val)
}, .id = NULL)

set.seed(123)
icc_results_null <- map_df(phenotype_cols, function(pheno_col) {
  df <- icc_data_full %>%
    select(participant_id, value = !!sym(pheno_col)) %>%
    na.omit()
  df$participant_id <- sample(df$participant_id)
  valid_subjects <- df %>% count(participant_id) %>% filter(n > 1) %>% pull(participant_id)
  if (length(valid_subjects) < 2) return(tibble(phenotype = pheno_col, icc = NA_real_))
  df_filtered <- df %>% filter(participant_id %in% valid_subjects)
  icc_val <- icc(df_filtered, model = "oneway", type = "consistency", unit = "single")$value
  tibble(phenotype = pheno_col, icc = icc_val)
}, .id = NULL)

icc_results_real$type <- "Real ICC"
icc_results_null$type <- "Null ICC (Shuffled)"
icc_plot_data <- bind_rows(icc_results_real, icc_results_null)

p_icc <- ggplot(icc_plot_data, aes(x = icc, y = reorder(phenotype, icc), color = type)) +
  geom_point(size = 4) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  facet_wrap(~type, scales = "free_x") +
  labs(
    title = "Phenotype Stability: Intra-class Correlation (ICC)",
    subtitle = "Comparing real stability vs. a null model (shuffled IDs)",
    x = "ICC Value",
    y = "Phenotype"
  ) +
  theme_bw() +
  theme(legend.position = "none")

ggsave(file.path(OUTPUT_DIR, "phenotype_stability_ICC.pdf"), p_icc, width = 12, height = 8)
print("    - ICC plot saved.")

# --- 9. Create a Tarball of All Plots ---
print("--- Creating a tarball of the results directory ---")
tar_file_name <- file.path(dirname(OUTPUT_DIR), "results_test_plots.tar.gz")
tryCatch({
  tar(tarfile = tar_file_name, files = OUTPUT_DIR, compression = 'gzip')
  print(paste("Successfully created tarball:", tar_file_name))
}, error = function(e) {
  print(paste("Error creating tarball:", e$message))
})

print("--- All plots generated successfully. Script finished. ---")