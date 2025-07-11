# ==============================================================================
# R SCRIPT FOR BATCH PHENOTYPE GENERATION & ICC CALCULATION
#
# Description:
# This script is designed for efficient, non-interactive batch processing.
# It performs the following key steps:
#   1. Identifies participants with repeat measurements from the covariates file.
#   2. Filters for measurements that pass a moderately strict per-sample quality control.
#   3. Defines a master function to process a single participant's data from
#      loading through to full phenotype calculation.
#   4. Iterates over the QC-passed, repeat-measurement participants.
#   5. Saves the final table of all derived phenotypes to a gzipped CSV file.
#   6. Runs a diagnostic check to report on sample size for different analyses.
#   7. Calculates the Intraclass Correlation Coefficient (ICC) for each
#      phenotype at multiple levels of adjustment to assess stability.
#   8. Performs a permutation test on the ICCs as a negative control.
#
# ==============================================================================


# --- 1. SETUP: LOAD LIBRARIES AND DEFINE CONSTANTS ---

# Ensure all necessary packages are installed and loaded.
packages <- c("tidyverse", "arrow", "data.table", "lubridate", "RcppRoll", "broom", "irr")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

# --- Configuration ---
# Define file paths for all necessary inputs.
COVARIATE_FILE_PATH   <- "combined_covariates.csv.gz"
TEMPERATURE_DATA_PATH <- "combined_temperature_data.parquet"
LOCATION_FILE_PATH    <- "data_participant.tsv"
OUTPUT_FILE           <- "repeat_participant_metrics.csv.gz"


# --- 2. IDENTIFY REPEAT PARTICIPANTS & APPLY QC ---

cat("Identifying participants with repeat measurements and applying moderate QC...\n")

# Load the covariates file which links the complex 'id' to the stable 'participant_id'
covars <- fread(COVARIATE_FILE_PATH)

# Find the first 300 participant_ids that have more than one measurement
repeat_participant_ids <- covars %>%
  count(participant_id, name = "n_repeats") %>%
  filter(n_repeats > 1) %>%
  slice_head(n = 300) %>%
  pull(participant_id)

# Define a threshold for calibration error (remove worst 2%)
cal_err_threshold <- quantile(covars$calErrAfter_mg, 0.98, na.rm = TRUE)
cat(sprintf("Calibration error threshold (98th percentile): %.2f mg\n", cal_err_threshold))

# Get the specific 'id' values for these repeat participants, applying the new QC criteria.
ids_to_process <- covars %>%
  filter(
    participant_id %in% repeat_participant_ids,
    goodWear == 1,
    goodCal == 1,
    wearDays >= 5, # Wear day requirement
    calErrAfter_mg < cal_err_threshold # Exclude top 2% of calibration errors
  )

cat(sprintf("Identified %d participants with repeats. After moderate QC, processing %d measurements.\n",
            length(repeat_participant_ids), nrow(ids_to_process)))


# --- 3. PRE-LOAD AND PREPARE EXTERNAL DATA (FOR EFFICIENCY) ---

cat("Pre-loading and preparing external location and temperature data...\n")

# Load external data files ONCE to avoid reading them inside the loop.
location <- as.tibble(fread(LOCATION_FILE_PATH))
temps_at_location <- fread("https://github.com/gushamilton/ukb_temp_epochs/raw/refs/heads/main/data/ukb_temps_per_date.tsv.gz")

# Clean location data
location_clean <- location %>%
  select(eid, assessment_centre = p54_i0)

# Clean and de-duplicate ambient temperature data
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

# Create a lookup table for participant location details
participant_location_lookup <- ids_to_process %>%
  mutate(eid = as.integer(participant_id)) %>%
  left_join(location_clean, by = "eid") %>%
  select(id, participant_id, assessment_centre)


# --- 4. DEFINE THE MASTER PROCESSING FUNCTION ---

# This single function encapsulates all steps to process one participant 'id'.
process_participant <- function(target_id, participant_info, ambient_temp_data) {
  
  # Helper function to calculate FFT metrics
  get_first_peak <- function(data_vec) {
    if (length(na.omit(data_vec)) < 20) return(tibble(peak_freq = NA, peak_amp = NA))
    vec_adjusted <- data_vec - mean(data_vec, na.rm = TRUE); vec_adjusted[is.na(vec_adjusted)] <- 0
    N <- length(vec_adjusted); ft <- fft(vec_adjusted); freq <- (0:(N-1)) / N; amplitude <- Mod(ft) / (N/2)
    half_indices <- 2:floor(N/2); amp_half <- amplitude[half_indices]; freq_half <- freq[half_indices]
    peak_indices <- which(diff(sign(diff(amp_half))) == -2) + 1
    first_peak_idx <- if (length(peak_indices) == 0) which.max(amp_half) else peak_indices[1]
    tibble(peak_freq = freq_half[first_peak_idx], peak_amp = amp_half[first_peak_idx])
  }
  
  # Helper function to calculate cosinor & other metrics
  calculate_other_metrics <- function(df) {
    # Distribution-based phenotypes
    p5 <- quantile(df$temp_final, 0.05, na.rm = TRUE)
    p95 <- quantile(df$temp_final, 0.95, na.rm = TRUE)
    temp_range <- p95 - p5
    relative_amplitude <- (p95 - p5) / (p95 + p5)
    
    # Model 1: Unadjusted Cosinor
    unadj_model <- lm(temp_final ~ cos(2*pi*hour/24) + sin(2*pi*hour/24), data = df)
    coeffs_unadj <- coef(unadj_model)
    if (any(is.na(coeffs_unadj)) || length(coeffs_unadj) < 3) { mesor <- NA; amplitude <- NA; acrophase_h <- NA
    } else {
      mesor <- coeffs_unadj[1]; amplitude <- sqrt(coeffs_unadj[2]^2 + coeffs_unadj[3]^2); acrophase_rad <- atan2(coeffs_unadj[3], coeffs_unadj[2]); acrophase_h <- (-acrophase_rad * (24 / (2 * pi))) %% 24
    }
    
    # Model 2: Activity-Adjusted Cosinor
    df_act <- df %>% mutate(enmo_centered = enmo_final - mean(enmo_final, na.rm = TRUE))
    adj_model <- lm(temp_final ~ enmo_centered + cos(2*pi*hour/24) + sin(2*pi*hour/24), data = df_act)
    coeffs_adj <- coef(adj_model)
    if (any(is.na(coeffs_adj)) || length(coeffs_adj) < 4) { mesor_adj <- NA; amplitude_adj <- NA; acrophase_h_adj <- NA; activity_effect_slope <- NA
    } else {
      mesor_adj <- coeffs_adj[1]; amplitude_adj <- sqrt(coeffs_adj[3]^2 + coeffs_adj[4]^2); acrophase_rad_adj <- atan2(coeffs_adj[4], coeffs_adj[3]); acrophase_h_adj <- (-acrophase_rad_adj * (24 / (2 * pi))) %% 24
      activity_effect_slope <- coeffs_adj[2]
    }
    
    # Model 3: Ambient Temperature-Adjusted Cosinor
    df_amb <- df %>% mutate(ambient_temp_centered = ambient_temp - mean(ambient_temp, na.rm = TRUE))
    if (all(is.na(df_amb$ambient_temp_centered))) {
      mesor_ambient <- NA; amplitude_ambient <- NA; acrophase_h_ambient <- NA; environmental_sensitivity <- NA; residual_variability <- NA
    } else {
      ambient_model <- lm(temp_final ~ ambient_temp_centered + cos(2*pi*hour/24) + sin(2*pi*hour/24), data = df_amb)
      coeffs_ambient <- coef(ambient_model)
      if (any(is.na(coeffs_ambient)) || length(coeffs_ambient) < 4) { mesor_ambient <- NA; amplitude_ambient <- NA; acrophase_h_ambient <- NA; environmental_sensitivity <- NA; residual_variability <- NA
      } else {
        mesor_ambient <- coeffs_ambient[1]; amplitude_ambient <- sqrt(coeffs_ambient[3]^2 + coeffs_ambient[4]^2); acrophase_rad_ambient <- atan2(coeffs_ambient[4], coeffs_ambient[3]); acrophase_h_ambient <- (-acrophase_rad_ambient * (24 / (2 * pi))) %% 24
        environmental_sensitivity <- coeffs_ambient[2]
        residual_variability <- sd(residuals(ambient_model), na.rm = TRUE)
      }
    }
    
    return(tibble(mesor, amplitude, acrophase_h, mesor_adj, amplitude_adj, acrophase_h_adj, mesor_ambient, amplitude_ambient, acrophase_h_ambient, temp_range, relative_amplitude, activity_effect_slope, environmental_sensitivity, residual_variability))
  }
  
  # --- Main function logic starts here ---
  cat(sprintf("Processing ID: %s\n", target_id))
  
  # Load data for the specific ID
  temp_dataset <- open_dataset(TEMPERATURE_DATA_PATH)
  raw_data <- temp_dataset %>%
    filter(id == target_id) %>%
    select(id, time, temp, enmoTrunc) %>%
    collect() %>%
    mutate(time = as_datetime(time))
    
  if (nrow(raw_data) == 0) return(NULL) # Skip if no data found
  
  # Pre-processing
  processed_data <- raw_data %>%
    group_by(id) %>%
    arrange(time) %>%
    mutate(
      rolling_mean = RcppRoll::roll_mean(temp, n = 21, fill = NA, align = "center"),
      rolling_sd = RcppRoll::roll_sd(temp, n = 21, fill = NA, align = "center"),
      z_score = (temp - rolling_mean) / rolling_sd,
      temp_z_filtered = if_else(abs(z_score) > 3, NA_real_, temp)
    ) %>%
    mutate(temp_smoothed = RcppRoll::roll_median(temp_z_filtered, n = 7, fill = NA, align = "center")) %>%
    ungroup()

  # Downsampling
  downsampled_data <- processed_data %>%
    filter(!is.na(temp_smoothed)) %>%
    mutate(time_5min = floor_date(time, unit = "5 minutes")) %>%
    group_by(id, time_5min) %>%
    summarise(
      temp_final = mean(temp_smoothed, na.rm = TRUE),
      enmo_final = mean(enmoTrunc, na.rm = TRUE),
      hour = first(hour(time)),
      .groups = 'drop'
    )
    
  # Link ambient temperature
  location_info <- participant_info %>% filter(id == target_id)
  
  data_with_ambient <- downsampled_data %>%
    mutate(date = as_date(time_5min), assessment_centre = location_info$assessment_centre) %>%
    left_join(ambient_temp_data, by = c("assessment_centre", "date"))
    
  # Calculate cosinor and other metrics
  model_metrics <- calculate_other_metrics(data_with_ambient)
  
  # Calculate FFT metrics
  fft_temp_metrics <- get_first_peak(processed_data$temp_smoothed) %>% rename_with(~paste0("fft_", .x, "_temp"))
  fft_enmo_metrics <- get_first_peak(processed_data$enmoTrunc) %>% rename_with(~paste0("fft_", .x, "_enmo"))
  
  # Calculate R-squared
  r2_fit <- lm(temp_smoothed ~ enmoTrunc, data = processed_data)
  r_squared_activity <- summary(r2_fit)$r.squared
  
  # Get representative date and mean ambient temp for this measurement period
  mean_date <- as_date(mean(data_with_ambient$time_5min, na.rm = TRUE))
  mean_ambient_temp_week <- mean(data_with_ambient$ambient_temp, na.rm = TRUE)
  
  # Combine all metrics
  final_metrics <- bind_cols(model_metrics, fft_temp_metrics, fft_enmo_metrics, tibble(r_squared_activity), tibble(mean_date), tibble(mean_ambient_temp_week))
  
  # Add identifiers and return
  final_metrics %>%
    add_column(id = target_id, .before = 1) %>%
    add_column(participant_id = location_info$participant_id, .before = 1)
}


# --- 5. EXECUTE BATCH PROCESSING ---

cat("\n--- Starting Batch Processing ---\n")

# Use map to apply the processing function to each ID and row-bind the results.
safe_process <- safely(process_participant)
all_results <- ids_to_process$id %>%
  map(~ safe_process(.x, participant_info = participant_location_lookup, ambient_temp_data = temps_at_location_clean))

# Separate successful results from errors
successful_results <- all_results %>% map("result") %>% compact() %>% bind_rows()
errors <- all_results %>% map("error") %>% compact()

if (length(errors) > 0) {
  cat(sprintf("\nWarning: Encountered %d errors during processing.\n", length(errors)))
}

cat("\n--- Batch Processing Complete ---\n")


# --- 6. SAVE RESULTS ---

cat(sprintf("Saving %d successful results to %s\n", nrow(successful_results), OUTPUT_FILE))
fwrite(successful_results, file = OUTPUT_FILE)


# --- 7. DIAGNOSTIC: CHECK SAMPLE SIZES FOR ICC ---

cat("\n--- Diagnostic Check for ICC Sample Sizes ---\n")

# Count participants with at least two measurements in the final dataset
n_total_repeats <- successful_results %>%
  count(participant_id) %>%
  filter(n > 1) %>%
  nrow()

# Count participants with at least two measurements that have valid ambient temp data
n_ambient_repeats <- successful_results %>%
  filter(!is.na(mean_ambient_temp_week)) %>%
  count(participant_id) %>%
  filter(n > 1) %>%
  nrow()

cat(sprintf("Total participants with >= 2 repeats: %d\n", n_total_repeats))
cat(sprintf("Participants with >= 2 repeats AND valid ambient temp data: %d\n", n_ambient_repeats))
cat("The difference explains the lower 'n' for ambient-adjusted ICCs.\n")


# --- 8. CALCULATE INTRACLASS CORRELATION COEFFICIENT (ICC) ---

cat("\n--- Calculating ICC for Repeat Measurements ---\n")

# Helper function to perform ICC calculation and formatting
calculate_icc <- function(data, group_col, time_col, value_col) {
  wide_data <- data %>%
    select({{group_col}}, {{time_col}}, {{value_col}}) %>%
    pivot_wider(names_from = {{time_col}}, values_from = {{value_col}}, names_prefix = "time_")
  
  icc_values <- wide_data %>% select(starts_with("time_"))
  
  if (ncol(icc_values) < 2) return(tibble(phenotype = value_col, n = NA, icc = NA, lower_ci = NA, upper_ci = NA, p_value = NA))
  
  icc_result <- tryCatch(
    icc(icc_values, model = "twoway", type = "agreement", unit = "single"),
    error = function(e) NULL
  )
  
  if (!is.null(icc_result)) {
    tibble(
      phenotype = value_col,
      n = icc_result$subjects,
      icc = icc_result$value,
      lower_ci = icc_result$lbound,
      upper_ci = icc_result$ubound,
      p_value = icc_result$p.value
    )
  } else {
    tibble(phenotype = value_col, n = NA, icc = NA, lower_ci = NA, upper_ci = NA, p_value = NA)
  }
}

# Prepare data for ICC calculation
icc_data <- successful_results %>%
  group_by(participant_id) %>%
  mutate(time_point = row_number()) %>%
  ungroup()

# Get list of numeric phenotypes to test
phenotype_cols <- icc_data %>%
  select(where(is.numeric), -participant_id, -time_point, -mean_ambient_temp_week) %>%
  colnames()

# --- 8a. Calculate TRUE ICC ---
cat("\n--- True ICC Results (Stability of Raw Phenotypes) ---\n")
true_icc_summary <- map_df(phenotype_cols, ~calculate_icc(icc_data, "participant_id", "time_point", .x))
print(true_icc_summary %>% arrange(desc(icc)))

# --- 8b. Calculate Seasonally-Adjusted ICC ---
cat("\n--- Seasonally-Adjusted ICC Results ---\n")
seasonally_adjusted_data <- icc_data %>%
  mutate(doy = yday(mean_date))
for (pheno in phenotype_cols) {
  model <- lm(get(pheno) ~ sin(2*pi*doy/365) + cos(2*pi*doy/365), data = seasonally_adjusted_data, na.action = na.exclude)
  seasonally_adjusted_data[[paste0(pheno, "_seas_adj")]] <- residuals(model)
}
adjusted_phenotype_cols <- paste0(phenotype_cols, "_seas_adj")
seasonally_adjusted_icc_summary <- map_df(adjusted_phenotype_cols, ~calculate_icc(seasonally_adjusted_data, "participant_id", "time_point", .x))
print(seasonally_adjusted_icc_summary %>% arrange(desc(icc)))


# --- 8c. Calculate Fully-Adjusted ICC (Season + Mean Ambient Temp) ---
cat("\n--- Fully-Adjusted (Season + Ambient Temp) ICC Results ---\n")
fully_adjusted_data <- seasonally_adjusted_data %>%
  mutate(mean_ambient_temp_week = scale(mean_ambient_temp_week)) # Scale for model stability

for (pheno in phenotype_cols) {
  seas_adj_pheno <- paste0(pheno, "_seas_adj")
  full_adj_pheno <- paste0(pheno, "_full_adj")
  
  if (seas_adj_pheno %in% names(fully_adjusted_data)) {
    model <- lm(get(seas_adj_pheno) ~ mean_ambient_temp_week, data = fully_adjusted_data, na.action = na.exclude)
    fully_adjusted_data[[full_adj_pheno]] <- residuals(model)
  }
}
fully_adjusted_phenotype_cols <- paste0(phenotype_cols, "_full_adj")
fully_adjusted_icc_summary <- map_df(fully_adjusted_phenotype_cols, ~calculate_icc(fully_adjusted_data, "participant_id", "time_point", .x))
print(fully_adjusted_icc_summary %>% arrange(desc(icc)))


# --- 8d. Calculate PERMUTED (NULL) ICC ---
cat("\n--- Permuted (Null) ICC Results (Negative Control) ---\n")
permuted_icc_data <- icc_data %>%
  group_by(time_point) %>%
  mutate(participant_id = sample(participant_id)) %>%
  ungroup()
permuted_icc_summary <- map_df(phenotype_cols, ~calculate_icc(permuted_icc_data, "participant_id", "time_point", .x))
print(permuted_icc_summary %>% arrange(desc(icc)))


cat("\n--- SCRIPT FINISHED SUCCESSFULLY ---\n")
