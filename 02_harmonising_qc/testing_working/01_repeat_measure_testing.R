# Title: Analyze Stability of Circadian Rhythms in Repeat Participants (Batch Processing)
# Author: Gemini
# Date: 2025-07-03
# Description: This script loads temperature data for participants with multiple wears
#              from a large Parquet file. It processes the data in batches to manage
#              memory, applies rigorous cleaning (outlier filtering, smoothing),
#              downsamples to 5-minute epochs, and runs circadian models to derive
#              phenotypes. The final output is a tidy data frame containing model
#              parameters for each valid wear, ready for stability analysis (ICC/Pearson's R).
#
# Key Decisions based on user feedback:
# - QC: Filters for wears >= 4 days.
# - FFT: Calculated on the final, smoothed, 5-minute epoch data to align with cosinor models.
# - Seasonality: A 'season' variable is derived and included in the final output.

# --- 1. Setup & Configuration ---

# Load all necessary libraries. Ensure they are installed.
# install.packages(c("arrow", "dplyr", "data.table", "tidyr", "purrr", "lubridate", "zoo", "RcppRoll"))
library(arrow)
library(dplyr)
library(data.table)
library(tidyr)
library(purrr)
library(lubridate)
library(zoo)
library(RcppRoll)

# --- Define Global Variables & Paths ---
# Please update these paths to match your local environment
COVARIATE_FILE_PATH <- "/mnt/project/ukb_temp_final/combined_covariates.csv.gz"
TEMPERATURE_DATA_PATH <- "/mnt/project/ukb_temp_final/combined_temperature_data.parquet"
OUTPUT_DIR <- "/mnt/project/ukb_temp_final/results"

# Create the output directory if it doesn't exist
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Define the number of participants to process in each batch.
# Adjust this based on your workstation's available RAM.
BATCH_SIZE_PARTICIPANTS <- 100

# --- 2. Load Covariates and Select Participants ---

print("--- Starting: Loading and filtering covariate data ---")

# Load the gzipped CSV file of covariates using the efficient fread function
covars <- fread(COVARIATE_FILE_PATH)

# Identify participant_ids that have more than one entry (i.e., repeat measurements)
participant_counts <- covars %>% count(participant_id)
repeat_participant_ids <- participant_counts %>% filter(n > 1) %>% pull(participant_id)

# Filter the main covariates file to keep all entries for these repeat participants
repeat_covars <- covars %>%
  filter(participant_id %in% repeat_participant_ids)

print(paste0("Found ", length(repeat_participant_ids), " participants with repeat measurements."))

# Apply Quality Control (QC) to the repeat participant data.
# We only want to model data from high-quality wears with sufficient duration.
qc_passed_repeats <- repeat_covars %>%
  filter(goodWear == 1, goodCal == 1, wearDays >= 4) # User-defined QC: wearDays >= 4

# --- Add Seasonality ---
# Function to determine season from a date (Northern Hemisphere meteorological seasons)
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

# Ensure startTime is a date object and derive the season for each wear
qc_passed_repeats <- qc_passed_repeats %>%
  mutate(
    startTime = as_datetime(startTime),
    season = get_season(startTime)
  )

# Get the list of unique measurement IDs ('id') that passed QC
qc_ids <- unique(qc_passed_repeats$id)

print(paste0("After QC (wearDays >= 4), analyzing ", length(qc_ids), " individual measurements from ",
            length(unique(qc_passed_repeats$participant_id)), " participants."))
print("--- Finished: Covariate processing complete ---")


# --- 3. Open Parquet Dataset ---

# Open the Parquet dataset. This is memory-efficient as it doesn't load the whole file.
print("--- Starting: Opening Parquet dataset ---")
temp_dataset <- open_dataset(TEMPERATURE_DATA_PATH)
print("Parquet dataset opened. Raw data will be processed in batches.")


# --- 4. Modular Data Processing & Modeling Functions ---

# Function to load raw data for a given set of IDs from the Parquet dataset
load_raw_data_for_ids <- function(ids_to_load, temp_dataset_ref) {
  temp_dataset_ref %>%
    filter(id %in% ids_to_load) %>%
    collect() %>%
    rename(event_timestamp = time) %>%
    mutate(event_timestamp = as_datetime(event_timestamp))
}

# Comprehensive function to preprocess raw data and downsample
preprocess_and_downsample <- function(raw_data_df) {
  raw_data_df %>%
    group_by(id) %>%
    # Step A: Z-score Filtering on raw 30s data to remove local spikes
    mutate(
      rolling_mean = RcppRoll::roll_mean(temp, n = 121, fill = NA, align = "center"), # 60.5 min window
      rolling_sd = RcppRoll::roll_sd(temp, n = 121, fill = NA, align = "center"),
      z_score = (temp - rolling_mean) / rolling_sd,
      temp_z_filtered = if_else(abs(z_score) > 3, NA_real_, temp)
    ) %>%
    # Step B: Smoothing the Z-score filtered data using a rolling median
    mutate(
      temp_smoothed_30s = RcppRoll::roll_median(
        temp_z_filtered,
        n = 21, # 10.5 min window
        fill = NA,
        align = "center"
      )
    ) %>%
    # Step C: Calculate activity-adjusted temperature residuals at the 30s level
    mutate(
      temp_residuals_from_activity_30s = {
        # Check for sufficient data and variance to fit a linear model
        if (sum(!is.na(temp) & !is.na(enmoTrunc)) >= 10 && var(enmoTrunc, na.rm = TRUE) > 1e-9) {
          residuals(lm(temp ~ enmoTrunc, data = cur_data(), na.action = na.omit))
        } else {
          NA_real_
        }
      }
    ) %>%
    ungroup() %>%
    # Step D: Downsample to 5-minute epochs
    filter(!is.na(temp_smoothed_30s)) %>% # Remove NAs before aggregation
    mutate(timestamp_5min = floor_date(event_timestamp, "5 minutes")) %>%
    group_by(id, timestamp_5min) %>%
    summarise(
      temp_final = mean(temp_smoothed_30s, na.rm = TRUE),
      temp_adjusted_final = mean(temp_residuals_from_activity_30s, na.rm = TRUE),
      enmo_final = mean(enmoTrunc, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    rename(timestamp = timestamp_5min)
}

# Function to fit cosinor model and extract parameters
calculate_cosinor_params <- function(df, temp_col_name) {
  # Ensure data is a dataframe and sorted by time
  df <- as.data.frame(df) %>% arrange(timestamp)
  
  # Create a time variable in hours from the start of the measurement
  df$time_hours <- as.numeric(difftime(df$timestamp, min(df$timestamp), units = "hours"))
  
  # Use the specified temperature column
  df$temperature <- df[[temp_col_name]]
  
  # Filter out NA values for modeling
  df_model <- df %>% filter(!is.na(temperature) & !is.na(time_hours))
  
  # Check for sufficient data and variance to fit the model
  if (nrow(df_model) < 20 || max(df_model$time_hours, na.rm=TRUE) < 24 || var(df_model$temperature, na.rm = TRUE) < 1e-6) {
    return(tibble(mesor = NA_real_, amplitude = NA_real_, acrophase_hr = NA_real_))
  }
  
  # Fit a standard 24-hour cosinor model
  cosinor_model <- tryCatch({
      lm(temperature ~ sin(2 * pi * time_hours / 24) + cos(2 * pi * time_hours / 24), data = df_model)
  }, error = function(e) NULL)
  
  if (is.null(cosinor_model) || length(coef(cosinor_model)) < 3) {
    return(tibble(mesor = NA_real_, amplitude = NA_real_, acrophase_hr = NA_real_))
  }
  
  # Extract coefficients
  intercept <- coef(cosinor_model)["(Intercept)"]
  beta_sin <- coef(cosinor_model)["sin(2 * pi * time_hours / 24)"]
  beta_cos <- coef(cosinor_model)["cos(2 * pi * time_hours / 24)"]
  
  if (any(is.na(c(intercept, beta_sin, beta_cos)))) {
    return(tibble(mesor = NA_real_, amplitude = NA_real_, acrophase_hr = NA_real_))
  }
  
  # Calculate metrics
  mesor <- intercept
  amplitude <- sqrt(beta_sin^2 + beta_cos^2)
  acrophase_rad <- atan2(-beta_sin, beta_cos) # Note: traditional acrophase uses -sin
  acrophase_hr <- (acrophase_rad %% (2 * pi)) / (2 * pi) * 24 # Convert to 0-24h scale
  
  return(tibble(mesor = mesor, amplitude = amplitude, acrophase_hr = acrophase_hr))
}

# Function to calculate FFT metrics from the 5-min epoch data
calculate_fft_params <- function(df, temp_col_name) {
  temp_series <- df[[temp_col_name]]
  if(all(is.na(temp_series))) return(tibble(fft_dominant_power = NA_real_, fft_dominant_period_hr = NA_real_))
  
  temp_series <- temp_series[!is.na(temp_series)]
  n <- length(temp_series)
  if (n <= 2) return(tibble(fft_dominant_power = NA_real_, fft_dominant_period_hr = NA_real_))
  
  # Perform FFT
  fft_result <- fft(temp_series - mean(temp_series))
  power <- Mod(fft_result[2:floor(n/2)])^2 / n
  
  # Calculate corresponding periods. Sampling interval is 5 mins = 1/12 hour.
  sampling_rate <- 1 / (5/60) # samples per hour
  frequencies <- (1:(floor(n/2)-1)) * sampling_rate / n
  periods_hr <- 1 / frequencies
  
  if(length(power) == 0) return(tibble(fft_dominant_power = NA_real_, fft_dominant_period_hr = NA_real_))
  
  dominant_power <- max(power, na.rm = TRUE)
  dominant_period_hr <- periods_hr[which.max(power)]
  
  return(tibble(fft_dominant_power = dominant_power, fft_dominant_period_hr = dominant_period_hr))
}


# Wrapper function to run all models on a single wear's data
run_all_models_on_wear <- function(df) {
  # Model 1: Standard circadian parameters
  cosinor_std <- calculate_cosinor_params(df, "temp_final")
  names(cosinor_std) <- paste0(names(cosinor_std), "_std")
  
  # Model 2: Activity-adjusted circadian parameters
  cosinor_adj <- calculate_cosinor_params(df, "temp_adjusted_final")
  names(cosinor_adj) <- paste0(names(cosinor_adj), "_adj")
  
  # Model 3: FFT on standard temperature
  fft_std <- calculate_fft_params(df, "temp_final")
  names(fft_std) <- paste0(names(fft_std), "_std")
  
  # Combine all results into a single-row tibble
  bind_cols(cosinor_std, cosinor_adj, fft_std)
}


# --- 5. Main Batch Processing Engine ---

print("--- Starting: Main batch processing and modeling loop ---")

# Get unique participant IDs that passed QC and have repeats
repeat_qc_participant_ids <- unique(qc_passed_repeats$participant_id)

# Create chunks of participant IDs to process in batches
participant_id_chunks <- split(repeat_qc_participant_ids, ceiling(seq_along(repeat_qc_participant_ids) / BATCH_SIZE_PARTICIPANTS))

# Use map_dfr to process each chunk and combine results into a single data frame
model_results <- map_dfr(seq_along(participant_id_chunks), function(i) {
  current_participant_ids_batch <- participant_id_chunks[[i]]
  
  print(paste0(">>> Processing Batch ", i, " of ", length(participant_id_chunks), " <<<"))
  
  # Get all individual wear IDs ('id') for the current batch of participants
  current_wear_ids_batch <- qc_passed_repeats %>%
    filter(participant_id %in% current_participant_ids_batch) %>%
    pull(id)
  
  # 1. Load raw data for this batch of wear IDs
  print(paste0("  1. Loading raw data for ", length(current_wear_ids_batch), " wears..."))
  raw_data_batch <- load_raw_data_for_ids(current_wear_ids_batch, temp_dataset)
  
  if (nrow(raw_data_batch) == 0) {
      print("  No data loaded for this batch, skipping.")
      return(NULL)
  }
  
  # 2. Preprocess the raw data and downsample to 5-minute epochs
  print("  2. Preprocessing data (filtering, smoothing, downsampling)...")
  processed_data_batch <- preprocess_and_downsample(raw_data_batch)
  
  # 3. Nest the processed data by wear 'id' and join with covariates
  print("  3. Nesting data and joining covariates...")
  nested_data_batch <- processed_data_batch %>%
    group_by(id) %>%
    nest() %>%
    left_join(qc_passed_repeats, by = "id")
    
  # 4. Run all models on each measurement in the batch
  print(paste0("  4. Running models on ", nrow(nested_data_batch), " individual measurements..."))
  
  # Use another map to apply the modeling function to each nested data frame
  results_for_batch <- nested_data_batch %>%
    mutate(model_output = map(data, run_all_models_on_wear)) %>%
    select(-data) %>% # Remove the large list-column of raw data
    unnest(model_output) # Unnest the single-row tibble of model results
  
  print(paste0("  Batch ", i, " complete."))
  return(results_for_batch)
})

print("--- Finished: All batches processed. Modeling complete. ---")


# --- 6. Final Analysis and Output ---

print("--- Starting: Final analysis and saving results ---")

# Display the first few rows of the full results table
print("--- Full Results Per Wear (Head) ---")
print(head(model_results))

# Group by participant_id to assess stability of metrics across multiple wears.
stability_analysis <- model_results %>%
  group_by(participant_id) %>%
  summarise(
    n_wears = n(),
    # Stability of standard Mesor
    mean_mesor_std = mean(mesor_std, na.rm = TRUE),
    sd_mesor_std = sd(mesor_std, na.rm = TRUE),
    cv_mesor_std = sd_mesor_std / mean_mesor_std,
    # Stability of standard Amplitude
    mean_amplitude_std = mean(amplitude_std, na.rm = TRUE),
    sd_amplitude_std = sd(amplitude_std, na.rm = TRUE),
    cv_amplitude_std = sd_amplitude_std / mean_amplitude_std,
    # Stability of activity-adjusted Mesor
    mean_mesor_adj = mean(mesor_adj, na.rm = TRUE),
    sd_mesor_adj = sd(mesor_adj, na.rm = TRUE),
    cv_mesor_adj = sd_mesor_adj / mean_mesor_adj
  ) %>%
  filter(n_wears >= 2) # Only include participants with at least 2 valid wears for stability analysis

# Display the first few rows of the stability analysis
print("--- Stability Analysis Across Wears (Head) ---")
print(head(stability_analysis))

# Save the results to CSV files in the specified output directory
fwrite(model_results, file.path(OUTPUT_DIR, "repeat_participants_full_model_results.csv"))
fwrite(stability_analysis, file.path(OUTPUT_DIR, "repeat_participants_stability_analysis.csv"))

print(paste("Results successfully saved to:", OUTPUT_DIR))
print("--- Script finished. ---")
