# ==============================================================================
# R SCRIPT FOR LARGE-SCALE BATCH PHENOTYPE GENERATION
#
# Description:
# This script is designed for robust, parallelized, large-scale batch processing
# of the entire temperature dataset. It is intended to be run on a high-
# performance computing cluster.
#
# Key features:
#   1. No Sample QC: Processes all available participant records without applying
#      any quality control filters.
#   2. Block Processing & Restart: Iterates through the master list of participants
#      in user-defined blocks and can be restarted from any block.
#   3. Parallel Execution: Uses the 'furrr' package to process participants
#      within each block in parallel, maximizing CPU usage.
#   4. Comprehensive Phenotyping: Calculates the full suite of ~20 phenotypes
#      (cosinor, distributional, FFT, etc.) for each participant.
#   5. Automatic Upload: Saves results for each block and automatically uploads
#      them to a specified DNAnexus path.
#
# ==============================================================================


# --- 1. SETUP: LOAD LIBRARIES AND DEFINE CONSTANTS ---

# Ensure all necessary packages are installed and loaded.
packages <- c("tidyverse", "arrow", "data.table", "lubridate", "RcppRoll", "broom", "furrr")
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
OUTPUT_DIR            <- "phenotype_blocks" # Directory to save output files

# Define batch processing parameters
BLOCK_SIZE <- 500 # Process 100 participants per job/block. A smaller size is safer for memory.
START_BLOCK <- 1   # <<-- EDIT THIS VALUE TO RESTART FROM A SPECIFIC BLOCK

# Define DNAnexus upload path
DX_UPLOAD_PATH <- "/temp_ukb_cohorts" # <<-- EDIT THIS TO YOUR TARGET PROJECT/FOLDER

# Create the output directory
dir.create(OUTPUT_DIR, showWarnings = FALSE)

# Set up parallel processing plan
# This will use all available cores on the machine.
plan(multisession)
cat(sprintf("Parallel processing enabled with %d workers.\n", nbrOfWorkers()))


# --- 2. LOAD AND PREPARE MASTER & EXTERNAL DATA ---

cat("Loading master lists and external data...\n")

# Load the covariates file to get the full list of all measurement IDs
all_ids_to_process <- fread(COVARIATE_FILE_PATH)

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
participant_location_lookup <- all_ids_to_process %>%
  mutate(eid = as.integer(participant_id)) %>%
  left_join(location_clean, by = "eid") %>%
  select(id, participant_id, assessment_centre)


# --- 3. DEFINE THE MASTER PROCESSING FUNCTION ---

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
  
  # Combine all metrics
  final_metrics <- bind_cols(model_metrics, fft_temp_metrics, fft_enmo_metrics, tibble(r_squared_activity))
  
  # Add identifiers and return
  final_metrics %>%
    add_column(id = target_id, .before = 1) %>%
    add_column(participant_id = location_info$participant_id, .before = 1)
}


# --- 4. EXECUTE BATCH PROCESSING ---

# Determine the total number of blocks to process
total_ids <- nrow(all_ids_to_process)
n_blocks <- ceiling(total_ids / BLOCK_SIZE)

cat(sprintf("\nTotal participants to process: %d\n", total_ids))
cat(sprintf("Processing in %d blocks of size %d...\n", n_blocks, BLOCK_SIZE))

# Loop through the data in blocks, starting from START_BLOCK
for (i in START_BLOCK:n_blocks) {
  
  start_row <- (i - 1) * BLOCK_SIZE + 1
  end_row <- min(i * BLOCK_SIZE, total_ids)
  
  cat(sprintf("\n--- Starting Block %d of %d (Rows %d to %d) ---\n", i, n_blocks, start_row, end_row))
  
  # Get the list of IDs for the current block
  block_ids <- all_ids_to_process$id[start_row:end_row]
  
  # Use future_map_dfr for parallel processing within the block
  # The 'safely' wrapper prevents the entire block from failing if one ID has an error.
  safe_process <- safely(process_participant)
  
  block_results_list <- block_ids %>%
    future_map(~ safe_process(.x, participant_info = participant_location_lookup, ambient_temp_data = temps_at_location_clean), .progress = TRUE)
  
  # Separate successful results from errors for this block
  block_successful_results <- block_results_list %>% map("result") %>% compact() %>% bind_rows()
  block_errors <- block_results_list %>% map("error") %>% compact()
  
  if (length(block_errors) > 0) {
    cat(sprintf("\nWarning: Encountered %d errors during processing of this block.\n", length(block_errors)))
  }
  
  # --- Save the results and upload ---
  if (nrow(block_successful_results) > 0) {
    # Define the output filename for the current block
    output_filename <- file.path(OUTPUT_DIR, sprintf("phenotypes_%d-%d.tsv.gz", start_row, end_row))
    cat(sprintf("Saving %d successful results for this block to %s\n", nrow(block_successful_results), output_filename))
    fwrite(block_successful_results, file = output_filename, sep = "\t")
    
    # Construct and execute the dx upload command
    dx_command <- sprintf("dx upload %s --path %s", output_filename, DX_UPLOAD_PATH)
    cat("Uploading to DNAnexus with command:\n", dx_command, "\n")
    system(dx_command)
  } else {
    cat("No successful results to save for this block.\n")
  }
}

cat("\n--- ALL BATCH PROCESSING COMPLETE ---\n")
