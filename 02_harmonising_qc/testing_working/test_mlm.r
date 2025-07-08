# --- 1. Load Libraries ---
# Use pacman to load/install packages. It will install them if they are not already present.
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, arrow, RcppRoll, patchwork, lubridate, lme4)

# --- 2. Load and Prepare Data ---
# This script assumes the initial data extraction (Part 2 of the research plan)
# has already been completed, resulting in a parquet file.
file_path <- "/Users/fh6520/Downloads/temperature_epochs_first100.parquet"
if (!file.exists(file_path)) {
  stop("Data file not found. Please update the 'file_path' variable.")
}
d <- read_parquet(file_path)

# --- 3. Full Pre-processing Pipeline ---
# This section implements the core of the pre-processing strategy (Part 2.2 and 2.3 of the plan).
# The goal is to transform the raw, noisy 30-second data into a clean, model-ready time series.

cat("Processing all loaded participants...\n")

processed_data <- d %>%
  # Convert time column from character to a proper datetime object for time-based operations.
  mutate(time = as_datetime(time)) %>%
  # Group by each participant. All subsequent rolling window operations will be contained within each individual.
  group_by(id) %>%
  
  # Step A: Z-score Filtering to remove local, non-physiological spikes.
  # This is a robust outlier detection method that adapts to the local signal characteristics.
  # A rolling Z-score is preferred over a global one as it is not biased by overall trends.
  mutate(
    rolling_mean = RcppRoll::roll_mean(temp, n = 121, fill = NA, align = "center"),
    rolling_sd = RcppRoll::roll_sd(temp, n = 121, fill = NA, align = "center"),
    z_score = (temp - rolling_mean) / rolling_sd,
    temp_z_filtered = if_else(abs(z_score) > 3, NA_real_, temp)
  ) %>%
  
  # Step B: Smoothing the Z-score filtered data to remove high-frequency noise.
  # We use a rolling median because it is robust to any remaining outliers.
  # RcppRoll is used for high performance on large datasets.
  mutate(
    temp_smoothed = RcppRoll::roll_median(
      temp_z_filtered,
      n = 21,          # A 10-minute window (21 points * 30s/point ≈ 10 mins)
      fill = NA,
      align = "center"
    )
  ) %>%
  ungroup() # Ungroup when finished with per-individual calculations.

cat("Processing complete.\n")


# --- 4. Downsample Data to 5-Minute Epochs ---
# The primary purpose of this step is to reduce the data volume and computational load
# for the modeling stage, while retaining the core circadian signal revealed by smoothing.
cat("Downsampling data to 5-minute epochs...\n")

downsampled_data <- processed_data %>%
  # We can't use time points where the smoothed temperature is NA.
  filter(!is.na(temp_smoothed)) %>%
  group_by(id) %>%
  # Create a new time column representing the start of each 5-minute interval.
  mutate(time_5min = floor_date(time, unit = "5 minutes")) %>%
  # Group by the new 5-minute epoch to aggregate the points within it.
  group_by(id, time_5min) %>%
  # Calculate the mean of the now-clean, smoothed temperature for the final value.
  summarise(
    temp_final = mean(temp_smoothed, na.rm = TRUE),
    .groups = 'drop' # Drop grouping structure after summarising.
  )

cat("Downsampling complete. Data size reduced significantly.\n")


# --- 5. Fit Global Mixed-Effects Model (MLM) on Downsampled Data ---
# This implements the primary modeling strategy (Part 3.2 of the plan).
# The LMM allows us to estimate the average circadian rhythm for the population (fixed effects)
# while also estimating each individual's unique deviation from that average (random effects).
# This "borrows strength" across individuals, providing more robust estimates.
cat("Fitting global mixed-effects model (lmer)...\n")

# Prepare data for the MLM by creating harmonic time variables (cosine and sine).
model_data_all <- downsampled_data %>%
  mutate(
    hour_of_day = hour(time_5min) + minute(time_5min) / 60,
    time_radian = 2 * pi * hour_of_day / 24
  )

# Fit the MLM. The formula specifies that temp_final depends on the overall rhythm
# (sin + cos terms), and that each 'id' has its own random intercept (mean),
# and its own random slopes for the sin and cos terms (rhythm shape).
mlm_model <- lmer(temp_final ~ sin(time_radian) + cos(time_radian) + (1 + sin(time_radian) + cos(time_radian) | id), data = model_data_all)

cat("MLM fitting complete.\n")


# --- 6. Loop Through Four Random Participants for Visualization ---
# This section generates the multi-panel plots to visually inspect the results
# of each stage of the pipeline and demonstrate key concepts.
set.seed(42) # For reproducibility
random_ids <- sample(unique(processed_data$id), 4)

for (i in 1:length(random_ids)) {
  
  current_id <- random_ids[i]
  cat("\n--- Generating plot", i, "for participant:", current_id, "---\n")
  
  # Isolate data for the current participant for plotting.
  participant_data_vis <- processed_data %>% filter(id == current_id)
  model_data_single <- model_data_all %>% filter(id == current_id)
  
  # Fit a simple Linear Model (LM) for comparison. This model only sees this one person's data.
  lm_model <- lm(temp_final ~ sin(time_radian) + cos(time_radian), data = model_data_single)
  model_data_single$temp_predicted_lm <- predict(lm_model)
  
  # Get predictions from the global MLM. These predictions incorporate the random effects for this specific person.
  model_data_single$temp_predicted_mlm <- predict(mlm_model, newdata = model_data_single)
  
  # --- Create the Plots ---
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
    labs(title = "D: Models on Full Data", subtitle = "Note: With high-quality data, fits are similar", y = "Temp (°C)", x = "Time", color = "Model Type") +
    theme_bw() + theme(legend.position = "bottom")
  
  # --- Create Plot E for the first participant to demonstrate shrinkage ---
  # This plot shows why an MLM is superior when data is sparse or noisy.
  if (i == 1) {
    cat("--- Creating shrinkage demonstration plot (E) ---\n")
    # Artificially create sparse data by sampling 5% of the points.
    sparse_data <- model_data_single %>% sample_frac(0.05)
    
    # Fit a simple LM just to this sparse data. It will likely be a poor, over-fitted model.
    lm_sparse <- lm(temp_final ~ sin(time_radian) + cos(time_radian), data = sparse_data)
    model_data_single$temp_predicted_lm_sparse <- predict(lm_sparse, newdata = model_data_single)
    
    plot_e <- ggplot(model_data_single, aes(x = time_5min)) +
      geom_point(data = sparse_data, aes(y = temp_final), color = "gray70", alpha = 0.8) +
      geom_line(aes(y = temp_predicted_lm_sparse, color = "LM on Sparse Data"), linewidth = 1.2) +
      geom_line(aes(y = temp_predicted_mlm, color = "MLM Prediction"), linewidth = 1.2) +
      scale_color_manual(values = c("LM on Sparse Data" = "red", "MLM Prediction" = "purple")) +
      labs(title = "E: Shrinkage Demonstrated on Sparse (5%) Data", subtitle = "MLM 'borrows strength' and remains stable", y = "Temp (°C)", x = "Time", color = "Model Type") +
      theme_bw() + theme(legend.position = "bottom")
      
    # Arrange plots with the special 5th plot.
    final_plot <- (plot_a | plot_b) / (plot_c | plot_d) / plot_e +
      plot_annotation(title = paste("Preprocessing & Modeling Steps for Participant:", current_id))
    
  } else {
    # For other participants, use the original 4-panel layout.
    final_plot <- (plot_a | plot_b) / (plot_c | plot_d) +
      plot_annotation(title = paste("Preprocessing & Modeling Steps for Participant:", current_id))
  }
  
  file_name <- paste0("participant_", i, "_visualization.pdf")
  ggsave(file_name, final_plot, width = 12, height = 11)
  cat("Saved plot to", file_name, "\n")
}

