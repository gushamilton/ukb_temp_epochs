# ==============================================================================
# R SCRIPT TO GENERATE GWAS-READY INPUT FILES FOR REGENIE
#
# This script performs two main tasks:
#   1. Downloads and processes batched phenotype files, applying QC filters
#      and calculating the mean for individuals with repeat measurements.
#   2. Processes the master participant data to create a covariate file with
#      age, sex, genetic PCs, and derived season of assessment.
#
# IMPORTANT: This script filters phenotype data to include only participants
#            identified as 'Caucasian' based on the p22006 field.
#            The covariate file will include ALL Caucasian participants.
# ==============================================================================

# --- 1. SETUP: LOAD LIBRARIES & DEFINE CONSTANTS ------------------------------

# Ensure necessary packages are installed and loaded
packages <- c("tidyverse", "data.table", "lubridate")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, repos = "http://cran.us.r-project.org")
    library(pkg, character.only = TRUE)
  }
}

# --- Constants ---
DX_PROJECT_ID                   <- "project-J1GYbG0JQz0gxQ6yZf1GqYkb"
DX_PHENOTYPE_DIR                <- "/temp_ukb_cohorts/"
DX_COMBINED_COVARIATES_FILE     <- "/temp_ukb_cohorts/combined_covariates.csv.gz"
LOCAL_PHENOTYPE_DIR   <- file.path(getwd(), "phenotype_batches") # Absolute path for robustness
LOCAL_COMBINED_COVARIATES_FILE  <- "combined_covariates.csv.gz"
MASTER_COVARIATE_FILE           <- "data_participant.tsv"
FINAL_PHENO_FILE                <- "pheno.tsv"
FINAL_COVAR_FILE                <- "covars.tsv"

dir.create(LOCAL_PHENOTYPE_DIR, showWarnings = FALSE)

# --- Define a function to get season from date ---
get_season_vars <- function(input_date) {
  if (is.na(input_date)) return(list(season_cos = NA_real_, season_sin = NA_real_))
  day_of_year <- yday(as.Date(input_date))
  # Normalize to 2*pi range
  day_rad <- day_of_year * (2 * pi / 365.25)
  list(season_cos = cos(day_rad), season_sin = sin(day_rad))
}

# --- Phenotype Priority for Correlation Filtering ---
# When two phenotypes are highly correlated (>= CORRELATION_CUTOFF),
# the one appearing later in this list will be removed.
# Prioritize raw signals over adjusted, then other derived metrics.
PHENOTYPE_PRIORITY <- c(
  "mesor", "amp", "acrophase", # Raw circadian parameters
  "temp_range", "rel_amp",     # Distributional metrics
  "fft_peak_freq_temp", "fft_peak_amp_temp", # FFT for temperature
  "fft_peak_freq_enmo", "fft_peak_amp_enmo", # FFT for activity
  "r_squared_activity",         # Activity correlation
  "act_slope", "env_sens", "resid_sd", # Activity/environment sensitivity, residual variability
  "mesor_adj", "amp_adj", "acrophase_adj", # Activity-adjusted circadian parameters
  "mesor_amb", "amp_amb", "acrophase_amb"  # Ambient-adjusted circadian parameters
)

CORRELATION_CUTOFF <- 0.8 # Absolute correlation threshold for filtering


# --- 2. PHENOTYPES: DOWNLOAD, COMBINE, QC, AND PROCESS ------------------------

cat("--- Starting Phenotype Processing ---\n")

# --- Download all batched phenotype files from DNAnexus ---
cat("Finding and downloading phenotype files from DNAnexus...\n")
dx_find_cmd <- sprintf("dx find data --project %s --folder %s --name 'phenotypes_*.tsv.gz' --brief",
                       DX_PROJECT_ID, DX_PHENOTYPE_DIR)
pheno_files <- system(dx_find_cmd, intern = TRUE)

if (length(pheno_files) == 0) {
  stop("No phenotype files found in ", DX_PHENOTYPE_DIR)
}

for (file_path in pheno_files) {
  system(sprintf("dx download %s -o %s/", file_path, LOCAL_PHENOTYPE_DIR))
}
cat("Downloaded", length(pheno_files), "phenotype files.\n")


# --- Download combined_covariates.csv.gz (contains QC metrics) ---
cat("Downloading combined_covariates.csv.gz...\n")
system(sprintf("dx download %s -o %s", DX_COMBINED_COVARIATES_FILE, LOCAL_COMBINED_COVARIATES_FILE))
qc_data <- fread(LOCAL_COMBINED_COVARIATES_FILE)


# --- Read and combine all local phenotype files ---
local_files <- list.files(LOCAL_PHENOTYPE_DIR, pattern = "\\.tsv\\.gz$", full.names = TRUE)
all_phenos <- rbindlist(lapply(local_files, fread))

cat("Combined all phenotype files. Total rows:", nrow(all_phenos), "\n")


# --- Load master covariate file and filter by ancestry ---
cat("Loading master covariate file and filtering by ancestry...\n")
header_line <- fread(MASTER_COVARIATE_FILE, nrows = 1, header = FALSE)
covar_master <- fread(MASTER_COVARIATE_FILE, skip = 1)
setnames(covar_master, names(header_line), as.character(header_line[1,]))

# Filter for Caucasian participants (p22006 == "Caucasian")
caucasian_eids <- covar_master %>%
  filter(p22006 == "Caucasian") %>%
  pull(eid)

cat(sprintf("Filtered to %d Caucasian participants.\n", length(caucasian_eids)))


# --- Merge phenotypes with QC data and filter by Caucasian EIDs ---
merged_data <- merge(all_phenos, qc_data, by = c("id", "participant_id"), all.x = TRUE) %>%
  filter(participant_id %in% caucasian_eids)


# --- Apply QC Filters (Adjust these parameters as needed) ---
cat("Applying QC filters...\n")

# QC Parameters:
# MIN_WEAR_DAYS: Minimum number of wear days required.
# MAX_CAL_ERR_MG: Maximum calibration error (mg) allowed.
# MAX_CLIPS_AFTER: Maximum number of clips after calibration allowed.
# REQUIRE_GOOD_WEAR: TRUE to require good wear time (goodWear == 1).
# REQUIRE_GOOD_CAL: TRUE to require good calibration (goodCal == 1).

MIN_WEAR_DAYS   <- 3
MAX_CAL_ERR_MG  <- 1000 # A high initial threshold, adjust based on data distribution
MAX_CLIPS_AFTER <- 1000 # A high initial threshold, adjust based on data distribution
REQUIRE_GOOD_WEAR <- TRUE
REQUIRE_GOOD_CAL  <- FALSE

qc_filtered_phenos <- merged_data %>%
  filter(
    wearDays >= MIN_WEAR_DAYS,
    calErrAfter_mg < MAX_CAL_ERR_MG,
    clipsAfter < MAX_CLIPS_AFTER
  ) 

if (REQUIRE_GOOD_WEAR) {
  qc_filtered_phenos <- qc_filtered_phenos %>%
    filter(goodWear == 1)
}

if (REQUIRE_GOOD_CAL) {
  qc_filtered_phenos <- qc_filtered_phenos %>%
    filter(goodCal == 1)
}

cat(sprintf("Filtered to %d rows after QC. Unique participants: %d\n",
            nrow(qc_filtered_phenos), uniqueN(qc_filtered_phenos$participant_id)))


# --- Process repeats: group by participant and calculate mean ---
# Ensure all phenotype columns are numeric before summarizing
pheno_cols <- setdiff(names(qc_filtered_phenos), 
                      c("id", "participant_id", "goodWear", "goodCal", 
                        "calErrAfter_mg", "wearDays", "clipsAfter", "p22006",
                        "xOffset_g", "yOffset_g", "zOffset_g", "clipsBefore", 
                        "startTime", "endTime"))
qc_filtered_phenos[, (pheno_cols) := lapply(.SD, as.numeric), .SDcols = pheno_cols]

# Group by participant_id and calculate the mean for all other columns
mean_phenos <- qc_filtered_phenos[, lapply(.SD, mean, na.rm = TRUE),
                          by = .(participant_id),
                          .SDcols = pheno_cols]

cat("Averaged phenotypes for participants with repeat measures. Final rows:", nrow(mean_phenos), "\n")


# --- Apply Correlation Filtering ---
cat("Applying correlation filtering...\n")

# Function to perform correlation filtering
filter_correlated_phenos <- function(data, pheno_list, cutoff, priority_order) {
  current_phenos <- pheno_list
  
  # Loop until no highly correlated pairs remain
  repeat {
    cor_matrix <- abs(cor(data[, .SD, .SDcols = current_phenos], use = "pairwise.complete.obs"))
    diag(cor_matrix) <- 0 # Ignore self-correlation
    
    # Find the highest correlation above cutoff
    max_cor_val <- max(cor_matrix, na.rm = TRUE)
    
    if (is.na(max_cor_val) || max_cor_val < cutoff) {
      break # No more highly correlated pairs
    }
    
    # Find the pair with the highest correlation
    max_cor_idx <- which(cor_matrix == max_cor_val, arr.ind = TRUE)[1, ]
    pheno1 <- rownames(cor_matrix)[max_cor_idx[1]]
    pheno2 <- colnames(cor_matrix)[max_cor_idx[2]]
    
    # Determine which phenotype to remove based on priority
    # Lower index in priority_order means higher priority
    priority1 <- which(priority_order == pheno1)
    priority2 <- which(priority_order == pheno2)
    
    if (priority1 > priority2) {
      to_remove <- pheno1
    } else {
      to_remove <- pheno2
    }
    
    current_phenos <- setdiff(current_phenos, to_remove)
    cat(sprintf("  Removed '%s' due to high correlation (r=%.2f) with '%s'.\n", to_remove, max_cor_val, setdiff(c(pheno1, pheno2), to_remove)))
  }
  return(current_phenos)
}

final_pheno_list <- filter_correlated_phenos(mean_phenos, pheno_cols, CORRELATION_CUTOFF, PHENOTYPE_PRIORITY)

cat("✔ Final selected phenotypes after correlation filtering:\n")
cat(paste(final_pheno_list, collapse = ", "), "\n")

# --- Format for REGENIE: Add FID/IID and save ---
final_phenos <- mean_phenos %>%
  select(all_of(final_pheno_list), participant_id) %>%
  mutate(FID = participant_id, IID = participant_id) %>%
  select(FID, IID, everything(), -participant_id)

fwrite(final_phenos, FINAL_PHENO_FILE, sep = "\t", na = "NA")
cat("✔ Successfully created final phenotype file:", FINAL_PHENO_FILE, "\n")


# --- 3. COVARIATES: LOAD AND PROCESS ------------------------------------------

cat("\n--- Starting Covariate Processing ---\n")

# Filter master covariate file by Caucasian EIDs
final_covars <- covar_master %>%
  filter(p22006 == "Caucasian") %>%
  as_tibble() %>%
  select(
    eid = eid,
    sex = p31,
    age = p21022,
    assessment_date = p53_i0,
    PC1 = p22009_a1, PC2 = p22009_a2, PC3 = p22009_a3, PC4 = p22009_a4, PC5 = p22009_a5,
    PC6 = p22009_a6, PC7 = p22009_a7, PC8 = p22009_a8, PC9 = p22009_a9, PC10 = p22009_a10,
    p22006 # Keep p22006 for reference if needed, but not as a covariate
  ) %>%
  mutate(
    FID = eid,
    IID = eid,
    season_vars = map(assessment_date, get_season_vars),
    season_cos = map_dbl(season_vars, "season_cos"),
    season_sin = map_dbl(season_vars, "season_sin"),
    sex_numeric = if_else(sex == "Male", 1, if_else(sex == "Female", 0, NA_real_))
  ) %>%
  select(FID, IID, age, sex = sex_numeric, season_cos, season_sin, starts_with("PC"))

fwrite(final_covars, FINAL_COVAR_FILE, sep = "\t", na = "NA")
cat("✔ Successfully created final covariate file:", FINAL_COVAR_FILE, "\n")

cat("\n--- All processing complete. ---\n")