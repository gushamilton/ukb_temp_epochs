# ==============================================================================
# R SCRIPT FOR LARGE-SCALE BATCH PHENOTYPE GENERATION
#
# Updated 2025-07-10 – *simplified*: uploads are now fire-and-forget
#                     (no post-upload verification or state check)
#
#  • Parallel batch processing of UKB temperature records
#  • Plain dx upload that preserves filenames by writing to a folder ending “/”
#  • Results land in project-J1GYbG0JQz0gxQ6yZf1GqYkb:/temp_ukb_cohorts/
#  • Restartable by editing START_BLOCK
# ==============================================================================


# --- 1  SETUP: LOAD LIBRARIES & CONSTANTS ------------------------------------

packages <- c(
  "tidyverse", "arrow", "data.table", "lubridate", "RcppRoll",
  "broom", "furrr"
)
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, repos = "http://cran.us.r-project.org")
    library(pkg, character.only = TRUE)
  }
}

# ---- Paths and batch parameters ---------------------------------------------
COVARIATE_FILE_PATH   <- "combined_covariates.csv.gz"
TEMPERATURE_DATA_PATH <- "combined_temperature_data.parquet"
LOCATION_FILE_PATH    <- "data_participant.tsv"
OUTPUT_DIR            <- "phenotype_blocks"          # local directory for TSVs

BLOCK_SIZE   <- 1000     # participants per block
START_BLOCK  <- 1       # change to restart
DX_FOLDER    <- "/temp_ukb_cohorts/"                 # ***trailing slash required***

# ---- DNAnexus project --------------------------------------------------------
DX_PROJECT_ID <- Sys.getenv(
  "DX_PROJECT_ID",
  unset = "project-J1GYbG0JQz0gxQ6yZf1GqYkb"         # fallback for local tests
)

dir.create(OUTPUT_DIR, showWarnings = FALSE)

plan(multisession)
cat(sprintf("Parallel processing enabled with %d workers.\n", nbrOfWorkers()))


# --- 2  HELPER: MINIMAL UPLOADER ---------------------------------------------
#     • preserves original filename
#     • no verification beyond dx exit status
# -----------------------------------------------------------------------------

upload_file <- function(local_file, project_id, dest_folder) {
  remote_path <- sprintf("%s:%s%s",
                         project_id,
                         dest_folder,
                         basename(local_file))       # keeps local name
  dx_args <- c("upload", local_file, "--destination", remote_path)
  message("→ ", paste(c("dx", dx_args), collapse = " "))
  exit_code <- system2("dx", dx_args)
  if (exit_code != 0L)
    stop("dx upload failed (exit ", exit_code, ").")
  invisible(NULL)
}


# --- 3  LOAD MASTER & EXTERNAL DATA ------------------------------------------

cat("Loading master lists and external data…\n")

all_ids_to_process <- fread(COVARIATE_FILE_PATH)

location <- as_tibble(fread(LOCATION_FILE_PATH))
temps_at_location <- fread(
  "https://github.com/gushamilton/ukb_temp_epochs/raw/refs/heads/main/data/ukb_temps_per_date.tsv.gz"
)

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
    .groups = "drop"
  )

participant_location_lookup <- all_ids_to_process %>%
  mutate(eid = as.integer(participant_id)) %>%
  left_join(location_clean, by = "eid") %>%
  select(id, participant_id, assessment_centre)


# --- 4  MASTER PROCESSING FUNCTION -------------------------------------------

process_participant <- function(target_id, participant_info, ambient_temp_data) {

  # Helper: FFT peak ----------------------------------------------------------
  get_first_peak <- function(vec) {
    if (length(na.omit(vec)) < 20)
      return(tibble(peak_freq = NA, peak_amp = NA))
    adj <- vec - mean(vec, na.rm = TRUE); adj[is.na(adj)] <- 0
    N   <- length(adj)
    amp <- Mod(fft(adj)) / (N / 2)
    freq <- (0:(N - 1)) / N
    idx  <- 2:floor(N / 2)
    peak <- which(diff(sign(diff(amp[idx]))) == -2) + 1
    p <- if (length(peak)) peak[1] else which.max(amp[idx])
    tibble(peak_freq = freq[idx][p], peak_amp = amp[idx][p])
  }

  # Helper: cosinor + distribution metrics -----------------------------------
  calc_metrics <- function(df) {
    p5  <- quantile(df$temp_final, 0.05, na.rm = TRUE)
    p95 <- quantile(df$temp_final, 0.95, na.rm = TRUE)
    temp_range <- p95 - p5
    rel_amp    <- (p95 - p5) / (p95 + p5)

    cu <- coef(lm(temp_final ~ cos(2*pi*hour/24) + sin(2*pi*hour/24), data = df))
    if (any(is.na(cu)) || length(cu) < 3) {
      mesor <- amp <- acrophase <- NA
    } else {
      mesor     <- cu[1]
      amp       <- sqrt(cu[2]^2 + cu[3]^2)
      acrophase <- (-atan2(cu[3], cu[2]) * 24 / (2*pi)) %% 24
    }

    df_act <- df %>% mutate(enmo_c = enmo_final - mean(enmo_final, na.rm = TRUE))
    ca <- coef(lm(temp_final ~ enmo_c + cos(2*pi*hour/24) + sin(2*pi*hour/24), data = df_act))
    if (any(is.na(ca)) || length(ca) < 4) {
      mesor_adj <- amp_adj <- acrophase_adj <- act_slope <- NA
    } else {
      mesor_adj     <- ca[1]
      amp_adj       <- sqrt(ca[3]^2 + ca[4]^2)
      acrophase_adj <- (-atan2(ca[4], ca[3]) * 24 / (2*pi)) %% 24
      act_slope     <- ca[2]
    }

    df_amb <- df %>% mutate(amb_c = ambient_temp - mean(ambient_temp, na.rm = TRUE))
    if (all(is.na(df_amb$amb_c))) {
      mesor_amb <- amp_amb <- acrophase_amb <- env_sens <- resid_sd <- NA
    } else {
      cb <- coef(lm(temp_final ~ amb_c + cos(2*pi*hour/24) + sin(2*pi*hour/24), data = df_amb))
      if (any(is.na(cb)) || length(cb) < 4) {
        mesor_amb <- amp_amb <- acrophase_amb <- env_sens <- resid_sd <- NA
      } else {
        mesor_amb     <- cb[1]
        amp_amb       <- sqrt(cb[3]^2 + cb[4]^2)
        acrophase_amb <- (-atan2(cb[4], cb[3]) * 24 / (2*pi)) %% 24
        env_sens      <- cb[2]
        resid_sd      <- sd(residuals(lm(temp_final ~ amb_c + cos(2*pi*hour/24) + sin(2*pi*hour/24), data = df_amb)), na.rm = TRUE)
      }
    }

    tibble(
      mesor, amp, acrophase,
      mesor_adj, amp_adj, acrophase_adj,
      mesor_amb, amp_amb, acrophase_amb,
      temp_range, rel_amp,
      act_slope, env_sens, resid_sd
    )
  }

  # --- Load data -------------------------------------------------------------
  temp_ds <- open_dataset(TEMPERATURE_DATA_PATH)
  raw <- temp_ds %>%
    filter(id == target_id) %>%
    select(id, time, temp, enmoTrunc) %>%
    collect() %>%
    mutate(time = as_datetime(time))
  if (nrow(raw) == 0) return(NULL)

  # --- Pre-processing --------------------------------------------------------
  proc <- raw %>%
    arrange(time) %>%
    group_by(id) %>%
    mutate(
      rolling_mean = RcppRoll::roll_mean(temp, 21, fill = NA, align = "center"),
      rolling_sd   = RcppRoll::roll_sd(temp,   21, fill = NA, align = "center"),
      z            = (temp - rolling_mean) / rolling_sd,
      temp_filt    = if_else(abs(z) > 3, NA_real_, temp),
      temp_smooth  = RcppRoll::roll_median(temp_filt, 7, fill = NA, align = "center")
    ) %>%
    ungroup()

  down <- proc %>%
    filter(!is.na(temp_smooth)) %>%
    mutate(t5 = floor_date(time, "5 minutes")) %>%
    group_by(id, t5) %>%
    summarise(
      temp_final = mean(temp_smooth, na.rm = TRUE),
      enmo_final = mean(enmoTrunc,   na.rm = TRUE),
      hour       = first(hour(time)),
      .groups = "drop"
    )

  # --- Join ambient temp -----------------------------------------------------
  loc <- participant_info %>% filter(id == target_id)
  with_amb <- down %>%
    mutate(
      date = as_date(t5),
      assessment_centre = loc$assessment_centre
    ) %>%
    left_join(ambient_temp_data, by = c("assessment_centre", "date"))

  # --- Metrics ---------------------------------------------------------------
  metrics     <- calc_metrics(with_amb)
  fft_temp    <- get_first_peak(proc$temp_smooth) %>% rename_with(~paste0("fft_", .x, "_temp"))
  fft_enmo    <- get_first_peak(proc$enmoTrunc)   %>% rename_with(~paste0("fft_", .x, "_enmo"))
  r2_activity <- summary(lm(temp_smooth ~ enmoTrunc, data = proc))$r.squared

  bind_cols(
    tibble(id = target_id, participant_id = loc$participant_id),
    metrics, fft_temp, fft_enmo,
    tibble(r_squared_activity = r2_activity)
  )
}


# --- 5  EXECUTE BATCH PROCESSING ---------------------------------------------

total_ids <- nrow(all_ids_to_process)
n_blocks  <- ceiling(total_ids / BLOCK_SIZE)

cat(sprintf("\nTotal participants: %d — processing in %d blocks of %d.\n",
            total_ids, n_blocks, BLOCK_SIZE))

for (blk in START_BLOCK:n_blocks) {
  s <- (blk - 1) * BLOCK_SIZE + 1
  e <- min(blk * BLOCK_SIZE, total_ids)

  cat(sprintf("\n--- Block %d/%d (rows %d-%d) ---\n", blk, n_blocks, s, e))

  ids <- all_ids_to_process$id[s:e]
  safe_proc <- safely(process_participant)

  results <- ids %>%
    future_map(
      ~ safe_proc(.x,
                  participant_info  = participant_location_lookup,
                  ambient_temp_data = temps_at_location_clean),
      .progress = TRUE
    )

  ok  <- results %>% map("result") %>% compact() %>% bind_rows()
  err <- results %>% map("error")  %>% compact()

  if (length(err))
    cat(sprintf("⚠️  %d errors in this block (skipped).\n", length(err)))

  if (nrow(ok)) {
    out_file <- file.path(
      OUTPUT_DIR,
      sprintf("phenotypes_%d-%d.tsv.gz", s, e)
    )
    fwrite(ok, out_file, sep = "\t")

    # --- simple upload -------------------------------------------------------
    upload_file(
      local_file = out_file,
      project_id = DX_PROJECT_ID,
      dest_folder = DX_FOLDER
    )
    cat("✔ Uploaded", basename(out_file), "to",
        sprintf("%s:%s\n", DX_PROJECT_ID, DX_FOLDER))
  } else {
    cat("No valid results in this block.\n")
  }
}

cat("\n--- ALL BATCH PROCESSING COMPLETE ---\n")
