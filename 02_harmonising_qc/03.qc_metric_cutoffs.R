#!/usr/bin/env Rscript
## optimise_qc_allphenos.R ----------------------------------------------------
## * Computes ICC for *all* phenotypes after dropping |r| > 0.80 duplicates.
## * Uses ≤ max_repeats participants (default 400) for every ICC computation.
## * Single-core; final CSV sorted by eff_sample_size; one column per phenotype.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  lapply(c("data.table","tidyverse","lme4"), require, character.only = TRUE)
})

# ------------- user‐tunable params -----------------------------------------
max_repeats <- 400     # ≤ this many repeat participants for ICC calc
corr_cutoff <- 0.80    # drop phenotypes with |r| > cutoff
min_repeat_subs <- 20  # must have ≥ this many repeat subjects after QC

# ------------- file paths ---------------------------------------------------
COVAR_FILE <- "combined_covariates.csv.gz"
PHENO_DIR  <- "phenotype_blocks"
OUT_CSV    <- "qc_grid_results.csv"

# ------------- load data ----------------------------------------------------
covar <- fread(
  COVAR_FILE,
  select = c("id","participant_id",
             "goodWear","goodCal","calErrAfter_mg",
             "wearDays","clipsAfter")
)
ph_files <- list.files(PHENO_DIR, "^phenotypes_.*\\.tsv\\.gz$", full.names = TRUE)
stopifnot(length(ph_files) > 0)
phenos <- rbindlist(lapply(ph_files, fread), use.names = TRUE, fill = TRUE)

dt <- merge(phenos, covar, by = c("id","participant_id"))
rm(phenos, covar); gc()

# ------------- helper: ICC via mixed model ----------------------------------
icc_mixed <- function(data, ph) {
  d <- data[!is.na(get(ph)), .(participant_id, y = get(ph))]
  reps <- d[, .N, by = participant_id][N > 1, participant_id]
  if (length(reps) < min_repeat_subs) return(NA_real_)

  if (length(reps) > max_repeats)          # subsample for speed
    reps <- sample(reps, max_repeats)

  d <- d[participant_id %in% reps]

  as.numeric(
    tryCatch({
      m  <- lmer(y ~ 1 + (1|participant_id), data = d,
                 control = lmerControl(check.conv.singular = "ignore"))
      vb <- as.data.frame(VarCorr(m))$vcov[1]
      ve <- attr(VarCorr(m), "sc")^2
      vb / (vb + ve)
    }, error = function(e) NA_real_)
  )
}

# ------------- build loose baseline ----------------------------------------
baseline <- dt[ wearDays >= 3 & clipsAfter < 1000 ]
stopifnot(baseline[, .N, by = participant_id][N > 1, .N] >= min_repeat_subs)

num_cols <- names(baseline)[vapply(baseline, is.numeric, FALSE)]
drop_cols <- c("id","participant_id","goodWear","goodCal",
               "calErrAfter_mg","wearDays","clipsAfter")
ph_pool <- setdiff(num_cols, drop_cols)

# ---------- de-duplicate: remove |r| > corr_cutoff --------------------------
cor_mat <- abs(cor(baseline[, ..ph_pool], use = "pairwise.complete.obs"))
diag(cor_mat) <- 0
while (any(cor_mat > corr_cutoff, na.rm = TRUE)) {
  bad <- colnames(cor_mat)[which(cor_mat > corr_cutoff, arr.ind = TRUE)[1, 2]]
  ph_pool <- setdiff(ph_pool, bad)
  cor_mat <- abs(cor(baseline[, ..ph_pool], use = "pairwise.complete.obs"))
  diag(cor_mat) <- 0
}
cat("Phenotypes kept after |r| >", corr_cutoff, "filter:", length(ph_pool), "\n")

# ---------- QC grid ---------------------------------------------------------
cal_err_q <- quantile(dt$calErrAfter_mg,
                      probs = c(.90,.95,.975,.98), na.rm = TRUE)

qc_grid <- tidyr::crossing(
  require_goodWear = c(TRUE, FALSE),
  require_goodCal  = c(TRUE, FALSE),
  wearDays_min     = 3:7,
  calErrThr        = cal_err_q,
  clipsThr         = c(100, 500, 1000)
)

# ---------- evaluate grid ---------------------------------------------------
evaluate_qc <- function(row) {
  q <- as.list(row)
  fx <- dt[
    (goodWear == 1 | !q$require_goodWear) &
      (goodCal  == 1 | !q$require_goodCal) &
      wearDays >= q$wearDays_min &
      calErrAfter_mg < q$calErrThr &
      clipsAfter     < q$clipsThr
  ]
  n_ids <- uniqueN(fx$participant_id)
  if (n_ids < min_repeat_subs) return(NULL)

  iccs <- map_dbl(ph_pool, ~ icc_mixed(fx, .x))
  m_icc <- mean(iccs, na.rm = TRUE)

  out <- tibble(
    require_goodWear = q$require_goodWear,
    require_goodCal  = q$require_goodCal,
    wearDays_min     = q$wearDays_min,
    calErrThr        = q$calErrThr,
    clipsThr         = q$clipsThr,
    N_individuals    = n_ids,
    mean_ICC         = m_icc,
    eff_sample_size  = n_ids * m_icc
  )
  # add ICC_<phenotype> columns
  for (i in seq_along(ph_pool))
    out[[paste0("ICC_", ph_pool[i])]] <- iccs[i]
  out
}

results_raw <- purrr::map_dfr(seq_len(nrow(qc_grid)),
                              ~ evaluate_qc(qc_grid[.x, ]))

# ---------- z-scores and ordering ------------------------------------------
med_N   <- median(results_raw$N_individuals)
iqr_N   <- IQR(results_raw$N_individuals)
med_eN  <- median(results_raw$eff_sample_size)
iqr_eN  <- IQR(results_raw$eff_sample_size)

results <- results_raw %>%
  mutate(
    z_N    = (N_individuals   - med_N ) / iqr_N,
    z_effN = (eff_sample_size - med_eN) / iqr_eN
  ) %>%
  arrange(desc(eff_sample_size))

fwrite(results, OUT_CSV)
cat("\nBest QC rule (row 1 of CSV):\n")
print(results[1, 1:10], width = Inf)  # print first 10 cols for brevity
cat(sprintf("\nFull grid with %d phenotypes written to “%s”\n",
            length(ph_pool), OUT_CSV))
