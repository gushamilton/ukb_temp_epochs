#!/usr/bin/env Rscript
## optimise_qc.R -------------------------------------------------------------
## Single-core grid search that maximises   N × mean(ICC)   across QC rules
## and discards phenotypes whose |r| ≥ 0.90 in the baseline subset.
## --------------------------------------------------------------------------

suppressPackageStartupMessages({
  pkgs <- c("data.table", "tidyverse", "lme4")
  for (p in pkgs) if (!requireNamespace(p, quietly = TRUE))
      stop("Package ‘", p, "’ must be installed.")
  library(data.table)
  library(tidyverse)   # dplyr, tidyr, purrr
  library(lme4)        # lmer for NA-robust ICC
})

# ---------- paths -----------------------------------------------------------
COVAR_FILE  <- "combined_covariates.csv.gz"
PHENO_DIR   <- "phenotype_blocks"
OUT_CSV     <- "qc_grid_results.csv"

# ---------- load data -------------------------------------------------------
cat("→ reading covariates …\n")
covar <- fread(
  COVAR_FILE,
  select = c("id","participant_id",
             "goodWear","goodCal","calErrAfter_mg",
             "wearDays","clipsAfter")
)

cat("→ reading phenotype blocks …\n")
ph_files <- list.files(PHENO_DIR, "^phenotypes_.*\\.tsv\\.gz$", full.names = TRUE)
if (!length(ph_files)) stop("No phenotype files found in ‘phenotype_blocks/’.")
phenos <- rbindlist(lapply(ph_files, fread), use.names = TRUE, fill = TRUE)

dt <- merge(phenos, covar, by = c("id","participant_id"))
rm(phenos, covar); gc()

# ---------- helper: ICC via mixed model -------------------------------------
icc_mixed <- function(data, ph) {
  d <- data[!is.na(get(ph)), .(participant_id, y = get(ph))]
  reps <- d[, .N, by = participant_id][N > 1, participant_id]
  if (length(reps) < 10) return(NA_real_)
  d <- d[participant_id %in% reps]

  val <- tryCatch({
    m  <- lmer(y ~ 1 + (1|participant_id), data = d,
               control = lmerControl(check.conv.singular = "ignore"))
    vc <- as.data.frame(VarCorr(m))
    vb <- vc$vcov[vc$grp == "participant_id"]
    ve <- attr(VarCorr(m), "sc")^2
    vb / (vb + ve)
  }, error = function(e) NA_real_)
  as.numeric(val)
}

# ---------- 3 a. build loose baseline ---------------------------------------
cat("→ building loose-QC baseline …\n")
baseline <- dt[ wearDays >= 3 & clipsAfter < 1000 ]
n_rep <- baseline[, .N, by = participant_id][N > 1, .N]
if (n_rep < 20) stop("Too few repeat participants in baseline (", n_rep, ").")

num_cols <- names(baseline)[vapply(baseline, is.numeric, FALSE)]
drop_qc  <- c("id","participant_id",
              "goodWear","goodCal","calErrAfter_mg",
              "wearDays","clipsAfter")
pheno_pool <- setdiff(num_cols, drop_qc)

# ---------- 3 b. de-duplicate highly-correlated phenotypes ------------------
cat("→ removing highly-correlated phenotypes (|r| ≥ 0.90) …\n")
cor_mat <- abs(cor(baseline[, ..pheno_pool], use = "pairwise.complete.obs"))
diag(cor_mat) <- 0

while (TRUE) {
  mx <- which(cor_mat > 0.90, arr.ind = TRUE)
  if (nrow(mx) == 0) break
  # always drop the *second* column of the first offending pair
  to_drop <- colnames(cor_mat)[mx[1, 2]]
  pheno_pool <- setdiff(pheno_pool, to_drop)
  cor_mat <- abs(cor(baseline[, ..pheno_pool], use = "pairwise.complete.obs"))
  diag(cor_mat) <- 0
}
cat("   remaining phenotypes after pruning :", length(pheno_pool), "\n")

# ---------- 3 c. pick top-5 by ICC ------------------------------------------
cat("→ computing ICCs for baseline …\n")
baseline_icc <- map_dbl(pheno_pool, ~ icc_mixed(baseline, .x))
top5 <- pheno_pool[order(baseline_icc, decreasing = TRUE)][1:5]
if (anyNA(top5) || length(top5) < 5)
  stop("Unable to identify five phenotypes with finite ICCs.")
cat("   top-5 phenotypes :", paste(top5, collapse = ", "), "\n\n")

# ---------- 4. QC grid -------------------------------------------------------
cal_err_q <- quantile(dt$calErrAfter_mg,
                      probs = c(.90,.95,.975,.98), na.rm = TRUE)
qc_grid <- crossing(
  require_goodWear = c(TRUE, FALSE),
  require_goodCal  = c(TRUE, FALSE),
  wearDays_min     = 3:7,
  calErrThr        = cal_err_q,
  clipsThr         = c(100, 500, 1000)
)

# ---------- 5. grid evaluation ----------------------------------------------
cat(sprintf("→ evaluating %d QC combinations …\n", nrow(qc_grid)))
evaluate_one <- function(row) {
  q <- as.list(row)
  fx <- dt[
    (goodWear == 1 | !q$require_goodWear) &
      (goodCal  == 1 | !q$require_goodCal) &
      wearDays >= q$wearDays_min           &
      calErrAfter_mg < q$calErrThr         &
      clipsAfter       < q$clipsThr
  ]
  n_ids <- length(unique(fx$participant_id))
  if (n_ids < 20) return(NULL)

  iccs  <- map_dbl(top5, ~ icc_mixed(fx, .x))
  m_icc <- mean(iccs, na.rm = TRUE)

  tibble(
    require_goodWear = q$require_goodWear,
    require_goodCal  = q$require_goodCal,
    wearDays_min     = q$wearDays_min,
    calErrThr        = q$calErrThr,
    clipsThr         = q$clipsThr,
    N_individuals    = n_ids,
    mean_ICC         = m_icc,
    eff_sample_size  = n_ids * m_icc
  )
}

results <- purrr::map_dfr(1:nrow(qc_grid),
                          ~ evaluate_one(qc_grid[.x, ]))

# ---------- 6. write & show --------------------------------------------------
fwrite(results, OUT_CSV)
best <- results %>% arrange(desc(eff_sample_size)) %>% slice(1)

cat("\n================  best QC rule  ================\n")
print(best, width = Inf)
cat(sprintf("\nfull grid written to “%s”\n", OUT_CSV))
