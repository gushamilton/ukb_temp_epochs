#!/bin/bash

# This script runs step 1 of REGENIE, which builds the whole-genome regression
# model.

set -e

PROJECT_ID="project-J1GYbG0JQz0gxQ6yZf1GqYkb"
DX_GENOTYPE_DIR="/data/gt_genrel_block"
DX_PHENOTYPE_FILE_PATH="/temp_ukb_cohorts/gwas/pheno.tsv"
DX_COVARIATE_FILE_PATH="/temp_ukb_cohorts/gwas/covars.tsv"
DX_OUTPUT_DIR="/data/regenie_step1"
DATA_FIELD="ukb22418"

# The phenotype columns to be analyzed. This is a placeholder and should be
# replaced with the actual column names from your phenotype file.
PHENOTYPE_COLS="mesor,amp,acrophase,mesor_adj,amp_adj,acrophase_adj,mesor_amb,amp_amb,acrophase_amb,temp_range,rel_amp,act_slope,env_sens,resid_sd,fft_peak_freq_temp,fft_peak_amp_temp,fft_peak_freq_enmo,fft_peak_amp_enmo,r_squared_activity"

echo "Starting REGENIE step 1..."

run_regenie_cmd="regenie --step 1 --out regenie_step1_results \
 --bed ${DATA_FIELD}_allQC_mrg_prnd_cohort_lo38_v3 \
 --phenoFile pheno.tsv \
 --covarFile covars.tsv \
 --phenoColList ${PHENOTYPE_COLS} \
 --covarColList age,sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,season_cos,season_sin \
 --bsize 1000 --bt --loocv --gz --threads 16 "

dx run swiss-army-knife -iin="${DX_GENOTYPE_DIR}/${DATA_FIELD}_allQC_mrg_prnd_cohort_lo38_v3.bed" \
   -iin="${DX_GENOTYPE_DIR}/${DATA_FIELD}_allQC_mrg_prnd_cohort_lo38_v3.bim" \
   -iin="${DX_GENOTYPE_DIR}/${DATA_FIELD}_allQC_mrg_prnd_cohort_lo38_v3.fam"\
   -iin="${DX_PHENOTYPE_FILE_PATH}" \
   -iin="${DX_COVARIATE_FILE_PATH}" \
   -icmd="${run_regenie_cmd}" --tag="Regenie_Step1" --instance-type "mem1_ssd1_v2_x16" \
--priority "high" \
   --destination="${PROJECT_ID}:${DX_OUTPUT_DIR}/" --brief --yes;

echo "REGENIE step 1 script finished."