#!/bin/bash

# This script runs step 2 of REGENIE on WES data. It first generates a list of
# high-quality SNPs for each chromosome and then runs the association analysis
# for each chromosome in parallel.

set -e

PROJECT_ID="project-J1GYbG0JQz0gxQ6yZf1GqYkb"
DX_PHENOTYPE_FILE_PATH="/temp_ukb_cohorts/gwas/pheno.tsv"
DX_COVARIATE_FILE_PATH="/temp_ukb_cohorts/gwas/covars.tsv"
DX_REGENIE_STEP1_DIR="/data/regenie_step1"
DX_WES_DATA_DIR="/Bulk/Exome sequences/Population level exome OQFE variants, PLINK format - final release/"
DX_WES_OUTPUT_DIR="/data/wes_gwas_results"
DATA_FIELD="ukb23158"

# The phenotype columns to be analyzed. This is a placeholder and should be
# replaced with the actual column names from your phenotype file.
PHENOTYPE_COLS="mesor,amp,acrophase,mesor_adj,amp_adj,acrophase_adj,mesor_amb,amp_amb,acrophase_amb,temp_range,rel_amp,act_slope,env_sens,resid_sd,fft_peak_freq_temp,fft_peak_amp_temp,fft_peak_freq_enmo,fft_peak_amp_enmo,r_squared_activity"

echo "Starting REGENIE step 2 for WES data..."

# First, generate the SNP lists for each chromosome
echo "Generating QC pass SNP lists..."
for i in {1..22} X; do
    run_plink_wes="plink2 --bfile ${DATA_FIELD}_c${i}_b0_v1\
      --no-pheno --keep pheno.tsv \
      --maf 0.0005 --mac 20 --geno 0.1 --mind 0.1\
      --write-snplist --write-samples --no-id-header\
      --out WES_c${i}_snps_qc_pass"

    dx run swiss-army-knife -iin="${DX_WES_DATA_DIR}/${DATA_FIELD}_c${i}_b0_v1.bed" \
     -iin="${DX_WES_DATA_DIR}/${DATA_FIELD}_c${i}_b0_v1.bim" \
     -iin="${DX_WES_DATA_DIR}/${DATA_FIELD}_c${i}_b0_v1.fam"\
     -iin="${DX_PHENOTYPE_FILE_PATH}" \
     -icmd="${run_plink_wes}" --tag="WES_QC_chr${i}" --instance-type "mem1_ssd1_v2_x16"\
     --destination="${PROJECT_ID}:${DX_WES_OUTPUT_DIR}/" --brief --yes
done

# Now, run REGENIE step 2 for each chromosome
echo "Running REGENIE step 2..."
for chr in {1..22} X; do
  run_regenie_cmd="regenie --step 2 --out assoc.c${chr} \
    --bed ${DATA_FIELD}_c${chr}_b0_v1 \
    --phenoFile pheno.tsv --covarFile covars.tsv \
    --bt --approx --firth-se --firth --extract WES_c${chr}_snps_qc_pass.snplist \
    --phenoColList ${PHENOTYPE_COLS} --covarColList age,sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,season_cos,season_sin \
    --pred ${DX_REGENIE_STEP1_DIR}/regenie_step1_results_pred.list --bsize 200 \
    --pThresh 0.05 --minMAC 3 --threads 16 --gz"

  dx run swiss-army-knife -iin="${DX_WES_DATA_DIR}/${DATA_FIELD}_c${chr}_b0_v1.bed" \
   -iin="${DX_WES_DATA_DIR}/${DATA_FIELD}_c${chr}_b0_v1.bim" \
   -iin="${DX_WES_DATA_DIR}/${DATA_FIELD}_c${chr}_b0_v1.fam" \
   -iin="${DX_WES_OUTPUT_DIR}/WES_c${chr}_snps_qc_pass.snplist" \
   -iin="${DX_PHENOTYPE_FILE_PATH}" \
   -iin="${DX_COVARIATE_FILE_PATH}" \
   -iin="${DX_REGENIE_STEP1_DIR}/regenie_step1_results_pred.list" \
   -iin="${DX_REGENIE_STEP1_DIR}/regenie_step1_results_1.loco.gz" \
   -icmd="${run_regenie_cmd}" --tag="Regenie_Step2_WES_chr${chr}" --instance-type "mem1_ssd1_v2_x16"\
--priority "high"
   --destination="${PROJECT_ID}:${DX_WES_OUTPUT_DIR}/" --brief --yes
done

echo "REGENIE step 2 WES script finished."