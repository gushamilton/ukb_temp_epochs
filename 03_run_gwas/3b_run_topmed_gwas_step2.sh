#!/bin/bash

# This script runs step 2 of REGENIE on TOPMed data. It first generates a list of
# high-quality SNPs for each chromosome and then runs the association analysis
# for each chromosome in parallel.

set -e

PROJECT_ID="project-J1GYbG0JQz0gxQ6yZf1GqYkb"
DX_PHENOTYPE_FILE_PATH="/temp_ukb_cohorts/gwas/pheno.tsv"
DX_COVARIATE_FILE_PATH="/temp_ukb_cohorts/gwas/covars.tsv"
DX_REGENIE_STEP1_DIR="/data/regenie_step1"
DX_TOPMED_DATA_DIR="/Bulk/Imputation/Imputation from genotype (TOPmed)/"
DX_TOPMED_OUTPUT_DIR="/data/topmed_gwas_results"

# The phenotype columns to be analyzed. This is a placeholder and should be
# replaced with the actual column names from your phenotype file.
PHENOTYPE_COLS="mesor,amp,acrophase,mesor_adj,amp_adj,acrophase_adj,mesor_amb,amp_amb,acrophase_amb,temp_range,rel_amp,act_slope,env_sens,resid_sd,fft_peak_freq_temp,fft_peak_amp_temp,fft_peak_freq_enmo,fft_peak_amp_enmo,r_squared_activity"

echo "Starting REGENIE step 2 for TOPMed data..."

# First, generate the SNP lists for each chromosome
echo "Generating QC pass SNP lists..."
for i in {1..22} X; do
    run_plink_topmed="plink2 --bgen ukb21007_c${i}_b0_v1.bgen ref-first --sample ukb21007_c${i}_b0_v1.sample \
      --no-pheno --keep pheno.tsv \
      --maf 0.0005 --mac 20 --geno 0.1 --mind 0.1\ 
      --write-snplist --write-samples --no-id-header\ 
      --out TOPMED_c${i}_snps_qc_pass"

    dx run swiss-army-knife -iin="${DX_TOPMED_DATA_DIR}/ukb21007_c${i}_b0_v1.bgen" \
     -iin="${DX_TOPMED_DATA_DIR}/ukb21007_c${i}_b0_v1.sample"\ 
     -iin="${DX_PHENOTYPE_FILE_PATH}" \
     -icmd="${run_plink_topmed}" --tag="TOPMED_QC_chr${i}" --instance-type "mem1_ssd1_v2_x16"\
--priority "high"
     --destination="${PROJECT_ID}:${DX_TOPMED_OUTPUT_DIR}/" --brief --yes
done

# Now, run REGENIE step 2 for each chromosome
echo "Running REGENIE step 2..."
for chr in {1..22} X; do
  run_regenie_cmd="regenie --step 2 --out assoc.c${chr} \
    --bgen ukb21007_c${chr}_b0_v1.bgen --sample ukb21007_c${chr}_b0_v1.sample \
    --phenoFile pheno.tsv --covarFile covars.tsv \
    --bt --approx --extract TOPMED_c${chr}_snps_qc_pass.snplist \
    --phenoColList ${PHENOTYPE_COLS} --covarColList age,sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,season_cos,season_sin \
    --pred ${DX_REGENIE_STEP1_DIR}/regenie_step1_results_pred.list --bsize 200 \
    --pThresh 0.05 --minMAC 3 --threads 16 --gz"

  dx run swiss-army-knife -iin="${DX_TOPMED_DATA_DIR}/ukb21007_c${chr}_b0_v1.bgen" \
   -iin="${DX_TOPMED_DATA_DIR}/ukb21007_c${chr}_b0_v1.sample" \
   -iin="${DX_TOPMED_OUTPUT_DIR}/TOPMED_c${chr}_snps_qc_pass.snplist" \
   -iin="${DX_PHENOTYPE_FILE_PATH}" \
   -iin="${DX_COVARIATE_FILE_PATH}" \
   -iin="${DX_REGENIE_STEP1_DIR}/regenie_step1_results_pred.list" \
   -iin="${DX_REGENIE_STEP1_DIR}/regenie_step1_results_1.loco.gz" \
   -icmd="${run_regenie_cmd}" --tag="Regenie_Step2_TOPMED_chr${chr}" --instance-type "mem1_ssd1_v2_x16"\
--priority "high"
   --destination="${PROJECT_ID}:${DX_TOPMED_OUTPUT_DIR}/" --brief --yes
done

echo "REGENIE step 2 TOPMed script finished."