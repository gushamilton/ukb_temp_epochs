#!/bin/bash

# This script prunes the QC-ed genotype file for Linkage Disequilibrium (LD).

set -e

PROJECT_ID="project-J1GYbG0JQz0gxQ6yZf1GqYkb"
DX_PHENOTYPE_FILE_PATH="/temp_ukb_cohorts/gwas/covars.tsv"
DX_DATA_DIR="/data/gt_genrel_block"
DATA_FIELD="ukb22418"

echo "Step 1c: LD pruning..."
run_plink_prune="plink2 --bfile ${DATA_FIELD}_allQC_v2_merged \
 --indep-pairwise 1000 50 0.4  --out ukb-pruning; \
ls *bed; \
plink2 --bfile ${DATA_FIELD}_allQC_v2_merged --extract ukb-pruning.prune.in \
 --keep covars.tsv --make-bed --out ${DATA_FIELD}_allQC_v2_mrg_prun_cohort ; \
wc *.bim "

dx run swiss-army-knife -iin="${DX_DATA_DIR}/${DATA_FIELD}_allQC_v2_merged.bed" \
   -iin="${DX_DATA_DIR}/${DATA_FIELD}_allQC_v2_merged.bim" \
   -iin="${DX_DATA_DIR}/${DATA_FIELD}_allQC_v2_merged.fam" \
   -iin="${DX_PHENOTYPE_FILE_PATH}" \
   -icmd="${run_plink_prune}" --tag="Step1c_Prune" --instance-type "mem1_ssd1_v2_x16"
--priority "high" \
   --destination="${PROJECT_ID}:${DX_DATA_DIR}/" --brief --yes

echo "Step 1c finished."