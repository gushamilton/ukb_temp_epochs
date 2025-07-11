#!/bin/bash

# This script applies QC filters to the merged genotype file.

set -e

PROJECT_ID="project-J1GYbG0JQz0gxQ6yZf1GqYkb"
DX_DATA_DIR="/data/gt_genrel_block"
DATA_FIELD="ukb22418"

echo "Step 1b: QC filtering..."
run_plink_qc="plink2 --bfile ${DATA_FIELD}_all_v2_merged \
 --maf 0.01 --mac 20 --geno 0.1 \
 --mind 0.1 --make-bed \
 --out  ${DATA_FIELD}_allQC_v2_merged"

dx run swiss-army-knife -iin="${DX_DATA_DIR}/${DATA_FIELD}_all_v2_merged.bed" \
   -iin="${DX_DATA_DIR}/${DATA_FIELD}_all_v2_merged.bim" \
   -iin="${DX_DATA_DIR}/${DATA_FIELD}_all_v2_merged.fam" \
   -icmd="${run_plink_qc}" --tag="Step1b_QC" --instance-type "mem1_ssd1_v2_x16"
--priority "high" \
   --destination="${PROJECT_ID}:${DX_DATA_DIR}/" --brief --yes

echo "Step 1b finished."