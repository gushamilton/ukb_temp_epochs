#!/bin/bash

# This script converts the pruned PLINK file to VCF format.

set -e

PROJECT_ID="project-J1GYbG0JQz0gxQ6yZf1GqYkb"
DX_DATA_DIR="/data/gt_genrel_block"
DATA_FIELD="ukb22418"

echo "Step 1d: Convert to VCF..."
run_plink_tovcf="plink2 --bfile ${DATA_FIELD}_allQC_v2_mrg_prun_cohort --make-bpgen --out TEMP1 --merge-x --sort-vars n;\
plink2 --bpfile TEMP1 --export vcf bgz --out ukb_gt_p_temp ;\
rm TEMP1.* "

dx run swiss-army-knife -iin="${DX_DATA_DIR}/${DATA_FIELD}_allQC_v2_mrg_prun_cohort.bed" \
   -iin="${DX_DATA_DIR}/${DATA_FIELD}_allQC_v2_mrg_prun_cohort.bim" \
   -iin="${DX_DATA_DIR}/${DATA_FIELD}_allQC_v2_mrg_prun_cohort.fam"\
   -icmd="${run_plink_tovcf}" --tag="Step1d_ToVCF" --instance-type "mem1_ssd1_v2_x16"
--priority "high"

   --destination="${PROJECT_ID}:${DX_DATA_DIR}/" --brief --yes

echo "Step 1d finished."