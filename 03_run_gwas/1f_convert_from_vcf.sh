#!/bin/bash

# This script converts the liftover VCF file back to PLINK format.

set -e

PROJECT_ID="project-J1GYbG0JQz0gxQ6yZf1GqYkb"
DX_DATA_DIR="/data/gt_genrel_block"
DATA_FIELD="ukb22418"

echo "Step 1f: Convert back to plink..."
run_plink_fromvcf="plink2 --vcf ukb_gt_lo38_sort.vcf.gz --make-bed --out ukb_gt_lo38_temp --id-delim \"_\" \
  --indiv-sort f ${DATA_FIELD}_allQC_v2_mrg_prun_cohort.fam --ref-allele force ukb_gt_lo38_sort.vcf.gz 4 3 '#' \
  --alt1-allele force ukb_gt_lo38_sort.vcf.gz 5 3 '#' --allow-extra-chr ;\
plink2 --bfile ukb_gt_lo38_temp --make-bed --out ${DATA_FIELD}_allQC_mrg_prnd_cohort_lo38_v2 \
  --fam ${DATA_FIELD}_allQC_v2_mrg_prun_cohort.fam --allow-extra-chr ;\
rm ukb_gt_lo38_temp.* ;\
grep alt ukb22418_allQC_mrg_prnd_cohort_lo38_v2.bim | awk '{print $2}' > remove_losnps.txt; \
plink2 --bfile ${DATA_FIELD}_allQC_mrg_prnd_cohort_lo38_v2 --exclude remove_losnps.txt --make-bed \
  -out ${DATA_FIELD}_allQC_mrg_prnd_cohort_lo38_v3 --allow-extra-chr "

dx run swiss-army-knife -iin="${DX_DATA_DIR}/ukb_gt_lo38_sort.vcf.gz" \
   -iin="${DX_DATA_DIR}/${DATA_FIELD}_allQC_v2_mrg_prun_cohort.fam"\
   -icmd="${run_plink_fromvcf}" --tag="Step1f_FromVCF" --instance-type "mem2_ssd1_v2_x16"\n--priority "high"\
   --destination="${PROJECT_ID}:${DX_DATA_DIR}/" --brief --yes

echo "Step 1f finished."