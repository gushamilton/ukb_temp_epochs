#!/bin/bash

# This script merges per-chromosome genotype files into a single dataset.

set -e

PROJECT_ID="project-J1GYbG0JQz0gxQ6yZf1GqYkb"
DX_PHENOTYPE_FILE_PATH="/temp_ukb_cohorts/gwas/pheno.tsv"
DX_DATA_DIR="/data/gt_genrel_block"
DATA_FIELD="ukb22418"

echo "Step 1a: Merging genotype files..."
run_merge="cp /mnt/project/Bulk/Genotype\ Results/Genotype\ calls/ukb22418_c[1-9]*_b0_v2.{bed,bim,fam} . ; cp /mnt/project/Bulk/Genotype\ Results/Genotype\ calls/ukb22418_cX_b0_v2.{bed,bim,fam} . ; ls *.bed | sed -e 's/.bed//g'> files_to_merge.txt; plink --merge-list files_to_merge.txt --make-bed --out ukb22418_all_v2_merged; rm files_to_merge.txt; rm ukb22418_c*"

dx run swiss-army-knife -iin="${DX_PHENOTYPE_FILE_PATH}" \
   -icmd="${run_merge}" --tag="Step1a_Merge" --instance-type "mem1_ssd1_v2_x16" --priority "high" \
   --destination="${PROJECT_ID}:${DX_DATA_DIR}/" --brief --yes

echo "Step 1a finished."