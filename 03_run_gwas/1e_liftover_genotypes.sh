#!/bin/bash

# This script performs liftover of genotype data from GRCh37 to GRCh38.

set -e

PROJECT_ID="project-J1GYbG0JQz0gxQ6yZf1GqYkb"
DX_DATA_DIR="/data/gt_genrel_block"

echo "Step 1e: Liftover to GRCh38..."
run_liftover="wget https://github.com/broadinstitute/picard/releases/download/2.27.4/picard.jar ;\
wget https://raw.githubusercontent.com/broadinstitute/gatk/master/scripts/funcotator/data_sources/gnomAD/b37ToHg38.over.chain ;\
cp /mnt/project/Bulk/Exome\ sequences/Exome\ OQFE\ CRAM\ files/helper_files/GRCh38_full_analysis_set_plus_decoy_hla.* . ;\
java -jar picard.jar LiftoverVcf -I ukb_gt_p_temp.vcf.gz -O ukb_gt_lo38.vcf \n   -C b37ToHg38.over.chain --REJECT rejected_variants.vcf -R GRCh38_full_analysis_set_plus_decoy_hla.fa \n   --RECOVER_SWAPPED_REF_ALT true --DISABLE_SORT true ;\
bcftools sort -o ukb_gt_lo38_sort.vcf.gz -O z ukb_gt_lo38.vcf"

dx run swiss-army-knife -iin="${DX_DATA_DIR}/ukb_gt_p_temp.vcf.gz" \
   -icmd="${run_liftover}" --tag="Step1e_Liftover" --instance-type "mem2_ssd1_v2_x16"\n--priority "high"\
   --destination="${PROJECT_ID}:${DX_DATA_DIR}/" --brief --yes

echo "Step 1e finished."