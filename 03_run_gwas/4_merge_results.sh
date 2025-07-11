#!/bin/bash

# This script merges the per-chromosome results from the REGENIE analysis into
# a single file for each GWAS (WES and TOPMed).

set -e

PROJECT_ID="project-J1GYbG0JQz0gxQ6yZf1GqYkb"
DX_WES_OUTPUT_DIR="/data/wes_gwas_results"
DX_TOPMED_OUTPUT_DIR="/data/topmed_gwas_results"

echo "Starting results merging..."

# Merge WES results
echo "Merging WES results..."
merge_wes_cmd='out_file="assoc.regenie.merged.txt"\
     cp /mnt/project${DX_WES_OUTPUT_DIR}/*.regenie.gz .\
     gunzip *.regenie.gz\
\
echo -e "CHROM\tGENPOS\tID\tALLELE0\tALLELE1\tA1FREQ\tN\tTEST\tBETA\tSE\tCHISQ\tLOG10P\tEXTRA" > $out_file\
\
files="./*.regenie"\
for f in $files\
do\
   tail -n+2 $f | tr " " "\t" >> $out_file\
done\
\
rm *.regenie'

dx run swiss-army-knife -icmd="${merge_wes_cmd}" --tag="Merge_WES" --instance-type "mem1_ssd1_v2_x16"\
--priority "high"\
   --destination="${PROJECT_ID}:${DX_WES_OUTPUT_DIR}/" --brief --yes

# Merge TOPMed results
echo "Merging TOPMed results..."
merge_topmed_cmd='out_file="assoc.regenie.merged.txt"\
     cp /mnt/project${DX_TOPMED_OUTPUT_DIR}/*.regenie.gz .\
     gunzip *.regenie.gz\
\
echo -e "CHROM\tGENPOS\tID\tALLELE0\tALLELE1\tA1FREQ\tN\tTEST\tBETA\tSE\tCHISQ\tLOG10P\tEXTRA" > $out_file\
\
files="./*.regenie"\
for f in $files\
do\
   tail -n+2 $f | tr " " "\t" >> $out_file\
done\
\
rm *.regenie'

dx run swiss-army-knife -icmd="${merge_topmed_cmd}" --tag="Merge_TOPMed" --instance-type "mem1_ssd1_v2_x16"\
--priority "high"\
   --destination="${PROJECT_ID}:${DX_TOPMED_OUTPUT_DIR}/" --brief --yes

echo "Results merging finished."
