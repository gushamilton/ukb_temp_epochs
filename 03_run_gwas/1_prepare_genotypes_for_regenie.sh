#!/bin/bash

# This script prepares the genotype data for REGENIE by running a series of
# sub-scripts on the DNAnexus platform. It merges, QCs, prunes, and liftovers the
# genotype data to GRCh38.

set -e

echo "Starting genotype preparation..."

bash 03_run_gwas/1a_merge_genotypes.sh
bash 03_run_gwas/1b_qc_filter_genotypes.sh
bash 03_run_gwas/1c_ld_prune_genotypes.sh
bash 03_run_gwas/1d_convert_to_vcf.sh
bash 03_run_gwas/1e_liftover_genotypes.sh
bash 03_run_gwas/1f_convert_from_vcf.sh

echo "Genotype preparation script finished."