#!/bin/bash

# This script runs the R script to create the GWAS input files and then uploads
# them to the specified location on the DNAnexus platform.

set -e

DX_PROJECT_ID="project-J1GYbG0JQz0gxQ6yZf1GqYkb"
DX_PHENO_COVAR_DIR="/temp_ukb_cohorts/gwas"

echo "Running R script to create GWAS inputs..."
Rscript 03_run_gwas/create_gwas_inputs.R
echo "Finished creating local GWAS input files."

echo "Uploading phenotype and covariate files to DNAnexus..."

# Create the destination directory on DNAnexus if it doesn't exist
dx mkdir -p ${DX_PROJECT_ID}:${DX_PHENO_COVAR_DIR}

# Upload the files
dx upload pheno.tsv --path ${DX_PROJECT_ID}:${DX_PHENO_COVAR_DIR}/pheno.tsv
dx upload covars.tsv --path ${DX_PROJECT_ID}:${DX_PHENO_COVAR_DIR}/covars.tsv

echo "✔ Successfully uploaded files to ${DX_PHENO_COVAR_DIR}"
