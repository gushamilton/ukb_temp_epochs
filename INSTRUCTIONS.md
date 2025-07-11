# End-to-End GWAS Pipeline Instructions

This document provides a step-by-step guide to running the full REGENIE-based GWAS pipeline for both WES and TOPMed data.

---

### **1. Prerequisites**

Before you begin, please ensure you have the following set up:

1.  **DNAnexus CLI (`dx-toolkit`):** You must have the `dx-toolkit` installed on your local machine.
2.  **Logged In:** You must be logged into your DNAnexus account. You can do this by running `dx login`.
3.  **Project Selected:** You must have your UK Biobank project selected. Run `dx select` and choose the project `project-J1GYbG0JQz0gxQ6yZf1GqYkb` if it is not already selected.
4.  **Working Directory:** All commands in this guide must be run from your local project directory: `/Users/fh6520/R/temp_ukb`.
5.  **Participant Data:** The file `data_participant.tsv` must be present in your working directory (`/Users/fh6520/R/temp_ukb`).
6.  **R Environment:** You need a working installation of R on your local machine. The R script will attempt to install its own required packages (`tidyverse`, `data.table`, `lubridate`) if they are not already present.

---

### **2. Execution Workflow**

Follow these steps in order. Each `bash` command executes a script on your local machine, which in turn submits jobs to the DNAnexus platform.

#### **Step 0: Create and Upload Phenotype and Covariate Files**

This script first runs an R script on your local machine to process the raw phenotype batches and create the final `pheno.tsv` and `covars.tsv` files. It then automatically uploads these two files to the `/temp_ukb_cohorts/gwas/` directory in your DNAnexus project, making them ready for the GWAS pipeline.

*   **Where to run:** Your local terminal.
*   **Command:**
    ```bash
    bash 03_run_gwas/0_create_gwas_inputs.sh
    ```
*   **What it does:**
    1.  Downloads multiple phenotype batches from DNAnexus to your local machine.
    2.  Runs an R script locally to generate `pheno.tsv` and `covars.tsv`.
    3.  Uploads the final `pheno.tsv` and `covars.tsv` to DNAnexus.

#### **Step 1: Prepare Genotypes for REGENIE**

This script submits a series of jobs to DNAnexus to merge, QC, prune, and liftover the genotype data to prepare it for REGENIE.

*   **Where to run:** Your local terminal.
*   **Command:**
    ```bash
    bash 03_run_gwas/1_prepare_genotypes_for_regenie.sh
    ```
*   **What to expect:** This will submit 6 jobs to DNAnexus. You can monitor their progress in the DNAnexus web interface.

#### **Step 2: Run REGENIE Step 1**

This script runs the first major step of REGENIE, which builds the whole-genome regression model. This is a computationally intensive step.

*   **Where to run:** Your local terminal.
*   **Command:**
    ```bash
    bash 03_run_gwas/2_run_regenie_step1.sh
    ```
*   **What to expect:** This submits a single, long-running job to DNAnexus.

#### **Step 3: Run REGENIE Step 2 (Association Testing)**

These two scripts run the association tests for the WES and TOPMed data, respectively. They will submit one job per chromosome.

*   **Where to run:** Your local terminal.
*   **Commands:**
    ```bash
    # Run the WES GWAS
    bash 03_run_gwas/3a_run_wes_gwas_step2.sh

    # Run the TOPMed GWAS
    bash 03_run_gwas/3b_run_topmed_gwas_step2.sh
    ```
*   **What to expect:** Each script will submit many jobs to DNAnexus (one for each chromosome). This is the main GWAS analysis.

#### **Step 4: Merge Final Results**

This final script merges the per-chromosome results from both the WES and TOPMed analyses into single, easy-to-use summary files.

*   **Where to run:** Your local terminal.
*   **Command:**
    ```bash
    bash 03_run_gwas/4_merge_results.sh
    ```
*   **What to expect:** This submits two jobs to DNAnexus, one for each GWAS run. The final merged files will be named `assoc.regenie.merged.txt` and will be located in the `/data/wes_gwas_results/` and `/data/topmed_gwas_results/` directories in your project.

---

### **3. A Note on Spot vs. On-Demand Instances**

You asked if these scripts are continuable, and the answer is **no, they are not automatically resumable.** Each `dx run` command submits an independent job to the platform. If a job is terminated (e.g., because a spot instance is reclaimed), it will fail and will **not** restart on its own.

*   **REGENIE (Steps 2 & 3)** and the **Genotype Preparation (Step 1)** are long-running, computationally expensive jobs.
*   If one of these jobs fails due to a spot instance interruption, you will lose all progress for that job and will have to rerun it from the beginning, wasting significant time and money.

**Recommendation:** For reliability, it is **strongly recommended** that you use **on-demand instances** for this pipeline. The scripts are currently configured to use standard on-demand instance types. You should not change this unless you are willing to accept the risk of job failures and the need for manual intervention.
