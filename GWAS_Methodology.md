# Genome-Wide Association Study (GWAS) Methodology

This section details the comprehensive pipeline employed for conducting Genome-Wide Association Studies (GWAS) on derived human body temperature phenotypes within the UK Biobank cohort. The methodology encompasses phenotype derivation, rigorous genotype preparation, and association testing using state-of-the-art statistical software.

## 1. Sample Definition, Phenotype Derivation, and Quality Control

This study focused on participants identified as **Caucasian** based on the UK Biobank field `p22006` ("Ethnic background (inferred)"). Raw, high-frequency wrist-worn accelerometer data from these participants were then processed to derive stable thermoregulatory phenotypes. This process involved several key steps:

*   **Data Extraction and Consolidation:** Raw `.cwa` files were processed to extract 30-second epoch-level temperature and activity data. These batched outputs were then consolidated into a single dataset.
*   **Repeat Measure Handling:** For participants with multiple measurement periods, phenotypes were averaged across all available records to provide a single, representative value per individual. This approach aims to capture the stable, long-term physiological characteristics.
*   **Phenotype Definitions:** The following derived temperature phenotypes were utilized:
    *   `mesor`: Rhythm- and environment-adjusted mean temperature.
    *   `amp`: Amplitude of the circadian rhythm.
    *   `acrophase`: Timing of the circadian rhythm's peak.
    *   `mesor_adj`, `amp_adj`, `acrophase_adj`: Corresponding parameters adjusted for activity.
    *   `mesor_amb`, `amp_amb`, `acrophase_amb`: Corresponding parameters adjusted for ambient temperature.
    *   `temp_range`: Range between 5th and 95th percentile temperature.
    *   `rel_amp`: Relative amplitude.
    *   `act_slope`: Slope of temperature change with activity.
    *   `env_sens`: Sensitivity to ambient temperature.
    *   `resid_sd`: Standard deviation of model residuals (temperature lability).
    *   `fft_peak_freq_temp`, `fft_peak_amp_temp`: Peak frequency and amplitude from Fast Fourier Transform of temperature.
    *   `fft_peak_freq_enmo`, `fft_peak_amp_enmo`: Peak frequency and amplitude from Fast Fourier Transform of ENMO (activity).
    *   `r_squared_activity`: R-squared from temperature vs. activity regression.
*   **Quality Control (QC):** Rigorous QC filters were applied to the derived phenotypes to ensure data reliability. Participants were included only if their data met predefined thresholds for:
    *   Minimum wear time (e.g., >X days).
    *   Maximum calibration error (e.g., <Y mg).
    *   Maximum number of data clips after calibration (e.g., <Z clips).
    *   Requirement for good wear time and good calibration (binary flags).
*   **Software:** Phenotype derivation and QC were performed using **R** (version X.X.X) with packages including `tidyverse` (version X.X.X), `data.table` (version X.X.X), and `lubridate` (version X.X.X).

## 2. Genotype Data Preparation

UK Biobank genotype data, initially in GRCh37 build, underwent a multi-step preparation process to ensure suitability for GWAS analysis with REGENIE. This preparation was performed exclusively on the **Caucasian subset** of participants.

*   **Data Merging:** Per-chromosome genotype files were merged into a single dataset.
*   **Genotype QC:** Standard quality control filters were applied using **PLINK2** (version X.X):
    *   Minor Allele Frequency (MAF) > 0.01
    *   Minor Allele Count (MAC) > 20
    *   Genotype missingness per SNP < 0.1
    *   Individual missingness < 0.1
*   **Linkage Disequilibrium (LD) Pruning:** The QC-filtered genotype data was pruned for LD (e.g., `indep-pairwise 1000 50 0.4`) using **PLINK2** to generate a set of approximately independent markers. This pruned set is used by REGENIE's Step 1 to build the whole-genome regression model.
*   **Liftover to GRCh38:** The genotype data was lifted over from GRCh37 to GRCh38 using **Picard** (version X.X.X) and **bcftools** (version X.X.X). This involved converting PLINK format to VCF, performing the liftover, and converting back to PLINK format.

## 3. Genome-Wide Association Analysis with REGENIE

Association testing was performed using **REGENIE** (version X.X.X), a powerful and scalable software designed for large-scale biobank data, which inherently accounts for relatedness and population structure.

*   **REGENIE Step 1 (Whole-Genome Regression):** This step involved building a whole-genome regression model using the LD-pruned genotype data. This model accounts for cryptic relatedness and population structure via a Leave-One-Chromosome-Out (LOOCV) approach. The output of this step is a set of predicted phenotypes for each individual, which are then used in Step 2.
*   **REGENIE Step 2 (Association Testing):** In this step, association tests were performed for each phenotype against:
    *   **Whole Exome Sequencing (WES) Data:** Using the UK Biobank's exome sequence data (GRCh38).
    *   **TOPMed Imputed Data:** Using the UK Biobank's TOPMed imputed genotype data (GRCh38), which provides broader genomic coverage.
*   **Covariates:** All models included the following covariates to adjust for confounding factors:
    *   Age
    *   Sex (binary: Male/Female)
    *   First 10 Genetic Principal Components (PC1-PC10) to account for population structure.
    *   Season of assessment, modeled using sine and cosine terms (e.g., `sin(2*pi*day_of_year/365.25)` and `cos(2*pi*day_of_year/365.25)`) to capture the cyclical nature of seasonal variation.
*   **QC Filters (REGENIE):** Within REGENIE, additional variant-level QC was applied (e.g., `minMAC 3`, `pThresh 0.05`).

## 4. Results Merging

Following the per-chromosome association tests, the results were merged into single, comprehensive summary files for both WES and TOPMed GWAS. This involved concatenating the per-chromosome output files and ensuring consistent formatting.

---

**Note:** Please replace placeholder values (e.g., `X.X.X`, `X days`, `Y mg`, `Z clips`) with the specific versions and thresholds used in your analysis.
