# A Comprehensive Research Plan for the GWAS of Human Body Temperature
**Version: 2.0**
**Objective:** To define a robust, end-to-end pipeline for deriving stable, heritable, and biologically meaningful thermoregulatory phenotypes from noisy, high-frequency wrist-worn sensor data for use in a Genome-Wide Association Study (GWAS).

---

### **Part 1: The Biological and Genetic Foundation of Human Thermoregulation**

A comprehensive investigation into the genetic architecture of human body temperature requires a foundational understanding of the complex, multi-system physiology that governs thermal homeostasis. Human body temperature is not a static parameter but a dynamic phenotype resulting from the interplay between basal metabolic processes, endogenous circadian rhythms, and reactive responses to internal and external stimuli. A successful GWAS must be designed to dissect these constituent components, as they are likely governed by distinct genetic pathways.

The regulation of body temperature is orchestrated by a sophisticated control system centered in the preoptic area (POA) of the anterior hypothalamus. Superimposed on this homeostatic regulation is a robust, endogenous circadian rhythm, driven by a core set of "clock genes" (*CLOCK*, *BMAL1*, *PER*, *CRY*). This direct link provides a strong a priori hypothesis that genetic variants within this pathway are prime candidates for influencing circadian temperature parameters.

Beyond the core circadian machinery, specific genetic pathways are known to modulate the body's response to thermal challenges and inflammation. These include:
* **Thermogenesis & Energy Balance:** Uncoupling proteins (*UCP1–5*) in mitochondria.
* **Cold & Heat Stress Response:** Genes like *ACTN3* and those encoding heat shock proteins (HSPs).
* **Pyrogenic (Fever) Pathways:** Cytokines (*IL-1*, *IL-6*), their receptors, and genes involved in hereditary recurrent fever syndromes (*MEFV*, *TNFRSF1A*).

A naively defined phenotype, such as a simple average of all temperature readings, would conflate genetic signals from these distinct systems. Therefore, a primary objective of this research plan is to leverage advanced modeling techniques to derive multiple, specific phenotypes that can disentangle these components, allowing for a more nuanced and interpretable GWAS.

---

### **Part 2: Initial Data Extraction from Raw `.cwa` Files**

This initial stage converts the raw, proprietary binary accelerometer files (`.cwa`) from the UK Biobank into a structured, analysis-ready format. The provided bash script automates this process for large-scale execution in a cloud environment like DNAnexus.

#### **2.1. Process Overview**

The script parallelizes the processing of a master list of `.cwa` files. For each file, it executes a worker function that performs the core conversion using the `accelerometer` package.

#### **2.2. Key Operations & Pseudocode**

The central logic for processing a single participant file can be summarized with the following pseudocode:

FUNCTION process_one_file(raw_file_path):
// Identify participant ID from filename
id = get_id_from(raw_file_path)

// Use the 'accProcess' tool to convert the binary .cwa file.
// This generates detailed epoch-level data and a summary file.
accProcess(
input = raw_file_path,
output_folder = temporary_directory
)

// --- Process the generated epoch data ---
// Load the detailed time-series CSV
epoch_data = read_csv(f"{id}-epoch.csv.gz")

// Select only the necessary columns for our analysis
relevant_columns = ["time", "temp", "enmoTrunc", "dataErrors", "samples", "clipsBeforeCalibr", "clipsAfterCalibr"]
filtered_epoch_data = epoch_data[relevant_columns]

// Rename columns for clarity
filtered_epoch_data.rename(columns={"clipsBeforeCalibr": "clipsBefore", "clipsAfterCalibr": "clipsAfter"})

// Add participant ID to each row
filtered_epoch_data.insert_column(0, "id", id)

// Save the cleaned time-series data in an efficient format
filtered_epoch_data.to_parquet("epoch.part.parquet")

// --- Process the generated summary data ---
// Load the summary JSON file containing QC metrics
summary_json = load_json(f"{id}-summary.json")

// Extract key summary and QC metrics into a single row
summary_row = {
id: id,
goodWearTime: summary_json["quality-goodWearTime"],
goodCalibration: summary_json["quality-goodCalibration"],
wearTimeDays: summary_json["wearTime-overall(days)"],
// ... and other relevant summary statistics
}

// Save the summary row as a CSV
summary_row.to_csv("summary.part.csv")

END FUNCTION


#### **2.3. Final Outputs of this Stage**
This script aggregates the results from all participants into two key files:
1.  **`temperature_epochs.parquet`**: A single large table containing the time-series data (temperature, activity, etc.) in 30-second epochs for all participants.
2.  **`summary_covariates.csv.gz`**: A table where each row corresponds to a participant, containing crucial metadata and QC metrics (e.g., total wear time, calibration quality).

---

### **Part 3: A Rigorous Pipeline for Phenotype Derivation from Noisy Sensor Data**

The translation of raw, high-frequency, and noisy wrist-worn sensor data into stable, reliable phenotypes is the most critical and challenging phase of this study.

#### **Step 3.1: Data Loading and Quality Control Filtering**
1.  **Load Data:** Load the `temperature_epochs.parquet` and `summary_covariates.csv` files.
2.  **Merge:** Join the two tables on the participant `id`.
3.  **Initial QC:** Filter out entire participants based on the summary metrics. Define sensible thresholds to ensure data quality.
    * `wearTimeDays > 5` (Ensure sufficient duration)
    * `goodCalibration == TRUE` (Ensure sensor calibration was successful)
    * `calErrAfter_mg < 10` (Ensure calibration error is low)

#### **Step 3.2: Time-Series Cleaning and Downsampling**
The raw 30-second temperature data is extremely noisy. The goal is to smooth this data and reduce its frequency to make subsequent modeling computationally feasible and more robust.

* **Challenge:** The data contains non-physiological spikes from environmental contact (e.g., washing hands, holding a hot drink). A simple rolling average is highly sensitive to these outliers.
* **Recommended Approach (Balanced): Rolling Median & Mean**
    1.  **Outlier Removal:** For each participant, apply a rolling **median** filter to the 30-second `temp` data. A window of 5-15 minutes is appropriate. This effectively removes sharp spikes without distorting the underlying signal.
    2.  **Downsampling:** Aggregate the cleaned 30-second data into **5-minute epochs**. For each 5-minute window, calculate the **mean** of the cleaned temperature and the **mean** of the activity metric (`enmoTrunc`).
* **Alternative Approach (Higher Fidelity, Higher Compute): Smoothing Spline**
    1.  For each participant, fit a smoothing spline to the full 30-second temperature data.
    2.  Evaluate the fitted spline at regular 5-minute intervals.
    * **Comparison:** A spline more accurately represents the continuous physiological curve but is substantially more computationally expensive than the rolling window approach. For a dataset of this scale, the rolling median/mean approach offers a robust and practical trade-off.

#### **Step 3.3: Environmental & Behavioural Data Integration**
1.  **Ambient Temperature:** Merge the time-series data with external daily weather data (e.g., from Met Office) based on the date of each measurement. This is a critical covariate.
2.  **Sleep Period Identification:** Use the downsampled 5-minute `enmoTrunc` (activity) data to define sleep windows. A common algorithm is to find periods of at least 2-3 hours where activity remains below a low threshold. This creates a "sleep-only" subset of the temperature data, which is invaluable for isolating the core biological rhythm from daytime behavioural noise.

---

### **Part 4: Modeling and Quantifying Heritable Phenotypes**

Once the time-series data has been rigorously preprocessed, the next critical step is to convert it into a set of stable, quantitative traits that are suitable for genetic analysis.

#### **Step 4.1: Establishing Phenotype Reliability with the Intraclass Correlation Coefficient (ICC)**
A fundamental prerequisite for a successful GWAS is that the phenotype of interest must be a stable, repeatable characteristic of an individual.
* **Action:** Use the cohort of ~2,000 participants with repeat measurements.
* **Method:**
    1.  Calculate all phenotypes for this cohort at both Time Point 1 and Time Point 2, running the entire pipeline independently for each time point.
    2.  For each phenotype, calculate the **Intra-class Correlation Coefficient (ICC)** between the Time 1 and Time 2 values.
* **Interpretation:** The ICC quantifies the reliability of the phenotype.
    * `ICC > 0.7`: High reliability. An excellent trait for GWAS.
    * `ICC 0.5-0.7`: Moderate reliability. A viable trait.
    * `ICC < 0.5`: Low reliability. The phenotype is too unstable over time and should be interpreted with extreme caution or discarded.
* **Decision:** Prioritize the phenotypes with the highest ICC for the main GWAS. It is likely that the "sleep-only" phenotypes and the Mesor will be the most reliable.

#### **Step 4.2: Phenotype Extraction via Harmonic Linear Modeling**
For each participant, fit a harmonic linear model to their 5-minute time-series data to extract the key parameters of their thermoregulatory rhythm. It is recommended to fit this model to both the full 24-hour data and the "sleep-only" data to generate two sets of phenotypes.

* **Model Specification (using R syntax):**
    ```R
    # 'time_radian' is the time of day converted to radians (0 to 2*pi)
    # 'ambient_temp' is the external weather data
    model <- lm(temp ~ sin(time_radian) + cos(time_radian) + ambient_temp)
    ```
* **Derived Phenotypes:** From the fitted `model` object for each person, extract:
    1.  **Mesor (Intercept):** The `(Intercept)` coefficient. Represents the rhythm- and environment-adjusted mean temperature.
    2.  **Amplitude (Rhythm Strength):** `sqrt(coef(model)["sin(time_radian)"]^2 + coef(model)["cos(time_radian)"]^2)`.
    3.  **Acrophase (Rhythm Timing):** `atan2(-coef(model)["sin(time_radian)"], coef(model)["cos(time_radian)"])`.
    4.  **Environmental Sensitivity:** The `coef(model)["ambient_temp"]`.
    5.  **Residual Variability (Stability):** The standard deviation of the model's residuals (`sd(residuals(model))`).

---

### **Part 5: Genome-Wide Association Study (GWAS) Design and Interpretation**

With robust and reliable phenotypes, the final analytical stage is to perform the GWAS to identify associated genetic variants.

#### **Step 5.1: Genotype Data QC and Imputation**
A standard, rigorous QC protocol will be applied to the raw genotype data using PLINK:
* **Variant-level QC:** Remove SNPs with low call rate (<98%), low MAF (<1%), or significant deviation from HWE ($p < 1 \times 10^{-6}$).
* **Individual-level QC:** Remove individuals with low call rate (<98%), sex discrepancies, or excessive heterozygosity.
* **Imputation:** Impute QC'd genotypes to a comprehensive reference panel (e.g., TopMed).

#### **Step 5.2: Optimal GWAS Software and Model Selection**
* **Software Recommendation:** **REGENIE** is recommended for its speed, scalability, and robust control of population structure and relatedness in biobank-scale datasets.
* **Covariate Strategy:** A tiered approach is necessary to distinguish direct effects from mediation.
    * **Primary Discovery Model:** `Phenotype ~ SNP + Age + Sex + PCs + Tech_Covariates + Season`.
    * **Secondary Conditional Models:** For significant SNPs from the primary model, add potential mediators (e.g., BMI, physical activity, sleep duration) as covariates to test for attenuation of the genetic signal.

#### **Step 5.3: Post-GWAS Analysis and Biological Validation**
1.  **Functional Annotation:** Use tools like FUMA to annotate significant loci and identify potential causal genes and regulatory effects (e.g., eQTLs).
2.  **Gene-Set and Pathway Enrichment:** Use tools like MAGMA to test if the GWAS results are enriched for genes in a priori biological pathways (circadian, inflammatory, metabolic).
3.  **Genetic Correlation:** Use LD Score Regression (LDSC) to estimate the genetic correlation ($r_g$) between the derived temperature phenotypes and other complex traits (e.g., insomnia, BMI, inflammatory diseases).

---

### **Part 6: Synthesis and Future Directions**

This plan outlines a comprehensive strategy, but inherent limitations exist, primarily related to unmeasured confounders (diet, room temperature, acute illness) and potential heterogeneity in wearable devices. These limitations can be partially addressed by including device type as a covariate and by treating the model's residual variance as a novel phenotype of temperature "lability."

Successful completion of this GWAS will open up several advanced research avenues:
* **Gene-Environment (GxE) Interactions:** Testing if SNP effects vary by season or activity level.
* **Multi-Trait GWAS:** Jointly analyzing the correlated temperature phenotypes to increase power.
* **Mendelian Randomization (MR):** Using the novel genetic instruments for temperature to test its causal effects on disease outcomes, moving from association to potential causality.


