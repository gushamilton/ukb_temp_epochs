#!/usr/bin/env bash
set -euo pipefail

##############################################################################
#
# run_production_final.sh
#
# Processes UK Biobank Raw CWA files in parallel and outputs to the
# highly efficient Parquet format. This script is optimized for the full
# analysis run.
#
# FINAL OUTPUTS (UPLOADED TO DNANEXUS):
# 1. temperature_epochs.parquet (zstd compressed)
# 2. summary_covariates.csv.gz
#
##############################################################################

############### CONFIG #######################################################
# The hardcoded full path to the text file containing the list of all CWA files.
LIST_FILE_PATH="/mnt/project/raw_activity_files.txt"

# The DNAnexus folder where final results will be uploaded.
OUT_DX_FOLDER="/temp_ukb_results"

# Number of parallel jobs. Should match the vCPU count of your RAP instance.
# For a full run, consider a larger instance like mem1_ssd1_v2_x16 and set this to 16.
PARALLEL_JOBS=8

# A name for this run, used in output filenames. e.g., "full_cohort"
# For testing, we'll keep this as "first100".
RUN_NAME="first100"

# Path for the conda environment.
ENV_DIR="$HOME/accenv_parquet"
##############################################################################


echo "--- Script starting at $(date) ---"

# --- Validate that the list file exists before we begin ---
if [ ! -f "$LIST_FILE_PATH" ]; then
    echo "ERROR: The list of raw files was not found at the hardcoded location:" >&2
    echo "  ${LIST_FILE_PATH}" >&2
    echo "Please ensure the file exists at this path before running." >&2
    exit 1
fi


########################################
# 0. Build or activate conda env       #
########################################
if ! command -v conda &>/dev/null ; then
    echo "ERROR: conda not found." >&2; exit 1
fi
source "$(conda info --base)/etc/profile.d/conda.sh"

# Force a clean environment by removing any old one that might exist.
if [ -d "$ENV_DIR" ]; then
    echo "--- An old environment exists at $ENV_DIR. Removing it for a clean start..."
    rm -rf "$ENV_DIR"
fi

echo "--- Creating new conda environment at $ENV_DIR..."
conda create -y -p "$ENV_DIR" -c conda-forge python=3.9 openjdk pip pyarrow parallel jq zstd
conda activate "$ENV_DIR"
echo "--- Installing python packages..."
pip install --no-cache-dir accelerometer

export OPENBLAS_NUM_THREADS=1

echo "--- Compiling Java parser..."
python -c '
import importlib.util, pathlib, subprocess, glob, sys
try:
    import accelerometer
except ImportError:
    print("ERROR: Could not import accelerometer package.", file=sys.stderr); sys.exit(1)
jdir = pathlib.Path(accelerometer.__file__).parent / "java"
jar = next(jdir.glob("JTransforms*.jar"))
java_files = glob.glob(str(jdir / "*.java"))
subprocess.run(["javac", "-cp", str(jar), *java_files], check=True, capture_output=True)
'
echo "--- Environment setup complete."


########################################
# 1. Define per-file worker function   #
########################################
# A dedicated log file will be used to track temporary directories.
WORK_DIR_LOG=$(mktemp -p "${TMPDIR:-/tmp}" "work_dirs_XXXX.log")

process_one_file () {
    local RAW_FILE_PATH=$1
    local ID
    ID=$(basename "$RAW_FILE_PATH" .cwa)
    local WORK_DIR
    WORK_DIR=$(mktemp -d -p "${TMPDIR:-/tmp}" "${ID}_XXXX")

    # The entire worker function's output is redirected to a per-job log file
    # to prevent any stray output from being captured by `parallel`.
    {
        pushd "$WORK_DIR" >/dev/null

        accProcess "$RAW_FILE_PATH" \
        --activityClassification False \
        --deleteIntermediateFiles False \
        --outputFolder .

        # Python script to create intermediate parts (Parquet only)
        python -c '
import pandas as pd
import json
import sys

file_id_with_path = sys.argv[1]
file_id = file_id_with_path.split("/")[-1].replace(".cwa", "")

epoch_csv_path = f"{file_id}-epoch.csv.gz"
summary_json_path = f"{file_id}-summary.json"

# --- Process Epoch Data ---
epoch_cols = ["time", "temp", "enmoTrunc", "dataErrors", "samples", "clipsBeforeCalibr", "clipsAfterCalibr"]
df_epoch = pd.read_csv(epoch_csv_path, usecols=epoch_cols)
df_epoch.insert(0, "id", file_id)

df_epoch.rename(columns={
    "clipsBeforeCalibr": "clipsBefore",
    "clipsAfterCalibr": "clipsAfter"
}, inplace=True)

# OPTIMIZATION: Write only to Parquet format
df_epoch.to_parquet("epoch.part.parquet", engine="pyarrow", compression="zstd")

# --- Process Summary Data ---
with open(summary_json_path, "r") as f:
    summary_data = json.load(f)
summary_record = {
    "id": file_id, "goodWear": summary_data.get("quality-goodWearTime"),
    "goodCal": summary_data.get("quality-goodCalibration"),
    "calErrAfter_mg": summary_data.get("calibration-errsAfter(mg)"),
    "wearDays": summary_data.get("wearTime-overall(days)"),
    "xOffset_g": summary_data.get("calibration-xOffset(g)"),
    "yOffset_g": summary_data.get("calibration-yOffset(g)"),
    "zOffset_g": summary_data.get("calibration-zOffset(g)"),
    "clipsBefore": summary_data.get("clipsBeforeCalibration"),
    "clipsAfter": summary_data.get("clipsAfterCalibration"),
    "startTime": summary_data.get("file-startTime"),
    "endTime": summary_data.get("file-endTime"),
}
pd.DataFrame([summary_record]).to_csv("summary.part.csv", index=False)
' "$RAW_FILE_PATH"

        popd >/dev/null

    } > "${WORK_DIR}.log" 2>&1

    # Atomically append the path of the successfully created temp directory to the log file.
    echo "$WORK_DIR" >> "$WORK_DIR_LOG"
}
# Export the function and the variable separately.
export -f process_one_file
export WORK_DIR_LOG


########################################
# 2. Run processing in parallel        #
########################################
# Set to process the first 100 files for the pilot run.
head -n 10 "$LIST_FILE_PATH" > "/tmp/pilot_list.txt"
echo "--- Starting parallel processing for ${RUN_NAME} run..."
# No longer capturing stdout. Each job writes its temp dir path to a log file.
parallel -j "$PARALLEL_JOBS" --halt soon,fail=1 --eta process_one_file :::: /tmp/pilot_list.txt
echo "--- Parallel processing complete."


########################################
# 3. Aggregate results                 #
########################################
echo "--- Aggregating results from all completed jobs..."
# Create a single, robust output directory for the final aggregated files.
FINAL_OUTPUT_DIR=$(mktemp -d -p "${TMPDIR:-/tmp}" "final_results_XXXX")
FINAL_EPOCH_PARQUET="${FINAL_OUTPUT_DIR}/temperature_epochs_${RUN_NAME}.parquet"
FINAL_SUMMARY_CSV="${FINAL_OUTPUT_DIR}/summary_covariates_${RUN_NAME}.csv"

# Aggregate using Python for robustness
python -c '
import pandas as pd
import glob
import sys
import os

work_dirs_log_path = sys.argv[1]
final_parquet_path = sys.argv[2]
final_summary_path = sys.argv[3]

# Read the list of temporary directories from the log file
with open(work_dirs_log_path, "r") as f:
    work_dirs = [line.strip() for line in f if line.strip()]

# --- Aggregate Epoch Parquet files ---
parquet_parts = [p for d in work_dirs for p in glob.glob(os.path.join(d, "epoch.part.parquet"))]
if parquet_parts:
    pd.concat(
        [pd.read_parquet(p) for p in parquet_parts], ignore_index=True
    ).to_parquet(final_parquet_path, engine="pyarrow", compression="zstd")
    print(f"Aggregated Parquet data saved to {final_parquet_path}")

# --- Aggregate Summary CSV files ---
summary_parts = [p for d in work_dirs for p in glob.glob(os.path.join(d, "summary.part.csv"))]
if summary_parts:
    pd.concat(
        [pd.read_csv(p) for p in summary_parts], ignore_index=True
    ).to_csv(final_summary_path, index=False)
    print(f"Aggregated summary data saved to {final_summary_path}")

' "$WORK_DIR_LOG" "$FINAL_EPOCH_PARQUET" "$FINAL_SUMMARY_CSV"

# Compress the final summary file
echo "--- Compressing aggregated summary file..."
gzip -f "$FINAL_SUMMARY_CSV"


########################################
# 4. Upload final results              #
########################################
echo "--- Uploading final result files to DNAnexus folder: ${OUT_DX_FOLDER}"
dx mkdir -p "$OUT_DX_FOLDER"

dx upload "${FINAL_EPOCH_PARQUET}" --path "${OUT_DX_FOLDER}/" --brief --wait
dx upload "${FINAL_SUMMARY_CSV}.gz" --path "${OUT_DX_FOLDER}/" --brief --wait
echo "--- Upload complete."


########################################
# 5. Cleanup                           #
########################################
echo "--- Cleaning up temporary directories..."
# Read the log file to get the list of directories to clean up
if [ -f "$WORK_DIR_LOG" ]; then
    while read -r wd; do
        if [ -d "$wd" ]; then
            rm -rf "$wd"
            # Also remove the per-job log file
            rm -f "${wd}.log"
        fi
    done < "$WORK_DIR_LOG"
    rm -f "$WORK_DIR_LOG"
fi
rm -rf "$FINAL_OUTPUT_DIR"
rm -f /tmp/pilot_list.txt

echo "--- Script finished successfully at $(date) ---"
