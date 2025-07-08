#!/usr/bin/env bash
set -euo pipefail

# ────────────────────────────────────────────────────────────────
#  This edition can be executed either
#    • as a DNAnexus applet   (env-vars $start_line / $end_line) or
#    • in a Swiss-Army-Knife  (two positional arguments)
# ----------------------------------------------------------------
#  Example SAK call for rows 1-100:
#      /mnt/project/run_production_final.sh 1 100
# ────────────────────────────────────────────────────────────────

# Accept CLI args when present (SAK) ─────────────────────────────
if [[ $# -ge 2 ]]; then
  export start_line="$1"
  export end_line="$2"
  shift 2
fi
: "${start_line:?need start_line}"
: "${end_line:?need end_line}"

# Allow the master list path to be overridden via env-var; otherwise default
: "${LIST_FILE_PATH:=/mnt/project/raw_activity_files.txt}"

##############################################################################
#
#  run_production_final.sh
#
#  Inputs
#     $start_line  (1-based, inclusive)
#     $end_line    (inclusive)
#
##############################################################################

############### SLOT-SPECIFIC CONFIG #########################################
START_LINE="$start_line"
END_LINE="$end_line"
RUN_NAME="rows_${START_LINE}_${END_LINE}"

# Auto-detect core count from the DNAnexus VM
PARALLEL_JOBS=${DX_WORKSPACE_VCPU_COUNT:-8}
##############################################################################

############### GLOBAL CONFIG ################################################
OUT_DX_FOLDER="/temp_ukb_results"      # where final files land in the project
ENV_DIR="$HOME/accenv_parquet"         # transient conda env
##############################################################################

echo "--- Script starting at $(date) for lines ${START_LINE}-${END_LINE} ---"
echo "LIST_FILE_PATH=$LIST_FILE_PATH  (override with  export LIST_FILE_PATH=… )"

if [ ! -f "$LIST_FILE_PATH" ]; then
    echo "ERROR: master list not found at ${LIST_FILE_PATH}" >&2
    exit 1
fi

# ---------- Conda environment ----------------------------------------------
if ! command -v conda &>/dev/null ; then
    echo "ERROR: conda not found." >&2
    exit 1
fi
source "$(conda info --base)/etc/profile.d/conda.sh"

rm -rf "$ENV_DIR" 2>/dev/null || true
echo "--- Creating conda env at $ENV_DIR ..."
conda create -y -p "$ENV_DIR" -c conda-forge python=3.9 openjdk pip \
              pyarrow parallel jq zstd
conda activate "$ENV_DIR"
pip install --no-cache-dir accelerometer
export OPENBLAS_NUM_THREADS=1

echo "--- Compiling accelerometer Java helper ..."
python - <<'PY'
import importlib.util, pathlib, subprocess, glob, sys
import accelerometer
jdir = pathlib.Path(accelerometer.__file__).parent / "java"
jar  = next(jdir.glob("JTransforms*.jar"))
subprocess.run(["javac", "-cp", str(jar), *glob.glob(str(jdir / "*.java"))],
               check=True)
PY

# ---------- Worker function -------------------------------------------------
WORK_DIR_LOG=$(mktemp -p "${TMPDIR:-/tmp}" "work_dirs_XXXX.log")

process_one_file () {
    local RAW_FILE_PATH=$1
    local ID=$(basename "$RAW_FILE_PATH" .cwa)
    local WORK_DIR
    WORK_DIR=$(mktemp -d -p "${TMPDIR:-/tmp}" "${ID}_XXXX")

    {
        cd "$WORK_DIR"

        accProcess "$RAW_FILE_PATH" \
          --activityClassification False \
          --deleteIntermediateFiles False \
          --outputFolder .

        python - <<'PY' "$RAW_FILE_PATH"
import pandas as pd, json, sys, os
raw = sys.argv[1]
fid = os.path.basename(raw).replace(".cwa", "")
epoch_csv = f"{fid}-epoch.csv.gz"
summary_js = f"{fid}-summary.json"

cols = ["time","temp","enmoTrunc","dataErrors","samples",
        "clipsBeforeCalibr","clipsAfterCalibr"]
df = pd.read_csv(epoch_csv, usecols=cols)
df.insert(0,"id",fid)
df.rename(columns={"clipsBeforeCalibr":"clipsBefore",
                   "clipsAfterCalibr":"clipsAfter"}, inplace=True)
df.to_parquet("epoch.part.parquet", engine="pyarrow", compression="zstd")

with open(summary_js) as f:
    s = json.load(f)
row = dict(id=fid,
           goodWear=s.get("quality-goodWearTime"),
           goodCal=s.get("quality-goodCalibration"),
           calErrAfter_mg=s.get("calibration-errsAfter(mg)"),
           wearDays=s.get("wearTime-overall(days)"),
           xOffset_g=s.get("calibration-xOffset(g)"),
           yOffset_g=s.get("calibration-yOffset(g)"),
           zOffset_g=s.get("calibration-zOffset(g)"),
           clipsBefore=s.get("clipsBeforeCalibration"),
           clipsAfter=s.get("clipsAfterCalibration"),
           startTime=s.get("file-startTime"),
           endTime=s.get("file-endTime"))
pd.DataFrame([row]).to_csv("summary.part.csv", index=False)
PY
    } > "${WORK_DIR}.log" 2>&1

    echo "$WORK_DIR" >> "$WORK_DIR_LOG"
}
export -f process_one_file
export WORK_DIR_LOG

# ---------- Sub-set the master list ----------------------------------------
SUBSET_LIST_FILE="/tmp/subset_list.txt"
sed -n "${START_LINE},${END_LINE}p" "$LIST_FILE_PATH" > "$SUBSET_LIST_FILE"

# ---------- Fan-out execution ----------------------------------------------
echo "--- Processing $(wc -l < "$SUBSET_LIST_FILE") files with $PARALLEL_JOBS threads ..."
parallel -j "$PARALLEL_JOBS" --halt soon,fail=1 --eta \
         process_one_file :::: "$SUBSET_LIST_FILE"
echo "--- Per-file jobs finished ---"

# ---------- Aggregation -----------------------------------------------------
FINAL_OUTPUT_DIR=$(mktemp -d -p "${TMPDIR:-/tmp}" "final_${RUN_NAME}_XXXX")
FINAL_EPOCH_PARQUET="${FINAL_OUTPUT_DIR}/temperature_epochs_${RUN_NAME}.parquet"
FINAL_SUMMARY_CSV="${FINAL_OUTPUT_DIR}/summary_covariates_${RUN_NAME}.csv"

python - <<'PY' "$WORK_DIR_LOG" "$FINAL_EPOCH_PARQUET" "$FINAL_SUMMARY_CSV"
import pandas as pd, glob, sys, os
log, out_pq, out_csv = sys.argv[1:]
with open(log) as f:
    work_dirs=[l.strip() for l in f if l.strip()]

parts=[p for d in work_dirs for p in glob.glob(os.path.join(d,"epoch.part.parquet"))]
if parts:
    pd.concat([pd.read_parquet(p) for p in parts], ignore_index=True).to_parquet(
        out_pq, engine="pyarrow", compression="zstd")
parts=[p for d in work_dirs for p in glob.glob(os.path.join(d,"summary.part.csv"))]
if parts:
    pd.concat([pd.read_csv(p) for p in parts], ignore_index=True).to_csv(
        out_csv, index=False)
PY

gzip -f "$FINAL_SUMMARY_CSV"

# ---------- Upload to DNAnexus --------------------------------------------
echo "--- Uploading to $OUT_DX_FOLDER ..."
dx mkdir -p "$OUT_DX_FOLDER"
dx upload "$FINAL_EPOCH_PARQUET"     --path "$OUT_DX_FOLDER/" --brief --wait
dx upload "${FINAL_SUMMARY_CSV}.gz"  --path "$OUT_DX_FOLDER/" --brief --wait

# ---------- Cleanup ---------------------------------------------------------
echo "--- Cleaning up ..."
if [ -f "$WORK_DIR_LOG" ]; then
  while read -r wd; do
    rm -rf "$wd" "${wd}.log"
  done < "$WORK_DIR_LOG"
  rm -f "$WORK_DIR_LOG"
fi
rm -rf "$FINAL_OUTPUT_DIR" "$SUBSET_LIST_FILE"

echo "--- Slice ${START_LINE}-${END_LINE} finished at $(date) ---"
