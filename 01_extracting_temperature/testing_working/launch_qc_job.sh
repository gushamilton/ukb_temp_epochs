#!/usr/bin/env bash

# Submits a single, on-demand job to process files 1-300 for QC purposes.

START_LINE=1
END_LINE=300
TEMPLATE_SCRIPT="run_production_1_100.sh"
DESTINATION_FOLDER="/temp_ukb_results"

# --- Safety Check ---
if [ ! -f "$TEMPLATE_SCRIPT" ]; then
    echo "ERROR: Template script '$TEMPLATE_SCRIPT' not found." >&2
    exit 1
fi

TEMP_SCRIPT_NAME="run_production_${START_LINE}_${END_LINE}.sh"

echo "Preparing and submitting ON-DEMAND job for files ${START_LINE} to ${END_LINE}..."

# 1. Create the temporary, parameterized script from the template
sed "s/START_LINE=.*/START_LINE=${START_LINE}/; s/END_LINE=.*/END_LINE=${END_LINE}/" "$TEMPLATE_SCRIPT" > "$TEMP_SCRIPT_NAME"

# 2. Upload the temporary script to the project
dx upload "$TEMP_SCRIPT_NAME" --path "$DESTINATION_FOLDER/" --brief

# 3. Run the job using the uploaded script
dx run app-swiss-army-knife \
    --instance-type "mem1_ssd1_v2_x8" \
    --destination "$DESTINATION_FOLDER" \
    -iimage="continuumio/miniconda3:latest" \
    -icmd="bash /mnt/project/${DESTINATION_FOLDER}/${TEMP_SCRIPT_NAME}" \
    -y --brief

# 4. Clean up the local temporary script file
rm "$TEMP_SCRIPT_NAME"

echo "Job for files ${START_LINE}-${END_LINE} has been submitted."
