#!/usr/bin/env bash

# Reads only the FIRST line of 'resubmit_list.txt' and launches a
# single job to test the submission process.

RESUBMIT_LIST="resubmit_list.txt"
TEMPLATE_SCRIPT="run_production_1_100.sh"
DESTINATION_FOLDER="/temp_ukb_results"

# --- Safety Checks ---
if [ ! -f "$RESUBMIT_LIST" ]; then
    echo "ERROR: The resubmit list '$RESUBMIT_LIST' was not found." >&2
    echo "Please run 'find_missing_jobs.sh' first." >&2
    exit 1
fi
if [ ! -f "$TEMPLATE_SCRIPT" ]; then
    echo "ERROR: Template script '$TEMPLATE_SCRIPT' not found." >&2
    exit 1
fi

echo "Found $(wc -l < "$RESUBMIT_LIST") jobs to resubmit. Submitting the first one as a test..."

# Read the resubmit list line by line
while read -r START END; do

    TEMP_SCRIPT_NAME="run_production_${START}_${END}.sh"

    echo "Resubmitting job for files ${START} to ${END}..."

    sed "s/START_LINE=.*/START_LINE=${START}/; s/END_LINE=.*/END_LINE=${END}/" "$TEMPLATE_SCRIPT" > "$TEMP_SCRIPT_NAME"

    dx upload "$TEMP_SCRIPT_NAME" --path "$DESTINATION_FOLDER/" --brief

    # The --priority flag is removed to run this as a normal on-demand job
    dx run app-swiss-army-knife \
        --instance-type "mem1_ssd1_v2_x8" \
        --destination "$DESTINATION_FOLDER" \
        -iimage="continuumio/miniconda3:latest" \
        -icmd="bash /mnt/project/${DESTINATION_FOLDER}/${TEMP_SCRIPT_NAME}" \
        -y --brief

    rm "$TEMP_SCRIPT_NAME"

    sleep 1

    # Exit the loop after the first iteration
    break
done < "$RESUBMIT_LIST"

echo "Test job has been submitted."
