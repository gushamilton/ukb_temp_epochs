#!/usr/bin/env bash

# Submits ~1202 jobs to DNAnexus, each processing a 96-file chunk.
# This high-volume approach maximizes the chance of completion for low-priority jobs.

CHUNK_SIZE=96
TOTAL_FILES=115355
TEMPLATE_SCRIPT="run_production_1_100.sh" 
DESTINATION_FOLDER="/temp_ukb_results"

if [ ! -f "$TEMPLATE_SCRIPT" ]; then
    echo "ERROR: Template script '$TEMPLATE_SCRIPT' not found." >&2
    exit 1
fi

JOB_COUNT=0
for (( i=0; i<TOTAL_FILES; i+=CHUNK_SIZE )); do
    JOB_COUNT=$((JOB_COUNT + 1))
    START=$((i + 1))
    END=$((i + CHUNK_SIZE))
    
    TEMP_SCRIPT_NAME="run_production_${START}_${END}.sh"

    echo "Preparing and submitting job ${JOB_COUNT} for files ${START} to ${END}..."

    # 1. Create the temporary, parameterized script from the template
    sed "s/START_LINE=.*/START_LINE=${START}/; s/END_LINE=.*/END_LINE=${END}/" "$TEMPLATE_SCRIPT" > "$TEMP_SCRIPT_NAME"

    # 2. Upload the temporary script to the project
    dx upload "$TEMP_SCRIPT_NAME" --path "$DESTINATION_FOLDER/" --brief

    # 3. Run the job using the uploaded script with low priority
    dx run app-swiss-army-knife \
        --instance-type "mem1_ssd1_v2_x8" \
        --destination "$DESTINATION_FOLDER" \
        -iimage="continuumio/miniconda3:latest" \
        -icmd="bash /mnt/project/${DESTINATION_FOLDER}/${TEMP_SCRIPT_NAME}" \
        --priority "low" \
        -y --brief

    # 4. Clean up the local temporary script file
    rm "$TEMP_SCRIPT_NAME"

    # Pause for 1 second to avoid overwhelming the DNAnexus API
    sleep 1
done

echo "All ${JOB_COUNT} jobs have been submitted."
