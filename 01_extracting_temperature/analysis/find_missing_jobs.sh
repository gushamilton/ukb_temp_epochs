#!/usr/bin/env bash

# Compares the expected job chunks against a list of completed
# files and generates a list of jobs to resubmit.

CHUNK_SIZE=96
TOTAL_FILES=115355
COMPLETED_FILES="files"
RESUBMIT_LIST="resubmit_list.txt"

# --- Safety Check ---
if [ ! -f "$COMPLETED_FILES" ]; then
    echo "ERROR: The input file '$COMPLETED_FILES' was not found." >&2
    exit 1
fi

# Clear any previous resubmit list
> "$RESUBMIT_LIST"

echo "Checking for missing file chunks..."

for (( i=0; i<TOTAL_FILES; i+=CHUNK_SIZE )); do
    START=$((i + 1))
    END=$((i + CHUNK_SIZE))
    
    # Create the pattern to search for, e.g., "rows_1_96"
    PATTERN="rows_${START}_${END}"

    # Search for the pattern in the list of completed files.
    # The -q flag makes grep quiet; we only care about its exit code.
    grep -q "$PATTERN" "$COMPLETED_FILES"
    
    # If grep exits with a non-zero status ($?), the pattern was not found.
    if [ $? -ne 0 ]; then
        echo "Missing chunk found: ${START}-${END}. Adding to resubmit list."
        # Add the start and end numbers to our resubmit file
        echo "${START} ${END}" >> "$RESUBMIT_LIST"
    fi
done

echo ""
echo "Scan complete."
echo "A list of jobs to rerun has been created at: $RESUBMIT_LIST"
