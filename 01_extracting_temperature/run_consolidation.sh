
#!/bin/bash

# --- Configuration ---
# Replace with the name of your conda environment
CONDA_ENV_NAME="your_conda_env"

# Replace with the path to your data and where you want to save the results
INPUT_DIR="/mnt/project/temp_ukb_results/"
OUTPUT_DIR="/opt/notebooks/temp_ukb_consolidated/"

# --- Script Start ---

# Activate the conda environment
source activate "$CONDA_ENV_NAME"

# Check if activation was successful
if [ $? -ne 0 ]; then
    echo "Error: Could not activate conda environment '$CONDA_ENV_NAME'"
    echo "Please make sure the environment exists and is accessible."
    exit 1
fi

# Install the required packages
pip install pandas pyarrow

# Run the consolidation script
python consolidate_data.py --input_dir "$INPUT_DIR" --output_dir "$OUTPUT_DIR"

# Deactivate the environment (optional)
conda deactivate

echo "Consolidation complete. Files saved in $OUTPUT_DIR"
