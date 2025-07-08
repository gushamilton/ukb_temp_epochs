
#THE SNAPSHOT NEEDED is r_env (on dna cloud workstation)
# Created snapshot: project-J1GYbG0JQz0gxQ6yZf1GqYkb:r_env (file-J1f2120JQz0Zxbb0j411Vx7F)

# option A – repeat the flag once per file (easiest)
dx run app-cloud_workstation --ssh \
  -isnapshot=r_env \
  -ifids=combined_temperature_data.parquet \
  -ifids=combined_covariates.csv.gz
# this runs once to make a job

# # Download the installer for Linux
# wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

# # Run the installer script
# bash Miniconda3-latest-Linux-x86_64.sh

# # Follow the prompts. It's safe to accept the defaults.
# # After installation, close and reopen your terminal for the changes to take effect.

# # Create the environment with all packages at once
# conda create -n r_env -c conda-forge r-base r-arrow r-dplyr r-data.table r-tidyr r-purrr r-lubridate r-zoo r-rcpproll

# Activate your new environment
conda activate r_env
