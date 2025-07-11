
#THE SNAPSHOT NEEDED is r_env (on dna cloud workstation)
# Created snapshot: project-J1GYbG0JQz0gxQ6yZf1GqYkb:r_env (file-J1f2120JQz0Zxbb0j411Vx7F)

# option A – repeat the flag once per file (easiest)
dx run app-cloud_workstation --ssh \
  -isnapshot=r_env \
    -imax_session_length="8h" \
    -imax_session_length="24h" \
  -ifids=file-J1Y7QQjJQz0ZZ12YzPf1jx7x \
  -ifids=file-J1fJby0JVf195JxJ73g571Zk \
  --instance-type mem1_ssd1_v2_x16 \
  -ifids=file-J1Y94JjJQz0pjJxJkg0v5KPq  -y --allow-ssh 0.0.0.0/16

DX_PROJECT_CONTEXT_ID=project-J1GYbG0JQz0gxQ6yZf1GqYkb
unset DX_WORKSPACE_ID
dx cd $DX_PROJECT_CONTEXT_ID:

# Download the installer for Linux
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

# Run the installer script
bash Miniconda3-latest-Linux-x86_64.sh

# Follow the prompts. It's safe to accept the defaults.
# After installation, close and reopen your terminal for the changes to take effect.



# Create the environment with all packages at once
conda create -n r_env -c conda-forge \
r-base \
r-arrow \
r-tidyverse \
r-data.table \
r-lubridate \
r-zoo \
r-rcpproll \
r-lme4 \
r-brms \
r-patchwork \
r-broom \
r-r.utils \
r-corrplot \
r-ggpubr \
r-cowplot \
r-ICC \
quarto

# Activate your new environment
conda activate r_env



conda activate r_env

dx download file-J1Y7QQjJQz0ZZ12YzPf1jx7x
dx download file-J1Y94JjJQz0pjJxJkg0v5KPq
dx download file-J1fJby0JVf195JxJ73g571Zk