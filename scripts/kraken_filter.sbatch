#!/usr/bin/env bash
#SBATCH --job-name=kraken_filter
#SBATCH --partition=icelake
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH --output=/rds/project/rds-FSWe0O6MFwc/users/zh384/logs/%x_%j.out
#SBATCH --error=/rds/project/rds-FSWe0O6MFwc/users/zh384/logs/%x_%j.err

set -euo pipefail

python3 /home/hb676/rds/hpc-work/probe_design/scripts/kraken_filter.py \
  --results_dir /rds/project/rds-FSWe0O6MFwc/users/hb676/kraken2 \
  --nodes       /home/hb676/rds/hpc-work/kraken2_db/nodes.dmp \
  --out_dir     /rds/project/rds-FSWe0O6MFwc/users/hb676/kraken2_filtered
