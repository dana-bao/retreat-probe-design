#!/usr/bin/env bash
#SBATCH --job-name=bamdam_filter
#SBATCH --partition=icelake
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH --output=/rds/project/rds-FSWe0O6MFwc/users/hb676/logs/%x_%j.out
#SBATCH --error=/rds/project/rds-FSWe0O6MFwc/users/hb676/logs/%x_%j.err

set -euo pipefail

# ---- EDIT: set MODE to "aedna" or "modern" ----
MODE="aedna"

if [[ "$MODE" == "modern" ]]; then
    LCA_DIR="/rds/project/rds-FSWe0O6MFwc/users/hb676/modern_bamdam_out_new"
    OUT_DIR="/rds/project/rds-FSWe0O6MFwc/users/hb676/modern_bamdam_filtered_new"
else
    LCA_DIR="/rds/project/rds-FSWe0O6MFwc/users/hb676/bamdam_out_new"
    OUT_DIR="/rds/project/rds-FSWe0O6MFwc/users/hb676/bamdam_filtered_new"
fi

NODES="/rds/project/rds-FSWe0O6MFwc/users/hb676/kraken2_taxonomy/nodes.dmp"

echo "MODE=$MODE  LCA_DIR=$LCA_DIR  OUT_DIR=$OUT_DIR"

python3 /home/hb676/rds/hpc-work/probe_design/scripts/bamdam_filter.py \
  --lca_dir "$LCA_DIR" \
  --nodes   "$NODES" \
  --out_dir "$OUT_DIR"
