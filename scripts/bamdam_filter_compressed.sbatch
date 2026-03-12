#!/usr/bin/env bash
#SBATCH --job-name=bamdam_filt_comp
#SBATCH --partition=icelake
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=8:00:00
#SBATCH --output=/rds/project/rds-FSWe0O6MFwc/users/hb676/logs/%x_%j.out
#SBATCH --error=/rds/project/rds-FSWe0O6MFwc/users/hb676/logs/%x_%j.err

set -euo pipefail

RDS="/rds/project/rds-FSWe0O6MFwc/users/hb676"
BAMDAM_DIR="${RDS}/bamdam_out_compressed"
OUT_DIR="${RDS}/bamdam_filtered_compressed"
NODES="${RDS}/kraken2_taxonomy/nodes.dmp"
SCRIPTS="/home/hb676/rds/hpc-work/probe_design/scripts"

mkdir -p "$OUT_DIR"

# Reads chunk .shrunk.lca files directly — no temp files written to disk
python3 "${SCRIPTS}/bamdam_filter_compressed.py" \
    --lca_dir "$BAMDAM_DIR" \
    --nodes   "$NODES" \
    --out_dir "$OUT_DIR"
