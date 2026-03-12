#!/usr/bin/env bash
#SBATCH --job-name=pipeline_viz_comp
#SBATCH --partition=icelake
#SBATCH --cpus-per-task=1
#SBATCH --mem=128G
#SBATCH --time=4:00:00
#SBATCH --output=/rds/project/rds-FSWe0O6MFwc/users/hb676/logs/%x_%j.out
#SBATCH --error=/rds/project/rds-FSWe0O6MFwc/users/hb676/logs/%x_%j.err

set -euo pipefail

RDS="/rds/project/rds-FSWe0O6MFwc/users/hb676"
SCRIPTS="/home/hb676/rds/hpc-work/probe_design/scripts"
SPECIES="${SCRIPTS}/species_in_both_dbs.txt"

python3 "${SCRIPTS}/pipeline_visualize_compressed.py" \
  --species              "$SPECIES" \
  --kraken_dir           "${RDS}/kraken2" \
  --kraken_filt_dir      "${RDS}/kraken2_filtered" \
  --ngslca_dir           "${RDS}/ngslca_out_compressed" \
  --ngslca_mod_dir       "${RDS}/modern_ngslca_out_compressed" \
  --bamdam_dir           "${RDS}/bamdam_out_compressed" \
  --bamdam_mod_dir       "${RDS}/modern_bamdam_out_compressed" \
  --bamdam_filt_dir      "${RDS}/bamdam_filtered_compressed" \
  --bamdam_mod_filt_dir  "${RDS}/modern_bamdam_filtered_compressed" \
  --eprobe_input_dir     "${RDS}/eprobe_input" \
  --modern_overlap_dir   "${RDS}/modern_overlap_counts" \
  --out_dir              "${RDS}/figures_compressed_new"
