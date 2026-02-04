#!/bin/bash
#SBATCH --job-name=aedna_tile
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=12:00:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

set -euo pipefail

# 1) Load python and activate venv
module load python/3.9.12/gcc/pdcqf4o5
source ~/venvs/biopython_py39/bin/activate

# 2) Inputs
LIST="assemblies.list"
ASSEMBLY="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST")"

# Adjust these paths to match your layout
METADATA="../data/all_availability2.tsv"
OUTDIR="../data/tiled_reads"

mkdir -p "$OUTDIR"

# 3) Run the script on one assembly
python3 aeDNA_simulation.py \
  --assembly "$ASSEMBLY" \
  --metadata "$METADATA" \
  --output_dir "$OUTDIR" \
  --deamination \
  --mutation
