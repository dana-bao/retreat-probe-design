#!/bin/bash
#SBATCH --job-name=aedna_tile
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

set -euo pipefail

module purge || true
source ~/venvs/biopython_py39/bin/activate
PY=~/venvs/biopython_py39/bin/python

LIST="${1:?Usage: sbatch --array=1-N slurm_tile.sbatch batches/batch_XXX.list}"

ASSEMBLY="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST" | tr -d '\r' || true)"
if [[ -z "${ASSEMBLY}" ]]; then
  echo "ERROR: No assembly found for task ${SLURM_ARRAY_TASK_ID} in ${LIST}"
  exit 1
fi

METADATA="../data/all_availability2.tsv"
OUTDIR="../data/tiled_reads"
mkdir -p "$OUTDIR" logs

echo "HOST=$(hostname)"
echo "PWD=$(pwd)"
echo "PY=$PY"
"$PY" -c "import sys; print('python:', sys.executable); print(sys.version)"
echo "ASSEMBLY=${ASSEMBLY}"

"$PY" aeDNA_simulation.py \
  --assembly "$ASSEMBLY" \
  --metadata "$METADATA" \
  --output_dir "$OUTDIR" \
  --deamination \
  --mutation
