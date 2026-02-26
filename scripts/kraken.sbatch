#!/bin/bash
#SBATCH --job-name=k2_1085736
#SBATCH --cpus-per-task=48
#SBATCH --mem=150G
#SBATCH --time=12:00:00
#SBATCH --output=/rds/project/rds-FSWe0O6MFwc/users/hb676/logs/k2_%A_%a.out
#SBATCH --error=/rds/project/rds-FSWe0O6MFwc/users/hb676/logs/k2_%A_%a.err
#SBATCH --array=0-16%2

set -euo pipefail

module purge
module load gcc/11
GCC_LIB=$(g++ -print-file-name=libstdc++.so.6)
export LD_LIBRARY_PATH="$(dirname "$GCC_LIB"):$LD_LIBRARY_PATH"
echo "Using libstdc++: $GCC_LIB"

# -------------------------
# Paths
# -------------------------
IN_DIR="/rds/project/rds-FSWe0O6MFwc/users/hb676/split_fasta/1085736"
OUT_DIR="/rds/project/rds-FSWe0O6MFwc/users/hb676/kraken2"
LOG_DIR="/rds/project/rds-FSWe0O6MFwc/users/hb676/logs"
DB="/home/hb676/rds/hpc-work/kraken2_db"

CONF="${CONF:-0.05}"

mkdir -p "$OUT_DIR" "$LOG_DIR"

# -------------------------
# Build file list
# -------------------------
mapfile -t FILES < <(
  find "$IN_DIR" -maxdepth 1 -type f -name "*.fasta" | sort
)

echo "Found ${#FILES[@]} input files"

if [ "${#FILES[@]}" -eq 0 ]; then
  echo "ERROR: No FASTA files found in $IN_DIR"
  exit 1
fi

TASK_ID="${SLURM_ARRAY_TASK_ID}"

if [ "$TASK_ID" -ge "${#FILES[@]}" ]; then
  echo "Task $TASK_ID exceeds file count ${#FILES[@]}; exiting."
  exit 0
fi

IN="${FILES[$TASK_ID]}"

if [ -z "$IN" ]; then
  echo "ERROR: IN is empty for task $TASK_ID"
  exit 1
fi

BASENAME="$(basename "$IN")"
STEM="${BASENAME%.*}"
CHUNK="1085736.${STEM}"

echo "Task ID: $TASK_ID"
echo "Input: $IN"
echo "Chunk: $CHUNK"
echo "START: $(date -Is)"

SECONDS=0

# -------------------------
# Run Kraken2 (timed)
# -------------------------
/usr/bin/time -v kraken2 \
  --db "$DB" \
  --threads "$SLURM_CPUS_PER_TASK" \
  --memory-mapping \
  --confidence "$CONF" \
  --report "$OUT_DIR/${CHUNK}.c${CONF}.report.txt" \
  --output "$OUT_DIR/${CHUNK}.c${CONF}.kraken.out" \
  "$IN"

rc=$?
echo "END: $(date -Is)"
echo "Elapsed_seconds=${SECONDS}"
echo "ExitCode=$rc"
