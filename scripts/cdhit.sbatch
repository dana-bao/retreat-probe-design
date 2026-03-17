#!/usr/bin/env bash
#SBATCH --job-name=cdhit_pooled
#SBATCH --partition=icelake
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/rds/project/rds-FSWe0O6MFwc/users/hb676/logs/%x_%j.out
#SBATCH --error=/rds/project/rds-FSWe0O6MFwc/users/hb676/logs/%x_%j.err

set -euo pipefail

module purge
module load cd-hit/4.6.8

IN_DIR="/rds/project/rds-FSWe0O6MFwc/users/hb676/deduplicate_input"
WORK_DIR="/rds/project/rds-FSWe0O6MFwc/users/hb676/cdhit_output"
mkdir -p "$WORK_DIR"

POOLED="${WORK_DIR}/pooled.fasta"
OUT="${WORK_DIR}/pooled_dedup"

# pool all species
echo "[$(date)] Pooling FASTA files..."
cat "${IN_DIR}"/*.fasta > "$POOLED"
echo "Total input: $(grep -c '^>' "$POOLED") sequences"

# run cd-hit-est on pooled file
echo "[$(date)] Running cd-hit-est..."
cd-hit-est \
    -i  "$POOLED" \
    -o  "$OUT" \
    -c  0.95 \
    -n  8 \
    -aS 1.0 \
    -r  1 \
    -d  0 \
    -T  "${SLURM_CPUS_PER_TASK}" \
    -M  0

echo "[$(date)] Done: $(grep -c '^>' "${OUT}.fasta") sequences after deduplication"
