#!/usr/bin/env bash
#SBATCH --job-name=compressbam_mod
#SBATCH --partition=icelake
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --array=0-6%7
#SBATCH --output=/rds/project/rds-FSWe0O6MFwc/users/zh384/logs/%x_%A_%a.out
#SBATCH --error=/rds/project/rds-FSWe0O6MFwc/users/zh384/logs/%x_%A_%a.err

# Compresses BAM headers (shards 00-80) for 7 modern species using compressbam.
# One array task per species; loops over all 81 shards internally.
# Output written as ${SAMPLE}_core_nt.${SHARD}.compressed.bam alongside the original.
# Skips shards where the compressed BAM already exists and is non-empty.
# Originals are never modified — replace manually after verifying output.

set -euo pipefail

BAMROOT="/rds/project/rds-FSWe0O6MFwc/users/hb676/modern_bowtie2_bam"
COMPRESSBAM="/rds/project/rds-FSWe0O6MFwc/users/hb676/scripts/compressbam"
THREADS="${SLURM_CPUS_PER_TASK:-4}"

SPECIES_LIST=(
    111071_modern
    1250597_modern
    344590_modern
    357475_modern
    39862_modern
    408822_modern
    42846_modern
)

N="${#SPECIES_LIST[@]}"
TASK="${SLURM_ARRAY_TASK_ID:-0}"

if [[ "$TASK" -ge "$N" ]]; then
    echo "TASK=$TASK out of range (N=$N) — exiting."
    exit 0
fi

SAMPLE="${SPECIES_LIST[$TASK]}"
SPECIES_PATH="${BAMROOT}/${SAMPLE}"
echo "[$(date)] TASK=$TASK  SAMPLE=$SAMPLE"
echo "HOST=$(hostname)  THREADS=$THREADS"

n_done=0
n_skip=0
n_fail=0

for SHARD_IDX in $(seq 0 80); do
    SHARD=$(printf "%02d" "$SHARD_IDX")
    input_bam="${SPECIES_PATH}/${SAMPLE}_core_nt.${SHARD}.bam"
    output_bam="${SPECIES_PATH}/${SAMPLE}_core_nt.${SHARD}.compressed.bam"

    if [[ ! -s "$input_bam" ]]; then
        continue   # BAM not present for this shard, skip silently
    fi

    if [[ -s "$output_bam" ]]; then
        (( n_skip++ )) || true
        continue   # already compressed in a previous run
    fi

    echo "  [$(date)] shard $SHARD: compressing $input_bam"

    if "$COMPRESSBAM" --threads "$THREADS" --input "$input_bam" --output "$output_bam"; then
        if [[ -s "$output_bam" ]]; then
            (( n_done++ )) || true
        else
            echo "  ERROR: empty output for shard $SHARD" >&2
            rm -f "$output_bam"
            (( n_fail++ )) || true
        fi
    else
        echo "  ERROR: compressbam failed for shard $SHARD (exit $?)" >&2
        rm -f "$output_bam"
        (( n_fail++ )) || true
    fi
done

echo "[$(date)] SAMPLE=$SAMPLE done: ${n_done} compressed, ${n_skip} skipped (already done), ${n_fail} failed"

if [[ "$n_fail" -gt 0 ]]; then
    exit 1
fi
