#!/usr/bin/env bash
#SBATCH --job-name=bt2_idx_shard787980
#SBATCH --partition=icelake
#SBATCH --cpus-per-task=32
#SBATCH --mem=256G
#SBATCH --time=12:00:00
#SBATCH --array=0-2
#SBATCH --output=/rds/project/rds-FSWe0O6MFwc/users/hb676/logs/%x_%A_%a.out
#SBATCH --error=/rds/project/rds-FSWe0O6MFwc/users/hb676/logs/%x_%A_%a.err

# Task 0 -> shard 78 (first third of assemblies)
# Task 1 -> shard 79 (second third of assemblies)
# Task 2 -> shard 80 (final third of assemblies)

set -euo pipefail

module purge
module load bowtie/2.5.0

ASSEMBLY_ROOT="/home/hb676/rds/hpc-work/probe_design/data/assemblies/ncbi_dataset/data"
IDXDIR="/rds/project/rds-FSWe0O6MFwc/users/hb676/core_nt/bowtie2_by_vol"
TASK="${SLURM_ARRAY_TASK_ID:-0}"
SHARD=$(( 78 + TASK ))
PREFIX="${IDXDIR}/core_nt.${SHARD}"

mkdir -p "$IDXDIR"

# Skip if already built
if ls "${PREFIX}".*.bt2 1>/dev/null 2>&1 || ls "${PREFIX}".*.bt2l 1>/dev/null 2>&1; then
    echo "[$(date)] Index already exists at $PREFIX — skipping."
    exit 0
fi

# Collect all assembly FASTAs sorted
mapfile -t ALL_FNAS < <(find "$ASSEMBLY_ROOT" -name "*_genomic.fna" | sort)
N="${#ALL_FNAS[@]}"
echo "[$(date)] Found $N FASTA files total."

if [[ "$N" -eq 0 ]]; then
    echo "ERROR: No *_genomic.fna files found under $ASSEMBLY_ROOT" >&2
    exit 1
fi

# Split into three equal chunks
CHUNK=$(( (N + 2) / 3 ))   # ceiling division
START=$(( TASK * CHUNK ))
END=$(( START + CHUNK ))
[[ "$END" -gt "$N" ]] && END="$N"
FNAS=( "${ALL_FNAS[@]:$START:$(( END - START ))}" )
echo "[$(date)] Shard $(( 78 + TASK )): assemblies $START-$(( END - 1 )) of $N (${#FNAS[@]} files)"

# bowtie2-build accepts comma-separated list of FASTAs
FASTA_LIST=$(IFS=','; echo "${FNAS[*]}")

echo "[$(date)] Building index at $PREFIX ..."
echo "HOST=$(hostname)  THREADS=${SLURM_CPUS_PER_TASK:-32}"

bowtie2-build --large-index --threads "${SLURM_CPUS_PER_TASK:-32}" \
    "$FASTA_LIST" "$PREFIX"

echo "[$(date)] finished bowtie2-build"

if ls "${PREFIX}".*.bt2l 1>/dev/null 2>&1 || ls "${PREFIX}".*.bt2 1>/dev/null 2>&1; then
    echo "[$(date)] OK: index files created at $PREFIX"
else
    echo "ERROR: No index files found for prefix $PREFIX" >&2
    exit 3
fi

echo "[$(date)] done."
