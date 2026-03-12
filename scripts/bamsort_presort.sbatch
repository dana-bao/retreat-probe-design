#!/usr/bin/env bash
#SBATCH --job-name=presort
#SBATCH --partition=icelake
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=36:00:00
#SBATCH --array=0-199%10
#SBATCH --output=/rds/project/rds-FSWe0O6MFwc/users/zh384/logs/%x_%A_%a.out
#SBATCH --error=/rds/project/rds-FSWe0O6MFwc/users/zh384/logs/%x_%A_%a.err

set -euo pipefail
export LC_ALL=C.UTF-8
export LANG=C.UTF-8

# ---- EDIT ----
SPECIES_ID="1398205"
BAMROOT="/rds/project/rds-FSWe0O6MFwc/users/hb676/bowtie2_bam"
PRESORT_ROOT="/rds/project/rds-FSWe0O6MFwc/users/hb676/bowtie2_presorted"

module purge
module use /usr/local/software/spack/spack-modules/icelake-20220710/linux-rocky8-icelake
module load samtools/1.14/gcc/v46lwk2d

THREADS="${SLURM_CPUS_PER_TASK:-8}"

SPECIES_PATH="$BAMROOT/$SPECIES_ID"
mapfile -t BAMS < <(find "$SPECIES_PATH" -name "*.bam" | sort)
N="${#BAMS[@]}"
echo "Found N=$N shard BAMs for species $SPECIES_ID"

TASK="${SLURM_ARRAY_TASK_ID:-0}"
if [[ "$TASK" -ge "$N" ]]; then
    echo "TASK=$TASK out of range (N=$N) -> exiting cleanly."
    exit 0
fi

BAM="${BAMS[$TASK]}"
echo "[$(date)] TASK=$TASK BAM=$BAM"

OUTDIR="$PRESORT_ROOT/$SPECIES_ID"
mkdir -p "$OUTDIR"

SHARD_OUT="$OUTDIR/shard_$(printf '%04d' "$TASK").sorted.sam.gz"

if [[ -s "$SHARD_OUT" ]] && gzip -t "$SHARD_OUT" 2>/dev/null; then
    echo "Already done: $SHARD_OUT — skipping."
    exit 0
fi

TMPDIR="$OUTDIR/tmp_presort_${SLURM_JOB_ID}_${TASK}"
mkdir -p "$TMPDIR"
trap 'rm -rf "$TMPDIR"' EXIT

echo "[$(date)] Sorting shard -> $SHARD_OUT"
samtools view -@ "$((THREADS - 2))" "$BAM" \
    | sort -k1,1 -T "$TMPDIR" -S 50G \
    | gzip -1 > "$SHARD_OUT"

echo "[$(date)] Done: $SHARD_OUT"
