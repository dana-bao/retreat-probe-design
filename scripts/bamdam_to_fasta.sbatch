#!/usr/bin/env bash
#SBATCH --job-name=bamdam_to_fasta
#SBATCH --partition=icelake
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --array=0-200%20
#SBATCH --output=/rds/project/rds-FSWe0O6MFwc/users/hb676/logs/%x_%A_%a.out
#SBATCH --error=/rds/project/rds-FSWe0O6MFwc/users/hb676/logs/%x_%A_%a.err

set -euo pipefail

# ---- tools ----
module purge
module use /usr/local/software/spack/spack-modules/icelake-20220710/linux-rocky8-icelake
module load samtools/1.14/gcc/v46lwk2d

# ---- paths ----
BAMDAM_OUTROOT="/rds/project/rds-FSWe0O6MFwc/users/hb676/bamdam_out"
FASTA_OUTROOT="/rds/project/rds-FSWe0O6MFwc/users/hb676/bamdam_fasta"

mkdir -p "$FASTA_OUTROOT"

# ---- build species list from completed shrunk BAMs ----
mapfile -t SPECIES_LIST < <(
    find "$BAMDAM_OUTROOT" -maxdepth 2 -name "*.shrunk.bam" -size +0c \
        | sed 's|.*/\([^/]*\)/.*|\1|' \
        | sort -u
)

N="${#SPECIES_LIST[@]}"
echo "Found N=$N species with completed shrunk BAMs"

TASK="${SLURM_ARRAY_TASK_ID:-0}"
if [[ "$TASK" -ge "$N" ]]; then
    echo "TASK=$TASK out of range (N=$N); exiting."
    exit 0
fi

SPECIES_ID="${SPECIES_LIST[$TASK]}"
echo "[$(date)] TASK=$TASK  SPECIES_ID=$SPECIES_ID"

SHRUNK_BAM="${BAMDAM_OUTROOT}/${SPECIES_ID}/${SPECIES_ID}.shrunk.bam"
OUT_FASTA="${FASTA_OUTROOT}/${SPECIES_ID}.shrunk.fasta"

if [[ ! -s "$SHRUNK_BAM" ]]; then
    echo "ERROR: shrunk BAM not found or empty: $SHRUNK_BAM"
    exit 1
fi

# ---- skip if already done ----
if [[ -s "$OUT_FASTA" ]]; then
    echo "Output already exists ($OUT_FASTA). Skipping."
    exit 0
fi

echo "[$(date)] Extracting FASTA from $SHRUNK_BAM"
samtools fasta \
    -F 2304 \
    -@ "${SLURM_CPUS_PER_TASK:-4}" \
    "$SHRUNK_BAM" \
    > "$OUT_FASTA"

echo "[$(date)] Done: $OUT_FASTA"
echo "Reads extracted: $(grep -c '^>' "$OUT_FASTA")"
