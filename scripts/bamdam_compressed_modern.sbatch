#!/usr/bin/env bash
#SBATCH --job-name=bamdam_comp_mod
#SBATCH --partition=icelake
#SBATCH --cpus-per-task=36
#SBATCH --mem=16G
#SBATCH --time=36:00:00
#SBATCH --array=0-27
#SBATCH --output=/rds/project/rds-FSWe0O6MFwc/users/zh384/logs/%x_%A_%a.out
#SBATCH --error=/rds/project/rds-FSWe0O6MFwc/users/zh384/logs/%x_%A_%a.err

set -euo pipefail

# ---- tools ----
module purge
module use /usr/local/software/spack/spack-modules/icelake-20220710/linux-rocky8-icelake
module load samtools/1.14/gcc/v46lwk2d

# ---- PATHS ----
SORTEDROOT="/rds/project/rds-FSWe0O6MFwc/users/hb676/modern_bowtie2_sorted_compressed"
LCA_OUTROOT="/rds/project/rds-FSWe0O6MFwc/users/hb676/modern_ngslca_out_compressed"
BAMDAM_OUTROOT="/rds/project/rds-FSWe0O6MFwc/users/hb676/modern_bamdam_out_compressed"
BAMDAM="/rds/project/rds-FSWe0O6MFwc/users/zw442/.envs/dnadamage/bin/bamdam"

STRANDEDNESS="ss"
UPTO="genus"

# ---- Hardcoded: 7 modern species × 4 chunks = 28 tasks ----
SPECIES_LIST=(
    111071_modern
    1250597_modern
    344590_modern
    357475_modern
    39862_modern
    408822_modern
    42846_modern
)

# ---- build item list: {species_id}_chunk{00-03} for each species ----
ITEM_LIST=()
for sid in "${SPECIES_LIST[@]}"; do
    for i in 00 01 02 03; do
        ITEM_LIST+=("${sid}_chunk${i}")
    done
done

N="${#ITEM_LIST[@]}"
echo "Total items: $N"

TASK="${SLURM_ARRAY_TASK_ID:-0}"
if [[ "$TASK" -ge "$N" ]]; then
    echo "TASK=$TASK out of range (N=$N); exiting."
    exit 0
fi

ITEM="${ITEM_LIST[$TASK]}"
echo "[$(date)] TASK=$TASK  ITEM=$ITEM"

SPECIES_ID=$(echo "$ITEM" | sed 's/_chunk[0-9]*//')
CHUNK_TAG=$(echo "$ITEM" | grep -oP 'chunk\d+')

OUTDIR="${BAMDAM_OUTROOT}/${SPECIES_ID}"
mkdir -p "$OUTDIR"

CHUNK_SHRUNK_BAM="${OUTDIR}/${ITEM}.shrunk.bam"
CHUNK_SHRUNK_LCA="${OUTDIR}/${ITEM}.shrunk.lca"

if [[ -s "$CHUNK_SHRUNK_BAM" ]]; then
    echo "Already done ($CHUNK_SHRUNK_BAM). Skipping."
    exit 0
fi

CHUNK_LCA="${LCA_OUTROOT}/${ITEM}.lca"
if [[ ! -f "$CHUNK_LCA" ]]; then
    echo "ERROR: LCA file not found: $CHUNK_LCA"
    exit 1
fi

# Check LCA has data (not just a header)
DATA_LINES=$(grep -cv '^#' "$CHUNK_LCA" 2>/dev/null || true)
if [[ "$DATA_LINES" -eq 0 ]]; then
    echo "WARNING: $CHUNK_LCA has no data lines — creating empty outputs."
    BAM_FOR_HDR="${SORTEDROOT}/${SPECIES_ID}/${SPECIES_ID}.${CHUNK_TAG}.sorted.bam"
    samtools view -b -H "$BAM_FOR_HDR" > "$CHUNK_SHRUNK_BAM"
    touch "$CHUNK_SHRUNK_LCA"
    echo "[$(date)] Done (empty): $CHUNK_SHRUNK_BAM"
    exit 0
fi

CHUNK_SRC="${SORTEDROOT}/${SPECIES_ID}/${SPECIES_ID}.${CHUNK_TAG}.sorted.bam"
if [[ ! -s "$CHUNK_SRC" ]]; then
    echo "ERROR: no BAM found: $CHUNK_SRC"
    exit 1
fi

echo "[$(date)] bamdam shrink: $ITEM"
"$BAMDAM" shrink \
    --in_bam  "$CHUNK_SRC" \
    --in_lca  "$CHUNK_LCA" \
    --out_bam "$CHUNK_SHRUNK_BAM" \
    --out_lca "$CHUNK_SHRUNK_LCA" \
    --stranded "$STRANDEDNESS" \
    --upto "$UPTO"

echo "[$(date)] Done: $CHUNK_SHRUNK_BAM"
