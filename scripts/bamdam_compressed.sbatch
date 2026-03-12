#!/usr/bin/env bash
#SBATCH --job-name=bamdam_comp
#SBATCH --partition=icelake
#SBATCH --cpus-per-task=36
#SBATCH --mem=16G
#SBATCH --time=36:00:00
#SBATCH --array=0-303%20
#SBATCH --output=/rds/project/rds-FSWe0O6MFwc/users/zh384/logs/%x_%A_%a.out
#SBATCH --error=/rds/project/rds-FSWe0O6MFwc/users/zh384/logs/%x_%A_%a.err

set -euo pipefail

# ---- tools ----
module purge
module use /usr/local/software/spack/spack-modules/icelake-20220710/linux-rocky8-icelake
module load samtools/1.14/gcc/v46lwk2d

# ---- PATHS ----
SORTEDROOT="/rds/project/rds-FSWe0O6MFwc/users/hb676/bowtie2_sorted_compressed"
LCA_OUTROOT="/rds/project/rds-FSWe0O6MFwc/users/hb676/ngslca_out_compressed"
BAMDAM_OUTROOT="/rds/project/rds-FSWe0O6MFwc/users/hb676/bamdam_out_compressed"
BAMDAM="/rds/project/rds-FSWe0O6MFwc/users/zw442/.envs/dnadamage/bin/bamdam"

STRANDEDNESS="ss"
UPTO="genus"

# ---- Hardcoded species list (76 species with completed 4-chunk ngsLCA output) ----
SPECIES_LIST=(
    1085736 111071 1118155 1250597 127563 129212 1293966 134639 135145
    1351830 139363 1398205 140813 1483842 1532265 1534596 1553662 157178
    157520 164537 1715899 175055 1753982 184104 185980 186720 1875977
    1875987 187635 1886436 1963269 199609 226588 229869 242528 284366
    284592 286097 287373 2886538 2905844 29771 2996 299780 2998212 309942
    32119 327266 328905 344590 346724 350623 36063 370026 39862 408822
    410632 428162 42846 441507 466727 475014 4959 50452 525771 56490
    596920 62722 630096 651987 707889 77067 83410 85710 88642 984906
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
