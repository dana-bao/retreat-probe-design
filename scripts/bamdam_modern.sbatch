#!/usr/bin/env bash
#SBATCH --job-name=bamdam_modern
#SBATCH --partition=icelake
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=36:00:00
#SBATCH --array=0-9
#SBATCH --output=/rds/project/rds-FSWe0O6MFwc/users/zh384/logs/%x_%A_%a.out
#SBATCH --error=/rds/project/rds-FSWe0O6MFwc/users/zh384/logs/%x_%A_%a.err

set -euo pipefail

# ---- tools ----
module purge
module use /usr/local/software/spack/spack-modules/icelake-20220710/linux-rocky8-icelake
module load samtools/1.14/gcc/v46lwk2d

# ---- PATHS ----
SORTEDROOT="/rds/project/rds-FSWe0O6MFwc/users/hb676/modern_bowtie2_sorted"
LCA_OUTROOT="/rds/project/rds-FSWe0O6MFwc/users/hb676/modern_ngslca_out_new"
BAMDAM_OUTROOT="/rds/project/rds-FSWe0O6MFwc/users/hb676/modern_bamdam_out_new"
BAMDAM="/rds/project/rds-FSWe0O6MFwc/users/zw442/.envs/dnadamage/bin/bamdam"

STRANDEDNESS="ss"
UPTO="genus"

# ---- hardcoded item list (6 full-LCA species + 4 chunks for 357475_modern) ----
ITEM_LIST=(
    111071_modern
    1250597_modern
    344590_modern
    39862_modern
    408822_modern
    42846_modern
    357475_modern_chunk00
    357475_modern_chunk01
    357475_modern_chunk02
    357475_modern_chunk03
)

TASK="${SLURM_ARRAY_TASK_ID:-0}"
N="${#ITEM_LIST[@]}"
if [[ "$TASK" -ge "$N" ]]; then
    echo "TASK=$TASK out of range (N=$N); exiting."
    exit 0
fi

ITEM="${ITEM_LIST[$TASK]}"
echo "[$(date)] TASK=$TASK  ITEM=$ITEM"

# ---- helper: convert SAM.gz to queryname-sorted BAM ----
CLEANUP_BAMS=()
convert_if_needed() {
    local label="$1"
    local src="$2"
    if [[ "$src" == *.sam.gz ]]; then
        local tmp="${OUTDIR}/${label}.tmp.bam"
        local hdr="${OUTDIR}/${label}.tmp.hdr"
        echo "[$(date)] Converting $label SAM.gz -> BAM"
        samtools view -b -@ "${SLURM_CPUS_PER_TASK:-4}" -o "$tmp" "$src"
        samtools view -H "$tmp" \
            | sed '/^@HD/ { /SO:/ { s/SO:[^\t]*/SO:queryname/; b }; s/$/\tSO:queryname/ }' \
            > "$hdr"
        samtools reheader "$hdr" "$tmp" > "${tmp%.bam}.rh.bam"
        mv "${tmp%.bam}.rh.bam" "$tmp"
        rm -f "$hdr"
        CONVERTED_BAM="$tmp"
        CLEANUP_BAMS+=("$tmp")
    else
        CONVERTED_BAM="$src"
    fi
}

# ============================================================
# CHUNK ITEM: process a single chunk independently
# ============================================================
if [[ "$ITEM" == *_chunk* ]]; then
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

    DATA_LINES=$(grep -cv '^#' "$CHUNK_LCA" 2>/dev/null || true)
    if [[ "$DATA_LINES" -eq 0 ]]; then
        echo "WARNING: $CHUNK_LCA has no data lines — creating empty outputs."
        samtools view -b -H \
            "$(find "${SORTEDROOT}/${SPECIES_ID}" \( \
                -name "${SPECIES_ID}.merged.sorted.${CHUNK_TAG}.sorted.sam.gz" \
                -o -name "${SPECIES_ID}.${CHUNK_TAG}.sorted.bam" \) \
            -size +0c | head -1)" > "$CHUNK_SHRUNK_BAM"
        touch "$CHUNK_SHRUNK_LCA"
        echo "[$(date)] Done (empty): $CHUNK_SHRUNK_BAM"
        exit 0
    fi

    CHUNK_SRC=""
    for candidate in \
        "${SORTEDROOT}/${SPECIES_ID}/${SPECIES_ID}.merged.sorted.${CHUNK_TAG}.sorted.sam.gz" \
        "${SORTEDROOT}/${SPECIES_ID}/${SPECIES_ID}.${CHUNK_TAG}.sorted.sam.gz" \
        "${SORTEDROOT}/${SPECIES_ID}/${SPECIES_ID}.${CHUNK_TAG}.sorted.bam"; do
        if [[ -s "$candidate" ]]; then
            CHUNK_SRC="$candidate"
            break
        fi
    done
    if [[ -z "$CHUNK_SRC" ]]; then
        echo "ERROR: no SAM/BAM found for $ITEM under $SORTEDROOT"
        exit 1
    fi

    convert_if_needed "$ITEM" "$CHUNK_SRC"

    echo "[$(date)] bamdam shrink: $ITEM"
    "$BAMDAM" shrink \
        --in_bam  "$CONVERTED_BAM" \
        --in_lca  "$CHUNK_LCA" \
        --out_bam "$CHUNK_SHRUNK_BAM" \
        --out_lca "$CHUNK_SHRUNK_LCA" \
        --stranded "$STRANDEDNESS" \
        --upto "$UPTO"

    [[ ${#CLEANUP_BAMS[@]} -gt 0 ]] && rm -f "${CLEANUP_BAMS[@]}"
    echo "[$(date)] Done: $CHUNK_SHRUNK_BAM"

# ============================================================
# FULL ITEM: single merged BAM + single LCA
# ============================================================
else
    SPECIES_ID="$ITEM"
    OUTDIR="${BAMDAM_OUTROOT}/${SPECIES_ID}"
    mkdir -p "$OUTDIR"

    SHRUNK_BAM="${OUTDIR}/${SPECIES_ID}.shrunk.bam"
    SHRUNK_LCA="${OUTDIR}/${SPECIES_ID}.shrunk.lca"

    if [[ -s "$SHRUNK_BAM" ]]; then
        echo "Already done ($SHRUNK_BAM). Skipping."
        exit 0
    fi

    BAM=""
    for ext in merged.sorted.sam.gz merged.sorted.bam; do
        candidate="${SORTEDROOT}/${SPECIES_ID}/${SPECIES_ID}.${ext}"
        if [[ -s "$candidate" ]]; then
            BAM="$candidate"
            break
        fi
    done
    if [[ -z "$BAM" ]]; then
        echo "ERROR: no input BAM/SAM.gz found for $SPECIES_ID under $SORTEDROOT"
        exit 1
    fi

    LCA="${LCA_OUTROOT}/${SPECIES_ID}.lca"
    if [[ ! -s "$LCA" ]]; then
        echo "ERROR: LCA file not found or empty: $LCA"
        exit 1
    fi

    echo "BAM: $BAM"
    echo "LCA: $LCA"

    convert_if_needed "$SPECIES_ID" "$BAM"

    echo "[$(date)] bamdam shrink: $SPECIES_ID"
    "$BAMDAM" shrink \
        --in_bam  "$CONVERTED_BAM" \
        --in_lca  "$LCA" \
        --out_bam "$SHRUNK_BAM" \
        --out_lca "$SHRUNK_LCA" \
        --stranded "$STRANDEDNESS" \
        --upto "$UPTO"

    [[ ${#CLEANUP_BAMS[@]} -gt 0 ]] && rm -f "${CLEANUP_BAMS[@]}"
    echo "[$(date)] Done: $SHRUNK_BAM"
fi
