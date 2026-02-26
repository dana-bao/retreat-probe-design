#!/usr/bin/env bash
#SBATCH --job-name=bamdam_modern
#SBATCH --partition=icelake
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=36:00:00
#SBATCH --array=0-5
#SBATCH --output=/rds/project/rds-FSWe0O6MFwc/users/zh384/logs/%x_%A_%a.out
#SBATCH --error=/rds/project/rds-FSWe0O6MFwc/users/zh384/logs/%x_%A_%a.err

set -euo pipefail

# ---- tools ----
module purge
module use /usr/local/software/spack/spack-modules/icelake-20220710/linux-rocky8-icelake
module load samtools/1.14/gcc/v46lwk2d

# ---- PATHS ----
SORTEDROOT="/rds/project/rds-FSWe0O6MFwc/users/hb676/modern_bowtie2_sorted"
LCA_OUTROOT="/rds/project/rds-FSWe0O6MFwc/users/hb676/modern_ngslca_out"
BAMDAM_OUTROOT="/rds/project/rds-FSWe0O6MFwc/users/hb676/modern_bamdam_out"
BAMDAM="/rds/project/rds-FSWe0O6MFwc/users/zw442/.envs/dnadamage/bin/bamdam"

STRANDEDNESS="ss"
UPTO="genus"

# ---- species list (0-indexed, matches --array=0-5) ----
# 357475_modern excluded: ngsLCA still running
SPECIES_LIST=(
    111071_modern
    1250597_modern
    344590_modern
    39862_modern
    408822_modern
    42846_modern
    357475_modern
)

TASK="${SLURM_ARRAY_TASK_ID:-0}"
SPECIES_ID="${SPECIES_LIST[$TASK]}"

LCA="${LCA_OUTROOT}/${SPECIES_ID}.lca"

echo "[$(date)] TASK=$TASK  SPECIES_ID=$SPECIES_ID"
echo "LCA: $LCA"

# ---- locate input BAM ----
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

echo "BAM: $BAM"

OUTDIR="${BAMDAM_OUTROOT}/${SPECIES_ID}"
mkdir -p "$OUTDIR"

SHRUNK_BAM="${OUTDIR}/${SPECIES_ID}.shrunk.bam"
SHRUNK_LCA="${OUTDIR}/${SPECIES_ID}.shrunk.lca"

# ---- skip if already done ----
if [[ -s "$SHRUNK_BAM" ]]; then
    echo "Output already exists ($SHRUNK_BAM). Skipping."
    exit 0
fi

# ---- convert .sam.gz to .bam if needed (bamdam requires BAM) ----
INPUT_BAM="$BAM"
TMP_BAM=""
if [[ "$BAM" == *.sam.gz ]]; then
    TMP_BAM="${OUTDIR}/${SPECIES_ID}.tmp.bam"
    TMP_HDR="${OUTDIR}/${SPECIES_ID}.tmp.hdr"
    echo "[$(date)] Converting .sam.gz to BAM -> $TMP_BAM"
    samtools view -b -@ "${SLURM_CPUS_PER_TASK:-4}" -o "$TMP_BAM" "$BAM"
    echo "[$(date)] Fixing SO header tag -> queryname"
    samtools view -H "$TMP_BAM" \
        | sed '/^@HD/ { /SO:/ { s/SO:[^\t]*/SO:queryname/; b }; s/$/\tSO:queryname/ }' \
        > "$TMP_HDR"
    samtools reheader "$TMP_HDR" "$TMP_BAM" > "${TMP_BAM%.bam}.rh.bam"
    mv "${TMP_BAM%.bam}.rh.bam" "$TMP_BAM"
    rm -f "$TMP_HDR"
    INPUT_BAM="$TMP_BAM"
fi

# ---- bamdam shrink ----
echo "[$(date)] Running bamdam shrink for $SPECIES_ID"
"$BAMDAM" shrink \
    --in_bam "$INPUT_BAM" \
    --in_lca "$LCA" \
    --out_bam "$SHRUNK_BAM" \
    --out_lca "$SHRUNK_LCA" \
    --stranded "$STRANDEDNESS" \
    --upto "$UPTO"

# ---- clean up temp BAM if created ----
if [[ -n "$TMP_BAM" ]]; then
    rm -f "$TMP_BAM" "${TMP_BAM%.bam}.rh.bam" "${TMP_BAM%.bam}.hdr"
    echo "Removed temp files for $SPECIES_ID"
fi

echo "[$(date)] Done: $SHRUNK_BAM"
