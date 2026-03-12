#!/usr/bin/env bash
#SBATCH --job-name=ngslca_modern
#SBATCH --partition=icelake
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --array=0-6
#SBATCH --output=/rds/project/rds-FSWe0O6MFwc/users/hb676/logs/%x_%A_%a.out
#SBATCH --error=/rds/project/rds-FSWe0O6MFwc/users/hb676/logs/%x_%A_%a.err

set -euo pipefail

# ---- EDIT PATHS ----
SORTEDROOT="/rds/project/rds-FSWe0O6MFwc/users/hb676/modern_bowtie2_sorted"
OUTROOT="/rds/project/rds-FSWe0O6MFwc/users/hb676/modern_ngslca_out_new"

K2TAX="/rds/project/rds-FSWe0O6MFwc/users/hb676/kraken2_taxonomy"
NAMES="${K2TAX}/names.dmp"
NODES="${K2TAX}/nodes.dmp"
ACC2TAX="${K2TAX}/seqid2taxid.acc2tax"
if [[ ! -s "$ACC2TAX" ]]; then
    echo "Converting seqid2taxid.map to ngsLCA acc2tax format..."
    awk 'BEGIN{OFS="\t"}{n=split($1,a,"|"); acc=a[n]; print acc, acc, $2, 0}' "${K2TAX}/seqid2taxid.map" > "$ACC2TAX"
    echo "Done: $ACC2TAX"
fi

NGSLCA="/rds/project/rds-FSWe0O6MFwc/users/zw442/.envs/dnadamage/bin/ngsLCA"

mkdir -p "$OUTROOT"

# ---- pick input file by array index (regular + chunk files) ----
# If a species has chunk files, skip its full-length file
mapfile -t ALL_FILES < <(find "$SORTEDROOT" \( \
    -name "*.merged.sorted.bam" \
    -o -name "*.merged.sorted.sam.gz" \
    -o -name "*.chunk*.sorted.sam.gz" \
  \) -size +100c | sort)

declare -A HAS_CHUNKS
for f in "${ALL_FILES[@]}"; do
  if [[ "$(basename "$f")" == *.chunk*.sorted.sam.gz ]]; then
    HAS_CHUNKS["$(basename "$(dirname "$f")")"]="1"
  fi
done

BAM_FILES=()
for f in "${ALL_FILES[@]}"; do
  fname="$(basename "$f")"
  species_dir="$(basename "$(dirname "$f")")"
  if [[ "$fname" != *.chunk*.sorted.sam.gz ]] && [[ -n "${HAS_CHUNKS[$species_dir]+x}" ]]; then
    echo "Skipping full-length file (chunks exist for $species_dir): $f"
    continue
  fi
  BAM_FILES+=("$f")
done

N="${#BAM_FILES[@]}"
echo "Found N=$N completed input files"

TASK="${SLURM_ARRAY_TASK_ID:-0}"
if [[ "$TASK" -ge "$N" ]]; then
  echo "TASK=$TASK out of range (N=$N) -> exiting cleanly."
  exit 0
fi

BAM="${BAM_FILES[$TASK]}"
SPECIES_ID="$(basename "$(dirname "$BAM")")"
FNAME="$(basename "$BAM")"

# For chunk files derive a unique output name per chunk
if [[ "$FNAME" == *.chunk*.sorted.sam.gz ]]; then
  CHUNK=$(echo "$FNAME" | grep -oP 'chunk\d+')
  LCA_LABEL="${SPECIES_ID}_${CHUNK}"
else
  LCA_LABEL="${SPECIES_ID}"
fi

echo "[$(date)] TASK=$TASK SPECIES_ID=$SPECIES_ID LCA_LABEL=$LCA_LABEL BAM=$BAM"

# ---- skip if already done ----
LCA_OUT="$OUTROOT/${LCA_LABEL}"
if [[ -s "${LCA_OUT}.lca" ]]; then
  echo "LCA output already exists (${LCA_OUT}.lca). Skipping."
  exit 0
fi

# ---- run ngsLCA ----
echo "[$(date)] Running ngsLCA for $LCA_LABEL"
cd "$OUTROOT"
"$NGSLCA" \
    -simscorelow 0.95 \
    -simscorehigh 1.0 \
    -fix-ncbi 0 \
    -names "$NAMES" \
    -nodes "$NODES" \
    -acc2tax "$ACC2TAX" \
    -bam "$BAM" \
    -outnames "$LCA_OUT" 2>/dev/null

echo "[$(date)] Done: ${LCA_OUT}.lca"

