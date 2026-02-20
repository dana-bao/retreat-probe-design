#!/usr/bin/env bash
#SBATCH --job-name=ngslca
#SBATCH --partition=icelake
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=36:00:00
#SBATCH --array=0-43%10
#SBATCH --output=/rds/project/rds-FSWe0O6MFwc/users/hb676/logs/%x_%A_%a.out
#SBATCH --error=/rds/project/rds-FSWe0O6MFwc/users/hb676/logs/%x_%A_%a.err

set -euo pipefail

# ---- EDIT PATHS ----
SORTEDROOT="/rds/project/rds-FSWe0O6MFwc/users/hb676/bowtie2_sorted"
OUTROOT="/rds/project/rds-FSWe0O6MFwc/users/hb676/ngslca_out"       # TODO: set output dir

NAMES="/rds/project/rds-FSWe0O6MFwc/users/hb676/taxonomy/names.dmp"
NODES="/rds/project/rds-FSWe0O6MFwc/users/hb676/taxonomy/nodes.dmp"
ACC2TAX="/rds/project/rds-FSWe0O6MFwc/users/hb676/taxonomy/nucl_gb.accession2taxid"

NGSLCA="/rds/project/rds-FSWe0O6MFwc/users/zw442/.envs/dnadamage/bin/ngsLCA"

mkdir -p "$OUTROOT"

# ---- pick species by array index ----
mapfile -t BAM_FILES < <(find "$SORTEDROOT" \( -name "*.merged.sorted.bam" -o -name "*.merged.sorted.sam.gz" \) -size +100c | sort)
N="${#BAM_FILES[@]}"
echo "Found N=$N completed input files"

TASK="${SLURM_ARRAY_TASK_ID:-0}"
if [[ "$TASK" -ge "$N" ]]; then
  echo "TASK=$TASK out of range (N=$N) -> exiting cleanly."
  exit 0
fi

BAM="${BAM_FILES[$TASK]}"
SPECIES_ID="$(basename "$(dirname "$BAM")")"

echo "[$(date)] TASK=$TASK SPECIES_ID=$SPECIES_ID BAM=$BAM"

# ---- skip if already done ----
LCA_OUT="$OUTROOT/${SPECIES_ID}"
if [[ -s "${LCA_OUT}.lca" ]]; then
  echo "LCA output already exists (${LCA_OUT}.lca). Skipping."
  exit 0
fi

# ---- run ngsLCA ----
echo "[$(date)] Running ngsLCA for $SPECIES_ID"
cd "$OUTROOT"
"$NGSLCA" \
    -simscorelow 0.95 \
    -simscorehigh 1.0 \
    -fix-ncbi 0 \
    -names "$NAMES" \
    -nodes "$NODES" \
    -acc2tax "$ACC2TAX" \
    -bam "$BAM" \
    -outnames "$LCA_OUT"

echo "[$(date)] Done: ${LCA_OUT}.lca"
