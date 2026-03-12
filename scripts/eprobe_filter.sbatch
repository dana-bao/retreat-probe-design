#!/usr/bin/env bash
#SBATCH --job-name=eprobe_filter
#SBATCH --partition=icelake
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --array=0-200%20
#SBATCH --output=/rds/project/rds-FSWe0O6MFwc/users/hb676/logs/%x_%A_%a.out
#SBATCH --error=/rds/project/rds-FSWe0O6MFwc/users/hb676/logs/%x_%A_%a.err

set -euo pipefail

module purge
module load python/3.11.0-icl
export PATH="/home/hb676/.local/bin:$PATH"

# ---- paths ----
FASTA_DIR="/rds/project/rds-FSWe0O6MFwc/users/hb676/eprobe_input"
OUT_DIR="/rds/project/rds-FSWe0O6MFwc/users/hb676/eprobe_filtered"

mkdir -p "$OUT_DIR"

# ---- build list: priority 1) *_bamdam_sub300k.fasta  2) *_sub300k.fasta  3) *_bamdam.fasta  4) plain *.fasta ----
declare -A HAS_BAMDAM_SUB HAS_SUB HAS_BAMDAM

while IFS= read -r f; do
    sid="$(basename "$f" _bamdam_sub300k.fasta)"
    HAS_BAMDAM_SUB["$sid"]=1
done < <(find "$FASTA_DIR" -maxdepth 1 -name "*_bamdam_sub300k.fasta" -size +0c)

while IFS= read -r f; do
    sid="$(basename "$f" _sub300k.fasta)"
    HAS_SUB["$sid"]=1
done < <(find "$FASTA_DIR" -maxdepth 1 -name "*_sub300k.fasta" -size +0c \
             | grep -v '_bamdam_sub300k\.fasta$')

while IFS= read -r f; do
    sid="$(basename "$f" _bamdam.fasta)"
    HAS_BAMDAM["$sid"]=1
done < <(find "$FASTA_DIR" -maxdepth 1 -name "*_bamdam.fasta" -size +0c)

mapfile -t FASTA_LIST < <(
    {
        # Tier 1: bamdam_sub300k
        find "$FASTA_DIR" -maxdepth 1 -name "*_bamdam_sub300k.fasta" -size +0c
        # Tier 2: plain sub300k (no bamdam_sub300k for same sid)
        while IFS= read -r f; do
            sid="$(basename "$f" _sub300k.fasta)"
            [[ -z "${HAS_BAMDAM_SUB[$sid]+x}" ]] && echo "$f"
        done < <(find "$FASTA_DIR" -maxdepth 1 -name "*_sub300k.fasta" -size +0c \
                     | grep -v '_bamdam_sub300k\.fasta$')
        # Tier 3: bamdam (no bamdam_sub300k and no sub300k for same sid)
        while IFS= read -r f; do
            sid="$(basename "$f" _bamdam.fasta)"
            [[ -z "${HAS_BAMDAM_SUB[$sid]+x}" ]] && [[ -z "${HAS_SUB[$sid]+x}" ]] && echo "$f"
        done < <(find "$FASTA_DIR" -maxdepth 1 -name "*_bamdam.fasta" -size +0c)
        # Tier 4: plain (no bamdam_sub300k, no sub300k, no bamdam for same sid)
        while IFS= read -r f; do
            sid="$(basename "$f" .fasta)"
            [[ -z "${HAS_BAMDAM_SUB[$sid]+x}" ]] && [[ -z "${HAS_SUB[$sid]+x}" ]] && [[ -z "${HAS_BAMDAM[$sid]+x}" ]] && echo "$f"
        done < <(find "$FASTA_DIR" -maxdepth 1 -name "*.fasta" -size +0c \
                     | grep -v '_bamdam\.fasta$' | grep -v '_sub300k\.fasta$')
    } | sort
)

N="${#FASTA_LIST[@]}"
echo "Found N=$N non-empty FASTA files"

TASK="${SLURM_ARRAY_TASK_ID:-0}"
if [[ "$TASK" -ge "$N" ]]; then
    echo "TASK=$TASK out of range (N=$N); exiting."
    exit 0
fi

FASTA="${FASTA_LIST[$TASK]}"
SPECIES_ID="$(basename "$FASTA" .fasta)"
OUT_PREFIX="${OUT_DIR}/${SPECIES_ID}"

echo "[$(date)] TASK=$TASK  SPECIES_ID=$SPECIES_ID"

# ---- skip if already done ----
if [[ -s "${OUT_PREFIX}.filtered.fa" ]]; then
    echo "Output already exists (${OUT_PREFIX}.filtered.fa). Skipping."
    exit 0
fi

# ---- run eprobe biophysical filter ----
echo "[$(date)] Running eprobe util filter on $FASTA"
eprobe util filter \
    --input  "$FASTA" \
    --output "$OUT_PREFIX" \
    --step   biophysical \
    --threads "${SLURM_CPUS_PER_TASK:-4}"

echo "[$(date)] Done: ${OUT_PREFIX}.filtered.fa"
echo "Passed: $(grep -c '^>' "${OUT_PREFIX}.filtered.fa" || true)"
