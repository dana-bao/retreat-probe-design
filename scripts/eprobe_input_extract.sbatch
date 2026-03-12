#!/usr/bin/env bash
#SBATCH --job-name=eprobe_input
#SBATCH --partition=icelake
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=36:00:00
#SBATCH --array=0-79%20
#SBATCH --output=/rds/project/rds-FSWe0O6MFwc/users/zh384/logs/%x_%A_%a.out
#SBATCH --error=/rds/project/rds-FSWe0O6MFwc/users/zh384/logs/%x_%A_%a.err

set -euo pipefail

KRAKEN_FILT="/rds/project/rds-FSWe0O6MFwc/users/hb676/kraken2_filtered"
BAMDAM_FILT="/rds/project/rds-FSWe0O6MFwc/users/hb676/bamdam_filtered_compressed"
MODERN_TILED="/rds/project/rds-FSWe0O6MFwc/users/hb676/modern_tiled"
OUT_DIR="/rds/project/rds-FSWe0O6MFwc/users/hb676/eprobe_input"

mapfile -t SPECIES_LIST < <(
    ls "$BAMDAM_FILT"/*.shrunk.correct_genus.lca 2>/dev/null \
    | xargs -n1 basename \
    | sed 's/\.shrunk\.correct_genus\.lca//'
)

TASK="${SLURM_ARRAY_TASK_ID:-0}"
N="${#SPECIES_LIST[@]}"

if [[ "$TASK" -ge "$N" ]]; then
    echo "TASK=$TASK out of range (N=$N) — exiting cleanly."
    exit 0
fi

SID="${SPECIES_LIST[$TASK]}"
echo "[$(date)] TASK=$TASK  SPECIES_ID=$SID"

python3 /rds/project/rds-FSWe0O6MFwc/users/hb676/scripts/eprobe_input_extract.py \
    --species_id      "$SID" \
    --kraken_filt_dir "$KRAKEN_FILT" \
    --bamdam_filt_dir "$BAMDAM_FILT" \
    --modern_tiled_dir "$MODERN_TILED" \
    --out_dir         "$OUT_DIR"

echo "[$(date)] Done."
