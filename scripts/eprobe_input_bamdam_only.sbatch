#!/usr/bin/env bash
#SBATCH --job-name=eprobe_bamdam_only
#SBATCH --partition=icelake
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --array=0-6%7
#SBATCH --output=/rds/project/rds-FSWe0O6MFwc/users/hb676/logs/%x_%A_%a.out
#SBATCH --error=/rds/project/rds-FSWe0O6MFwc/users/hb676/logs/%x_%A_%a.err

set -euo pipefail

# Species whose overlap (kraken2 ∩ bamdam) is too low — use bamdam_filtered only
SPECIES_IDS=(
    56490
    346724
    370026
    2905844
    596920
    29771
    1483842
    707889
    187635
    229869
    88642
    1875987
    129212
    164537
    185980
    2998212
    139363
    32119
    328905
    1553662
    1963269
    226588
    2996
    466727
    651987
    62722
    1293966
    1534596
    1753982
    525771
    441507
    135145
    1875977
    2886538
    1886436
    83410
)

TASK="${SLURM_ARRAY_TASK_ID:-0}"
N="${#SPECIES_IDS[@]}"

if [[ "$TASK" -ge "$N" ]]; then
    echo "TASK=$TASK out of range (N=$N) — exiting cleanly."
    exit 0
fi

SID="${SPECIES_IDS[$TASK]}"
echo "[$(date)] TASK=$TASK  SPECIES_ID=$SID"

BAMDAM_FILT="/rds/project/rds-FSWe0O6MFwc/users/hb676/bamdam_filtered_compressed"
MODERN_TILED="/rds/project/rds-FSWe0O6MFwc/users/hb676/modern_tiled"
OUT_DIR="/rds/project/rds-FSWe0O6MFwc/users/hb676/eprobe_input"
SCRIPTS="/rds/project/rds-FSWe0O6MFwc/users/hb676/scripts"

python3 "${SCRIPTS}/eprobe_input_bamdam_only.py" \
    --species_id       "$SID" \
    --bamdam_filt_dir  "$BAMDAM_FILT" \
    --modern_tiled_dir "$MODERN_TILED" \
    --out_dir          "$OUT_DIR"

echo "[$(date)] Done."
