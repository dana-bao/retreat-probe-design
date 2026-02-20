#!/usr/bin/env bash
#SBATCH --job-name=merge_sort_names
#SBATCH --partition=icelake
#SBATCH --cpus-per-task=20
#SBATCH --mem=240G
#SBATCH --time=12:00:00
#SBATCH --array=0-81%10
#SBATCH --output=/rds/project/rds-FSWe0O6MFwc/users/hb676/logs/%x_%A_%a.out
#SBATCH --error=/rds/project/rds-FSWe0O6MFwc/users/hb676/logs/%x_%A_%a.err

set -euo pipefail
export LC_ALL=C.UTF-8
export LANG=C.UTF-8

# ---- EDIT PATHS ----
BAMROOT="/rds/project/rds-FSWe0O6MFwc/users/hb676/bowtie2_bam"   # contains species_id/ folders
OUTROOT="/rds/project/rds-FSWe0O6MFwc/users/hb676/bowtie2_sorted" # output destination
LOGROOT="/rds/project/rds-FSWe0O6MFwc/users/hb676/logs"

THREADS="${SLURM_CPUS_PER_TASK:-20}"
THREADS="${THREADS//[^0-9]/}"
[[ -n "$THREADS" ]] || { echo "ERROR: THREADS empty"; exit 10; }

# Per-thread memory for samtools sort
SORT_MEM_PER_THREAD="10G"

mkdir -p "$OUTROOT" "$LOGROOT"

# ---- tools ----
module purge
module use /usr/local/software/spack/spack-modules/icelake-20220710/linux-rocky8-icelake
module load samtools/1.14/gcc/v46lwk2d


echo "HOST=$(hostname)"
echo "samtools=$(command -v samtools)"
samtools --version | head -n 2
echo "THREADS=$THREADS  SORT_MEM_PER_THREAD=$SORT_MEM_PER_THREAD"

# ---- pick species folder by array index ----
mapfile -t SPECIES_DIRS < <(find "$BAMROOT" -mindepth 1 -maxdepth 1 -type d | sort)
N="${#SPECIES_DIRS[@]}"
echo "Found N=$N species folders under $BAMROOT"

TASK="${SLURM_ARRAY_TASK_ID:-0}"
if [[ "$TASK" -ge "$N" ]]; then
  echo "TASK=$TASK out of range (N=$N) -> exiting cleanly."
  exit 0
fi

SPECIES_PATH="${SPECIES_DIRS[$TASK]}"
SPECIES_ID="$(basename "$SPECIES_PATH")"

echo "[$(date)] TASK=$TASK SPECIES_ID=$SPECIES_ID SPECIES_PATH=$SPECIES_PATH"

# ---- gather input BAMs ----
shopt -s nullglob
BAMS=( "$SPECIES_PATH"/*.bam )
shopt -u nullglob

if [[ "${#BAMS[@]}" -eq 0 ]]; then
  echo "No BAMs found in $SPECIES_PATH; skipping."
  exit 0
fi

# ---- outputs ----
OUTDIR="$OUTROOT/$SPECIES_ID"
mkdir -p "$OUTDIR"

SORTED="$OUTDIR/${SPECIES_ID}.merged.sorted.sam.gz"

# Optional: skip if already done and valid
if [[ -s "$SORTED" ]] && gzip -t "$SORTED" 2>/dev/null; then
  echo "Output exists and is valid ($SORTED). Skipping."
  exit 0
fi

# ---- temp ----
TMPDIR="$OUTDIR/tmp_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$TMPDIR"
echo "Temp dir: $TMPDIR"
touch "$TMPDIR/.write_test" || { echo "TMPDIR not writable: $TMPDIR"; exit 20; }
rm -f "$TMPDIR/.write_test"
df -h "$TMPDIR" || true

echo "Input BAM count: ${#BAMS[@]}"

# ---- merge and sort by name (pipe to avoid large intermediate header) ----
echo "[$(date)] Merging and sorting -> $SORTED"
{
  samtools view -H "${BAMS[0]}"
  for bam in "${BAMS[@]}"; do
    samtools view "$bam"
  done
} | samtools sort -@ "$THREADS" -m "$SORT_MEM_PER_THREAD" -n -T "$TMPDIR/sorttmp" -O SAM - | \
  gzip -1 > "$SORTED"

rm -rf "$TMPDIR" || true

# ---- validate ----
echo "[$(date)] quickcheck"
gzip -t "$SORTED" && echo "SAM.gz ok"

echo "[$(date)] Done: $SORTED"
