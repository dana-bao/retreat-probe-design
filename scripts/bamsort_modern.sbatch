#!/usr/bin/env bash
#SBATCH --job-name=bamsort_modern
#SBATCH --partition=icelake
#SBATCH --cpus-per-task=20
#SBATCH --mem=240G
#SBATCH --time=36:00:00
#SBATCH --array=0-6
#SBATCH --output=/rds/project/rds-FSWe0O6MFwc/users/zh384/logs/%x_%A_%a.out
#SBATCH --error=/rds/project/rds-FSWe0O6MFwc/users/zh384/logs/%x_%A_%a.err

set -euo pipefail
export LC_ALL=C.UTF-8
export LANG=C.UTF-8

# ---- EDIT PATHS ----
BAMROOT="/rds/project/rds-FSWe0O6MFwc/users/hb676/modern_bowtie2_bam"
OUTROOT="/rds/project/rds-FSWe0O6MFwc/users/hb676/modern_bowtie2_sorted"

THREADS="${SLURM_CPUS_PER_TASK:-20}"
THREADS="${THREADS//[^0-9]/}"
[[ -n "$THREADS" ]] || { echo "ERROR: THREADS empty"; exit 10; }

SORT_MEM_PER_THREAD="5G"   # 10 threads x 5G = 50G, leaving ~190G for merge + OS

mkdir -p "$OUTROOT" /rds/project/rds-FSWe0O6MFwc/users/zh384/logs

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

echo "Input BAM count: ${#BAMS[@]}"

# ---- outputs ----
OUTDIR="$OUTROOT/$SPECIES_ID"
mkdir -p "$OUTDIR"

SORTED="$OUTDIR/${SPECIES_ID}.merged.sorted.sam.gz"

# Skip if valid output already exists
if [[ -s "$SORTED" ]] && gzip -t "$SORTED" 2>/dev/null && [[ $(stat -c%s "$SORTED") -gt 1048576 ]]; then
  echo "Output already exists ($SORTED). Skipping."
  exit 0
fi

# ---- temp ----
TMPDIR="$OUTDIR/tmp_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$TMPDIR"
echo "Temp dir: $TMPDIR"
touch "$TMPDIR/.write_test" || { echo "TMPDIR not writable: $TMPDIR"; exit 20; }
rm -f "$TMPDIR/.write_test"
df -h "$TMPDIR" || true

# ---- stream reads from each BAM, sort by name with GNU sort ----
echo "[$(date)] Merging and sorting -> $SORTED"
{
  samtools view -H "${BAMS[0]}"
  for bam in "${BAMS[@]}"; do
    samtools view -@ 2 "$bam"
  done | sort -k1,1 -T "$TMPDIR" --parallel="$THREADS" -S 200G
} | gzip -1 > "$SORTED"

rm -rf "$TMPDIR" || true

# ---- validate ----
echo "[$(date)] quickcheck"
gzip -t "$SORTED" && echo "SAM.gz ok"

echo "[$(date)] Done: $SORTED"
