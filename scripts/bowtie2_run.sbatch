#!/usr/bin/env bash
#SBATCH --job-name=bt2_run
#SBATCH --partition=icelake
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --array=0-6395%200
#SBATCH --output=/rds/project/rds-FSWe0O6MFwc/users/hb676/logs/%x_%A_%a.out
#SBATCH --error=/rds/project/rds-FSWe0O6MFwc/users/hb676/logs/%x_%A_%a.err

set -euo pipefail

export LC_ALL=C.UTF-8
export LANG=C.UTF-8

READDIR="/home/hb676/rds/rds-aegis-1-FSWe0O6MFwc/users/hb676/tiled_reads"
IDXDIR="/home/hb676/rds/rds-aegis-1-FSWe0O6MFwc/users/hb676/core_nt/bowtie2_by_vol"
BAMROOT="/rds/project/rds-FSWe0O6MFwc/users/hb676/bowtie2_bam"
READ_EXT="fasta"

mkdir -p "$BAMROOT" /rds/project/rds-FSWe0O6MFwc/users/hb676/logs

# ---- modules (MUST succeed) ----
module purge
module load bowtie/2.5.0
module load ceuadmin/samtools/1.20

echo "HOST=$(hostname)"
echo "bowtie2=$(command -v bowtie2)"
echo "samtools=$(command -v samtools)"
bowtie2 --version | head -n 2
samtools --version | head -n 2

# ---- threads: force numeric ----
THREADS="${SLURM_CPUS_PER_TASK:-16}"
THREADS="${THREADS//[^0-9]/}"
if [[ -z "$THREADS" ]]; then
  echo "ERROR: THREADS resolved to empty; SLURM_CPUS_PER_TASK='${SLURM_CPUS_PER_TASK:-unset}'" >&2
  exit 10
fi
echo "THREADS=$THREADS"

# ---- samples list ----
mapfile -t READS < <(find "$READDIR" -maxdepth 1 -type f -name "*.${READ_EXT}" | sort)
N="${#READS[@]}"
echo "Found N=$N ${READ_EXT} files under $READDIR"

if [[ "$N" -eq 0 ]]; then
  echo "ERROR: No *.${READ_EXT} files found in $READDIR" >&2
  exit 1
fi

NSHARDS=78
TASK="${SLURM_ARRAY_TASK_ID:-0}"
SAMPLE_IDX=$(( TASK / NSHARDS ))
SHARD_IDX=$(( TASK % NSHARDS ))
SHARD=$(printf "%02d" "$SHARD_IDX")

if [[ "$SAMPLE_IDX" -ge "$N" ]]; then
  echo "ERROR: SAMPLE_IDX=$SAMPLE_IDX out of range (N=$N). Fix --array range." >&2
  exit 2
fi

READFILE="${READS[$SAMPLE_IDX]}"
SAMPLE="$(basename "$READFILE" ".${READ_EXT}")"
REFBASE="${IDXDIR}/core_nt.${SHARD}"

BAMDIR="${BAMROOT}/${SAMPLE}"
mkdir -p "$BAMDIR"

ref_id="core_nt.${SHARD}"
output_bam="${BAMDIR}/${SAMPLE}_${ref_id}.bam"
logfile="${BAMDIR}/${SAMPLE}_${ref_id}.bowtie2.log"

echo "[$(date)] TASK=$TASK SAMPLE=$SAMPLE SHARD=$SHARD"
echo "READFILE=$READFILE"
echo "REFBASE=$REFBASE"
echo "output_bam=$output_bam"
echo "First line of READFILE: $(head -n 1 "$READFILE" || true)"

# index presence
ls -lh "${REFBASE}".1.bt2* | head -n 2

# ---- temp: DO NOT use /tmp ----
# Put temp under project filesystem to avoid node /tmp space issues.
TMPROOT="${BAMROOT}/_tmp_bt2"
mkdir -p "$TMPROOT"

JOBTMP="${TMPROOT}/${SLURM_JOB_ID}"
mkdir -p "$JOBTMP"

# Write a temp BAM, then validate, then atomic move.
tmp_bam="${JOBTMP}/${SAMPLE}_${ref_id}.bam.tmp"

# Optional: show available space for debugging
df -h "$TMPROOT" | head -n 2 || true

# ---- mapping ----
bowtie2 -p "$THREADS" -k 100 -x "$REFBASE" -f -U "$READFILE" --no-unal 2> "$logfile" \
| samtools view -@ "$THREADS" -b -o "$tmp_bam" -

# validate temp BAM before moving
samtools quickcheck -v "$tmp_bam"

mv -f "$tmp_bam" "$output_bam"
echo "[$(date)] Done: $output_bam"
