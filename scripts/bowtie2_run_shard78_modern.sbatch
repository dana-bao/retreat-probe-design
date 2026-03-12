#!/usr/bin/env bash
#SBATCH --job-name=bt2_shard787980_mod
#SBATCH --partition=icelake
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --time=12:00:00
#SBATCH --array=0-14%15
#SBATCH --output=/rds/project/rds-FSWe0O6MFwc/users/zh384/logs/%x_%A_%a.out
#SBATCH --error=/rds/project/rds-FSWe0O6MFwc/users/zh384/logs/%x_%A_%a.err

# Maps modern tiled reads against shards 78, 79, and 80.
# Array layout: task = sample_idx * 3 + shard_offset
#   shard_offset 0 -> shard 78, shard_offset 1 -> shard 79, shard_offset 2 -> shard 80
# Total tasks = N_samples * 3 (adjust --array if N_samples != 82)

set -euo pipefail

export LC_ALL=C.UTF-8
export LANG=C.UTF-8

READDIR="/rds/project/rds-FSWe0O6MFwc/users/hb676/modern_tiled"
IDXDIR="/rds/project/rds-FSWe0O6MFwc/users/hb676/core_nt/bowtie2_by_vol"
BAMROOT="/rds/project/rds-FSWe0O6MFwc/users/hb676/modern_bowtie2_bam"
READ_EXT="fasta"
N_SHARDS=3   # shards 78, 79, and 80

SPECIES_LIST=(
    344590_modern
    39862_modern
    357475_modern
    408822_modern
    42846_modern
)

mkdir -p "$BAMROOT" /rds/project/rds-FSWe0O6MFwc/users/hb676/logs

module purge
module load bowtie/2.5.0
module load ceuadmin/samtools/1.20

echo "HOST=$(hostname)"
bowtie2 --version | head -n 2
samtools --version | head -n 2

THREADS="${SLURM_CPUS_PER_TASK:-16}"
THREADS="${THREADS//[^0-9]/}"
[[ -n "$THREADS" ]] || { echo "ERROR: THREADS empty" >&2; exit 10; }
echo "THREADS=$THREADS"

# Decode task -> (sample, shard)
TASK="${SLURM_ARRAY_TASK_ID:-0}"
SAMPLE_IDX=$(( TASK / N_SHARDS ))
SHARD=$(( 78 + TASK % N_SHARDS ))

N="${#SPECIES_LIST[@]}"
if [[ "$SAMPLE_IDX" -ge "$N" ]]; then
    echo "SAMPLE_IDX=$SAMPLE_IDX out of range (N=$N) — exiting cleanly."
    exit 0
fi

REFBASE="${IDXDIR}/core_nt.${SHARD}"
if ! ls "${REFBASE}".1.bt2* 1>/dev/null 2>&1 && ! ls "${REFBASE}".1.bt2l 1>/dev/null 2>&1; then
    echo "ERROR: Index not found at $REFBASE" >&2
    exit 2
fi

SAMPLE="${SPECIES_LIST[$SAMPLE_IDX]}"
READFILE="${READDIR}/${SAMPLE}.${READ_EXT}"
if [[ ! -f "$READFILE" ]]; then
    echo "ERROR: read file not found: $READFILE" >&2
    exit 1
fi
BAMDIR="${BAMROOT}/${SAMPLE}"
mkdir -p "$BAMDIR"

output_bam="${BAMDIR}/${SAMPLE}_core_nt.${SHARD}.bam"
logfile="${BAMDIR}/${SAMPLE}_core_nt.${SHARD}.bowtie2.log"

echo "[$(date)] TASK=$TASK  SAMPLE=$SAMPLE  SHARD=$SHARD"
echo "READFILE=$READFILE"
echo "output_bam=$output_bam"

if [[ -s "$output_bam" ]]; then
    echo "Output exists; skipping."
    exit 0
fi

TMPROOT="${BAMROOT}/_tmp_bt2"
mkdir -p "$TMPROOT"
JOBTMP="${TMPROOT}/${SLURM_JOB_ID}_${TASK}"
mkdir -p "$JOBTMP"
tmp_bam="${JOBTMP}/${SAMPLE}_core_nt.${SHARD}.bam.tmp"

bowtie2 -p "$THREADS" -k 100 -x "$REFBASE" -f -U "$READFILE" --no-unal 2> "$logfile" \
    | samtools view -@ "$THREADS" -b -o "$tmp_bam" -

samtools quickcheck -v "$tmp_bam"
mv -f "$tmp_bam" "$output_bam"

rm -rf "$JOBTMP" || true
echo "[$(date)] Done: $output_bam"
