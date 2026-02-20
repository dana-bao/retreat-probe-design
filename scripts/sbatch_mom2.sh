#!/bin/bash
#SBATCH --job-name=kraken2_classify_corent
#SBATCH --output=logs_uniq/k2_%A_%a.out
#SBATCH --error=logs_uniq/k2_%A_%a.err
#SBATCH --array=1-10
#SBATCH --nodes=1
#SBATCH --nodelist=planex6
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=30G
#SBATCH --time=24:00:00
#SBATCH --partition=PlanEx2

# --------------------------
# Setup
# --------------------------

# Input files list (order must match array indices)
FILES=(
  129212.fasta.gz
  1753982.fasta.gz
  184104.fasta.gz
  284592.fasta.gz
  286097.fasta.gz
  29771.fasta.gz
  309942.fasta.gz
  42846.fasta.gz
  596920.fasta.gz
  85710.fasta.gz
)

# Pick the right file
FASTQ=${FILES[$SLURM_ARRAY_TASK_ID-1]}
BASENAME=$(basename "$FASTQ" .fasta.gz)

# Create directories if needed
mkdir -p results_uniq logs_uniq

echo "[$(date)] Starting Kraken2 classification on $FASTQ"
echo "Job ID: $SLURM_JOB_ID | Array Task: $SLURM_ARRAY_TASK_ID | Node: $HOSTNAME"

# --------------------------
# Run Kraken2
# --------------------------
CONFIDENCE=0.2
OUTPUT="results_uniq/${BASENAME}_task${SLURM_ARRAY_TASK_ID}.k2.$CONFIDENCE.core_nt.out"
REPORT="results_uniq/${BASENAME}_task${SLURM_ARRAY_TASK_ID}.k2.$CONFIDENCE.core_nt.report"

srun k2 classify \
  --db /dev/shm/aje39/kraken/core_nt \
  --memory-mapping \
  --threads 6 \
  --confidence $CONFIDENCE \
  --output "$OUTPUT" \
  --report "$REPORT" \
  --report-minimizer-data \
  "$FASTQ"

#  --report "$REPORT" \

echo "[$(date)] Finished $FASTQ"


# --------------------------
# Run Kraken2
# --------------------------
CONFIDENCE=0.1
OUTPUT="results_uniq/${BASENAME}_task${SLURM_ARRAY_TASK_ID}.k2.$CONFIDENCE.core_nt.out"
REPORT="results_uniq/${BASENAME}_task${SLURM_ARRAY_TASK_ID}.k2.$CONFIDENCE.core_nt.report"

srun k2 classify \
  --db /dev/shm/aje39/kraken/core_nt \
  --memory-mapping \
  --threads 6 \
  --confidence $CONFIDENCE \
  --output "$OUTPUT" \
  --report "$REPORT" \
  --report-minimizer-data \
  "$FASTQ"

#  --report "$REPORT" \

echo "[$(date)] Finished $FASTQ"


# --------------------------
# Run Kraken2
# --------------------------
CONFIDENCE=0.05
OUTPUT="results_uniq/${BASENAME}_task${SLURM_ARRAY_TASK_ID}.k2.$CONFIDENCE.core_nt.out"
REPORT="results_uniq/${BASENAME}_task${SLURM_ARRAY_TASK_ID}.k2.$CONFIDENCE.core_nt.report"

srun k2 classify \
  --db /dev/shm/aje39/kraken/core_nt \
  --memory-mapping \
  --threads 6 \
  --confidence $CONFIDENCE \
  --output "$OUTPUT" \
  --report "$REPORT" \
  --report-minimizer-data \
  "$FASTQ"

#  --report "$REPORT" \

echo "[$(date)] Finished $FASTQ"


# --------------------------
# Run Kraken2
# --------------------------
CONFIDENCE=0
OUTPUT="results_uniq/${BASENAME}_task${SLURM_ARRAY_TASK_ID}.k2.$CONFIDENCE.core_nt.out"
REPORT="results_uniq/${BASENAME}_task${SLURM_ARRAY_TASK_ID}.k2.$CONFIDENCE.core_nt.report"

srun k2 classify \
  --db /dev/shm/aje39/kraken/core_nt \
  --memory-mapping \
  --threads 6 \
  --confidence $CONFIDENCE \
  --output "$OUTPUT" \
  --report "$REPORT" \
  --report-minimizer-data \
  "$FASTQ"

#  --report "$REPORT" \

echo "[$(date)] Finished $FASTQ"
