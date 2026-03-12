#!/bin/bash
#SBATCH -p icelake
#SBATCH -c 32
#SBATCH --mem=110G
#SBATCH -t 12:00:00
#SBATCH -J bt2_core_nt
#SBATCH --array=0-77%6
#SBATCH -o bt2_core_nt.%A_%a.out
#SBATCH -e bt2_core_nt.%A_%a.err

set -euo pipefail

module purge
module load bowtie/2.5.0

BASE=/rds/project/rds-FSWe0O6MFwc/users/hb676/core_nt
FADIR="${BASE}/fasta_by_vol"
IDXDIR="${BASE}/bowtie2_by_vol"
mkdir -p "$IDXDIR"

vol=$(printf "%02d" "${SLURM_ARRAY_TASK_ID}")
prefix="${IDXDIR}/core_nt.${vol}"

# Accept either compressed or plain FASTA
if [[ -s "${FADIR}/core_nt.${vol}.fasta.gz" ]]; then
  fa_gz="${FADIR}/core_nt.${vol}.fasta.gz"
  INPUT_COMPRESSED=1
elif [[ -s "${FADIR}/core_nt.${vol}.fasta" ]]; then
  fa_gz="${FADIR}/core_nt.${vol}.fasta"
  INPUT_COMPRESSED=0
else
  echo "ERROR: No input FASTA found for vol=$vol in $FADIR" >&2
  exit 2
fi

echo "[$(date)] vol=$vol host=$(hostname) cpus=${SLURM_CPUS_PER_TASK:-32} mem=110G"
echo "[$(date)] input=$fa_gz"
echo "[$(date)] prefix=$prefix"

# Skip if already built
if ls "${prefix}".*.bt2 1>/dev/null 2>&1 || ls "${prefix}".*.bt2l 1>/dev/null 2>&1; then
  echo "[$(date)] Index exists; skipping."
  exit 0
fi

# Prefer node-local scratch if available
WORKROOT="/rds/project/rds-FSWe0O6MFwc/users/hb676/core_nt/tmp_bt2"
WORKDIR="${WORKROOT}/bt2_core_nt_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$WORKROOT"
mkdir -p "$WORKDIR"

fa="${WORKDIR}/core_nt.${vol}.fasta"

echo "[$(date)] workdir=$WORKDIR"
if [[ "$INPUT_COMPRESSED" -eq 1 ]]; then
  echo "[$(date)] decompressing to $fa ..."
  gzip -dc "$fa_gz" > "$fa"
else
  echo "[$(date)] input already uncompressed, copying to workdir ..."
  cp "$fa_gz" "$fa"
fi

# Sanity check (fast)
head -n 2 "$fa" | sed 's/.*/[$(date)] &/'

echo "[$(date)] starting bowtie2-build..."
bowtie2-build --large-index --threads "${SLURM_CPUS_PER_TASK:-32}" "$fa" "$prefix"

echo "[$(date)] finished bowtie2-build"

# Verify outputs
if ls "${prefix}".*.bt2l 1>/dev/null 2>&1 || ls "${prefix}".*.bt2 1>/dev/null 2>&1; then
  echo "[$(date)] OK: index files created."
else
  echo "ERROR: No index files found for prefix $prefix" >&2
  exit 3
fi

rm -rf "$WORKDIR"
echo "[$(date)] done."

