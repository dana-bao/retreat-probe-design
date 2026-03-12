#!/bin/bash
#SBATCH --job-name=ncbi_dl
#SBATCH --partition=icelake
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=/home/%u/rds/hpc-work/probe_design/logs/ncbi_dl-%j.out
#SBATCH --error=/home/%u/rds/hpc-work/probe_design/logs/ncbi_dl-%j.err

set -euo pipefail
export PATH=$HOME/bin:$PATH

BASE=$HOME/rds/hpc-work/probe_design
OUT=$BASE/data/assemblies
LIST=$BASE/metadata/assembly_accession.txt
LOGDIR=$BASE/logs

cd "$OUT"

datasets download genome accession --inputfile "$LIST" --dehydrated --filename assemblies.dehydrated.zip

unzip -t assemblies.dehydrated.zip >/dev/null
unzip -q assemblies.dehydrated.zip

datasets rehydrate --directory .