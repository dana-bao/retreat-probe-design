#!/usr/bin/env bash
#SBATCH --job-name=prisec_low
#SBATCH --partition=icelake
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=12:00:00
#SBATCH --array=0-20%5
#SBATCH --output=/rds/project/rds-FSWe0O6MFwc/users/hb676/logs/%x_%A_%a.out
#SBATCH --error=/rds/project/rds-FSWe0O6MFwc/users/hb676/logs/%x_%A_%a.err

set -euo pipefail

module purge
module use /usr/local/software/spack/spack-modules/icelake-20220710/linux-rocky8-icelake
module load samtools/1.14/gcc/v46lwk2d

# ---- species list: target_taxid  correct_taxid (organism_taxid from db_coverage.tsv) ----
SPECIES_IDS=(
    1483842  410632  56490   2886538 2905844
    370026   327266  286097  287373  346724
    229869   1553662 1875987 88642   525771
    707889   29771   187635  1875977 1293966
    1250597
)

CORRECT_TAXIDS=(
    414939   169227  370292  2741129 2905844
    370026   492533  8386    287373  866018
    229869   153785  177636  88642   1943797
    707889   29771   187635  1875977 1293975
    128015
)

TASK="${SLURM_ARRAY_TASK_ID:-0}"
N="${#SPECIES_IDS[@]}"

if [[ "$TASK" -ge "$N" ]]; then
    echo "TASK=$TASK out of range (N=$N) — exiting cleanly."
    exit 0
fi

SPECIES_ID="${SPECIES_IDS[$TASK]}"
CORRECT_TAXID="${CORRECT_TAXIDS[$TASK]}"

echo "[$(date)] TASK=$TASK  SPECIES_ID=$SPECIES_ID  CORRECT_TAXID=$CORRECT_TAXID"

BAMROOT="/rds/project/rds-FSWe0O6MFwc/users/hb676/bowtie2_bam"
SHARD_MAX=77   # inclusive; shards 00-77 only
THREADS="${SLURM_CPUS_PER_TASK:-8}"

NAMES="/rds/project/rds-FSWe0O6MFwc/users/hb676/kraken2_taxonomy/names.dmp"
NODES="/rds/project/rds-FSWe0O6MFwc/users/hb676/kraken2_taxonomy/nodes.dmp"
ACC2TAX="/rds/project/rds-FSWe0O6MFwc/users/hb676/kraken2_taxonomy/seqid2taxid.acc2tax"

OUTDIR="/rds/project/rds-FSWe0O6MFwc/users/hb676/off_target_analysis"
mkdir -p "$OUTDIR"

OUT_SUMMARY="${OUTDIR}/${SPECIES_ID}.primary_secondary_summary.txt"

# ---- skip if already done ----
if [[ -s "$OUT_SUMMARY" ]]; then
    echo "Already done: $OUT_SUMMARY — skipping."
    exit 0
fi

# ---- collect shards 00-77 ----
mapfile -t ALL_SHARDS < <(ls "$BAMROOT/$SPECIES_ID/"*.compressed.bam 2>/dev/null | sort)
if [[ "${#ALL_SHARDS[@]}" -eq 0 ]]; then
    echo "ERROR: no compressed BAMs found under $BAMROOT/$SPECIES_ID/"
    exit 1
fi
SHARDS=("${ALL_SHARDS[@]:0:$((SHARD_MAX + 1))}")
echo "Using ${#SHARDS[@]} shards (of ${#ALL_SHARDS[@]} total) piped to analysis"

samtools merge -n -@ "$THREADS" -f -O sam - "${SHARDS[@]}" \
  | python3 /rds/project/rds-FSWe0O6MFwc/users/hb676/scripts/primary_secondary_analysis.py \
      --sam           - \
      --target_taxid  "$CORRECT_TAXID" \
      --acc2tax       "$ACC2TAX" \
      --names         "$NAMES" \
      --nodes         "$NODES" \
      --out_summary   "$OUT_SUMMARY" \
      --top 20

echo "[$(date)] Done."
echo "Summary: $OUT_SUMMARY"
