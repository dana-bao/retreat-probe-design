#!/usr/bin/env bash
#SBATCH --job-name=ngslca_comp
#SBATCH --partition=icelake
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=36:00:00
#SBATCH --array=0-327%10
#SBATCH --output=/rds/project/rds-FSWe0O6MFwc/users/zh384/logs/%x_%A_%a.out
#SBATCH --error=/rds/project/rds-FSWe0O6MFwc/users/zh384/logs/%x_%A_%a.err

# ngsLCA for bamsort_compressed output (*.chunk*.sorted.bam files).
# Array upper bound = 82 species * 4 chunks = 328 tasks (0-327).
# Incomplete species (< 4 chunks) are naturally absent from the file list.
# Already-finished .lca files are skipped automatically.
# Resubmit safely when more bamsort_compressed species finish.

set -euo pipefail

module purge
module use /usr/local/software/spack/spack-modules/icelake-20220710/linux-rocky8-icelake
module load samtools/1.14/gcc/v46lwk2d

SORTEDROOT="/rds/project/rds-FSWe0O6MFwc/users/hb676/bowtie2_sorted_compressed"
OUTROOT="/rds/project/rds-FSWe0O6MFwc/users/hb676/ngslca_out_compressed"

K2TAX="/rds/project/rds-FSWe0O6MFwc/users/hb676/kraken2_taxonomy"
NAMES="${K2TAX}/names.dmp"
NODES="${K2TAX}/nodes.dmp"
ACC2TAX="${K2TAX}/seqid2taxid.acc2tax"
if [[ ! -s "$ACC2TAX" ]]; then
    echo "Converting seqid2taxid.map to ngsLCA acc2tax format..."
    awk 'BEGIN{OFS="\t"}{n=split($1,a,"|"); acc=a[n]; print acc, acc, $2, 0}' "${K2TAX}/seqid2taxid.map" > "$ACC2TAX"
    echo "Done: $ACC2TAX"
fi

NGSLCA="/rds/project/rds-FSWe0O6MFwc/users/zw442/.envs/dnadamage/bin/ngsLCA"

mkdir -p "$OUTROOT"

# ---- build file list: only species with all 4 chunks present AND passing quickcheck ----
mapfile -t BAM_FILES < <(
    for d in "$SORTEDROOT"/*/; do
        sp=$(basename "$d")
        chunks=( "$d"${sp}.chunk*.sorted.bam )
        [[ ${#chunks[@]} -eq 4 ]] || continue
        all_ok=true
        for c in "${chunks[@]}"; do
            samtools quickcheck "$c" 2>/dev/null || { all_ok=false; break; }
        done
        $all_ok && printf '%s\n' "${chunks[@]}"
    done | sort
)

N="${#BAM_FILES[@]}"
echo "Found N=$N chunk files across complete species (4-chunk sets only)"

TASK="${SLURM_ARRAY_TASK_ID:-0}"
if [[ "$TASK" -ge "$N" ]]; then
    echo "TASK=$TASK out of range (N=$N) -> exiting cleanly."
    exit 0
fi

BAM="${BAM_FILES[$TASK]}"
SPECIES_ID="$(basename "$(dirname "$BAM")")"
FNAME="$(basename "$BAM")"

CHUNK=$(echo "$FNAME" | grep -oP 'chunk\d+')
LCA_LABEL="${SPECIES_ID}_${CHUNK}"

echo "[$(date)] TASK=$TASK SPECIES_ID=$SPECIES_ID LCA_LABEL=$LCA_LABEL BAM=$BAM"

# ---- skip if already done ----
LCA_OUT="$OUTROOT/${LCA_LABEL}"
if [[ -s "${LCA_OUT}.lca" ]]; then
    echo "LCA output already exists (${LCA_OUT}.lca). Skipping."
    exit 0
fi

# ---- run ngsLCA ----
echo "[$(date)] Running ngsLCA for $LCA_LABEL"
cd "$OUTROOT"
"$NGSLCA" \
    -simscorelow 0.95 \
    -simscorehigh 1.0 \
    -fix-ncbi 0 \
    -names "$NAMES" \
    -nodes "$NODES" \
    -acc2tax "$ACC2TAX" \
    -bam "$BAM" \
    -outnames "$LCA_OUT" 2>/dev/null

echo "[$(date)] Done: ${LCA_OUT}.lca"
