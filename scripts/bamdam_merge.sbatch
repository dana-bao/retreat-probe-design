#!/usr/bin/env bash
#SBATCH --job-name=bamdam_merge
#SBATCH --partition=icelake
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=4:00:00
#SBATCH --output=/rds/project/rds-FSWe0O6MFwc/users/hb676/logs/%x_%j.out
#SBATCH --error=/rds/project/rds-FSWe0O6MFwc/users/hb676/logs/%x_%j.err

set -euo pipefail

module purge
module use /usr/local/software/spack/spack-modules/icelake-20220710/linux-rocky8-icelake
module load samtools/1.14/gcc/v46lwk2d

BAMDAM_OUTROOT="/rds/project/rds-FSWe0O6MFwc/users/hb676/bamdam_out_new"
LCA_OUTROOT="/rds/project/rds-FSWe0O6MFwc/users/hb676/ngslca_out_new"

# Find all species that have chunk shrunk outputs
mapfile -t CHUNK_SPECIES < <(
    find "$BAMDAM_OUTROOT" -name "*_chunk*.shrunk.bam" -size +0c \
        | sed 's|.*/||;s|_chunk[0-9]*\.shrunk\.bam$||' \
        | sort -u
)

echo "Found ${#CHUNK_SPECIES[@]} species with chunk outputs to merge"

for SPECIES_ID in "${CHUNK_SPECIES[@]}"; do
    OUTDIR="${BAMDAM_OUTROOT}/${SPECIES_ID}"
    SHRUNK_BAM="${OUTDIR}/${SPECIES_ID}.shrunk.bam"
    SHRUNK_LCA="${OUTDIR}/${SPECIES_ID}.shrunk.lca"

    if [[ -s "$SHRUNK_BAM" ]]; then
        echo "Already merged ($SHRUNK_BAM). Skipping."
        continue
    fi

    # Verify all expected chunks are present
    N_EXPECTED=$(find "$LCA_OUTROOT" -maxdepth 1 \
        -name "${SPECIES_ID}_chunk*.lca" -size +0c | wc -l)
    mapfile -t CHUNK_BAMS < <(
        find "$OUTDIR" -name "${SPECIES_ID}_chunk*.shrunk.bam" | sort
    )
    mapfile -t CHUNK_LCAS < <(
        find "$OUTDIR" -name "${SPECIES_ID}_chunk*.shrunk.lca" | sort
    )

    if [[ "${#CHUNK_BAMS[@]}" -lt "$N_EXPECTED" ]]; then
        echo "WARNING: $SPECIES_ID has ${#CHUNK_BAMS[@]}/${N_EXPECTED} chunk BAMs — skipping (not all chunks done yet)."
        continue
    fi

    echo "[$(date)] Merging ${#CHUNK_BAMS[@]} chunks for $SPECIES_ID"
    samtools merge -f "$SHRUNK_BAM" "${CHUNK_BAMS[@]}"
    cat "${CHUNK_LCAS[@]}" > "$SHRUNK_LCA"

    # Clean up chunk intermediates
    rm -f "${CHUNK_BAMS[@]}" "${CHUNK_LCAS[@]}"
    echo "[$(date)] Done: $SHRUNK_BAM"
done
