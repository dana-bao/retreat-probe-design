#!/usr/bin/env bash
#SBATCH --job-name=lca_split
#SBATCH --partition=icelake
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=10:00:00
#SBATCH --output=/rds/project/rds-FSWe0O6MFwc/users/zh384/logs/%x_%j.out
#SBATCH --error=/rds/project/rds-FSWe0O6MFwc/users/zh384/logs/%x_%j.err

set -euo pipefail
export LC_ALL=C

# ---- EDIT ----
SPECIES_ID="${SPECIES_ID:-1293966}"
N_CHUNKS="${N_CHUNKS:-4}"
SORTEDROOT="${SORTEDROOT:-/rds/project/rds-FSWe0O6MFwc/users/hb676/bowtie2_sorted}"
LCA_OUTROOT="${LCA_OUTROOT:-/rds/project/rds-FSWe0O6MFwc/users/hb676/ngslca_out_new}"

MERGED_LCA="${LCA_OUTROOT}/${SPECIES_ID}.lca"
CHUNK_PREFIX="${SORTEDROOT}/${SPECIES_ID}/${SPECIES_ID}.merged.sorted"

echo "[$(date)] lca_split: SPECIES_ID=$SPECIES_ID N_CHUNKS=$N_CHUNKS"

# ---- check inputs ----
if [[ ! -s "$MERGED_LCA" ]]; then
    echo "ERROR: merged LCA not found or empty: $MERGED_LCA"
    exit 1
fi
for i in $(seq 0 $((N_CHUNKS - 1))); do
    cf=$(printf "%s.chunk%02d.sorted.sam.gz" "$CHUNK_PREFIX" "$i")
    if [[ ! -s "$cf" ]]; then
        echo "ERROR: SAM chunk not found: $cf"
        exit 1
    fi
done

# ---- skip if all chunk LCAs already exist ----
all_done=true
for i in $(seq 0 $((N_CHUNKS - 1))); do
    cf="${LCA_OUTROOT}/${SPECIES_ID}_chunk$(printf '%02d' $i).lca"
    [[ -s "$cf" ]] || { all_done=false; break; }
done
if $all_done; then
    echo "All $N_CHUNKS chunk LCAs already exist. Skipping."
    exit 0
fi

TMPDIR_WORK=$(mktemp -d)
trap 'rm -rf "$TMPDIR_WORK"' EXIT

# ---- Step 1: scan SAM chunks in parallel to find boundaries ----
# Each chunk (except the first) contributes one boundary = first read name of that chunk.
for i in $(seq 1 $((N_CHUNKS - 1))); do
    chunk_sam=$(printf "%s.chunk%02d.sorted.sam.gz" "$CHUNK_PREFIX" "$i")
    (
        # pipefail disabled: awk exits early after first read, sending SIGPIPE to zcat
        set +o pipefail
        first=$(zcat "$chunk_sam" | awk '!/^@/{print $1; exit}')
        if [[ -z "$first" ]]; then
            echo "ERROR: no alignment lines in $chunk_sam" >&2; exit 1
        fi
        printf '%s\n' "$first" > "$TMPDIR_WORK/boundary_${i}"
        echo "  [$(date)] first read of chunk $(printf '%02d' $i): $first" >&2
    ) &
done
wait

# Collect boundaries in order
BOUNDARY_ARGS=()
for i in $(seq 1 $((N_CHUNKS - 1))); do
    b=$(< "$TMPDIR_WORK/boundary_${i}")
    BOUNDARY_ARGS+=(-v "b${i}=${b}")
done

# ---- Step 2: generate awk split script dynamically ----
AWK_SCRIPT="$TMPDIR_WORK/split.awk"

{
    printf 'BEGIN { chunk = 0 }\n'

    # Header lines (start with tab or #) go to all chunks
    printf '/^\t|^#/ {\n'
    for i in $(seq 0 $((N_CHUNKS - 1))); do
        printf '    print > out%d\n' "$i"
    done
    printf '    next\n}\n'

    # Data lines: advance chunk when read name crosses a boundary
    printf '{\n    rname = $1\n'
    for i in $(seq 0 $((N_CHUNKS - 2))); do
        printf '    if (chunk == %d && rname >= b%d) chunk = %d\n' "$i" "$((i+1))" "$((i+1))"
    done

    # Route to output file
    for i in $(seq 0 $((N_CHUNKS - 1))); do
        if [[ $i -eq 0 ]]; then
            printf '    if (chunk == 0) print > out0\n'
        else
            printf '    else if (chunk == %d) print > out%d\n' "$i" "$i"
        fi
    done
    printf '    n[chunk]++\n}\n'

    # Summary on exit
    printf 'END {\n    for (i=0; i<%d; i++)\n        printf "[lca_split] chunk %%02d: %%d lines\\n", i, n[i] | "cat >&2"\n}\n' "$N_CHUNKS"
} > "$AWK_SCRIPT"

# Build -v out0=... out1=... arguments
OUT_ARGS=()
for i in $(seq 0 $((N_CHUNKS - 1))); do
    outfile="${LCA_OUTROOT}/${SPECIES_ID}_chunk$(printf '%02d' $i).lca"
    OUT_ARGS+=(-v "out${i}=${outfile}")
done

echo "[$(date)] Splitting $MERGED_LCA into $N_CHUNKS chunks (awk)..."
awk "${OUT_ARGS[@]}" "${BOUNDARY_ARGS[@]}" -f "$AWK_SCRIPT" "$MERGED_LCA"

echo "[$(date)] Done. Output: ${LCA_OUTROOT}/${SPECIES_ID}_chunk*.lca"
