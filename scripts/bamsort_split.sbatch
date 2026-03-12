#!/usr/bin/env bash
#SBATCH --job-name=bamsort_split
#SBATCH --partition=icelake
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --output=/rds/project/rds-FSWe0O6MFwc/users/hb676/logs/%x_%j.out
#SBATCH --error=/rds/project/rds-FSWe0O6MFwc/users/hb676/logs/%x_%j.err

set -euo pipefail
export LC_ALL=C.UTF-8
export LANG=C.UTF-8

# ---- EDIT ----
SPECIES_ID="${SPECIES_ID:-134639}"
SORTEDROOT="${SORTEDROOT:-/rds/project/rds-FSWe0O6MFwc/users/hb676/bowtie2_sorted}"
N_CHUNKS=4
INFILE="$SORTEDROOT/$SPECIES_ID/${SPECIES_ID}.merged.sorted.sam.gz"
CHUNK_PREFIX="$SORTEDROOT/$SPECIES_ID/${SPECIES_ID}.merged.sorted"

if [[ ! -s "$INFILE" ]]; then
    echo "ERROR: input file not found or empty: $INFILE"
    exit 1
fi

# ---- skip if all chunks already exist and are valid ----
all_done=true
for i in $(seq 0 $((N_CHUNKS - 1))); do
    cf=$(printf "%s.chunk%02d.sorted.sam.gz" "$CHUNK_PREFIX" "$i")
    if ! [[ -s "$cf" ]] || ! gzip -t "$cf" 2>/dev/null; then
        all_done=false
        break
    fi
done
if $all_done; then
    echo "All $N_CHUNKS chunks already exist and are valid. Skipping."
    exit 0
fi

TMPDIR=$(mktemp -d "$SORTEDROOT/$SPECIES_ID/tmp_split_XXXXXX")
trap 'rm -rf "$TMPDIR"' EXIT

# ---- write inline splitter script ----
SPLIT_PY="$TMPDIR/sam_split.py"
cat > "$SPLIT_PY" << 'PYEOF'
"""Split a queryname-sorted SAM stream into N chunks at read-name boundaries."""
import argparse, subprocess, sys

ap = argparse.ArgumentParser()
ap.add_argument("--n_chunks", type=int, required=True)
ap.add_argument("--out_prefix", required=True)
ap.add_argument("--total_aln_lines", type=int, required=True)
args = ap.parse_args()

n = args.n_chunks
target = max(1, args.total_aln_lines // n)

header_lines = []
first_aln = None
for raw in sys.stdin.buffer:
    if raw.startswith(b"@"):
        header_lines.append(raw)
    else:
        first_aln = raw
        break

if first_aln is None:
    print("[sam_split] WARNING: no alignment lines found", file=sys.stderr)
    sys.exit(0)

header_bytes = b"".join(header_lines)

def open_chunk(i):
    path = f"{args.out_prefix}.chunk{i:02d}.sorted.sam.gz"
    proc = subprocess.Popen(["gzip", "-1"], stdin=subprocess.PIPE, stdout=open(path, "wb"))
    proc.stdin.write(header_bytes)
    return proc, path

chunk_idx = 0
proc, path = open_chunk(0)
print(f"[sam_split] chunk {chunk_idx:02d} -> {path}", file=sys.stderr, flush=True)

lines_in_chunk = 0
prev_name = None
line = first_aln

while line:
    name = line.split(b"\t", 1)[0]
    if chunk_idx < n - 1 and lines_in_chunk >= target and prev_name and name != prev_name:
        proc.stdin.close()
        proc.wait()
        chunk_idx += 1
        lines_in_chunk = 0
        proc, path = open_chunk(chunk_idx)
        print(f"[sam_split] chunk {chunk_idx:02d} -> {path}", file=sys.stderr, flush=True)
    proc.stdin.write(line)
    lines_in_chunk += 1
    prev_name = name
    line = sys.stdin.buffer.readline()

proc.stdin.close()
proc.wait()
print(f"[sam_split] done — wrote {chunk_idx + 1} chunk(s)", file=sys.stderr, flush=True)
PYEOF

# ---- count alignment lines ----
echo "[$(date)] Counting alignment lines in $INFILE ..."
TOTAL_LINES=$(zcat "$INFILE" | grep -cv "^@")
echo "  Total alignment lines: $TOTAL_LINES"

# ---- split ----
echo "[$(date)] Splitting into $N_CHUNKS chunks -> ${CHUNK_PREFIX}.chunk*.sorted.sam.gz"
zcat "$INFILE" | python3 "$SPLIT_PY" \
    --n_chunks        "$N_CHUNKS" \
    --out_prefix      "$CHUNK_PREFIX" \
    --total_aln_lines "$TOTAL_LINES"

# ---- validate ----
echo "[$(date)] Validating chunks..."
for i in $(seq 0 $((N_CHUNKS - 1))); do
    cf=$(printf "%s.chunk%02d.sorted.sam.gz" "$CHUNK_PREFIX" "$i")
    if [[ -s "$cf" ]]; then
        gzip -t "$cf" && echo "  chunk $i ok: $cf"
    else
        echo "  WARNING: chunk $i missing or empty: $cf"
    fi
done
echo "[$(date)] Done."
