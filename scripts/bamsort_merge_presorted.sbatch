#!/usr/bin/env bash
#SBATCH --job-name=merge_presorted
#SBATCH --partition=icelake
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=36:00:00
#SBATCH --output=/rds/project/rds-FSWe0O6MFwc/users/zh384/logs/%x_%j.out
#SBATCH --error=/rds/project/rds-FSWe0O6MFwc/users/zh384/logs/%x_%j.err

set -euo pipefail
export LC_ALL=C.UTF-8
export LANG=C.UTF-8

# ---- EDIT ----
SPECIES_ID="1398205"
BAMROOT="/rds/project/rds-FSWe0O6MFwc/users/hb676/bowtie2_bam"
PRESORT_ROOT="/rds/project/rds-FSWe0O6MFwc/users/hb676/bowtie2_presorted"
OUTROOT="/rds/project/rds-FSWe0O6MFwc/users/hb676/bowtie2_sorted"
N_CHUNKS=4

module purge
module use /usr/local/software/spack/spack-modules/icelake-20220710/linux-rocky8-icelake
module load samtools/1.14/gcc/v46lwk2d

PRESORT_DIR="$PRESORT_ROOT/$SPECIES_ID"
OUTDIR="$OUTROOT/$SPECIES_ID"
mkdir -p "$OUTDIR"

CHUNK_PREFIX="$OUTDIR/${SPECIES_ID}.merged.sorted"

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

# ---- gather pre-sorted shards ----
mapfile -t SORTED_SHARDS < <(find "$PRESORT_DIR" -name "shard_*.sorted.sam.gz" | sort)
N="${#SORTED_SHARDS[@]}"
echo "Found $N pre-sorted shards."
if [[ "$N" -eq 0 ]]; then
    echo "ERROR: no presorted shards found in $PRESORT_DIR"
    exit 1
fi

# ---- get SAM header from first BAM of species ----
mapfile -t BAMS < <(find "$BAMROOT/$SPECIES_ID" -name "*.bam" | sort)
if [[ "${#BAMS[@]}" -eq 0 ]]; then
    echo "ERROR: no BAMs found in $BAMROOT/$SPECIES_ID"
    exit 1
fi

# ---- set up temp dir and cleanup ----
TMPDIR=$(mktemp -d "$PRESORT_DIR/tmp_merge_XXXXXX")
BGPIDS=()
FIFOS=()

cleanup() {
    kill "${BGPIDS[@]}" 2>/dev/null || true
    rm -rf "$TMPDIR"
}
trap cleanup EXIT

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

# ---- count total alignment lines across all shards (parallel) ----
echo "[$(date)] Counting alignment lines across $N shards..."
COUNT_PIDS=()
COUNT_FILES=()
for i in "${!SORTED_SHARDS[@]}"; do
    cf="$TMPDIR/cnt_${i}"
    zcat "${SORTED_SHARDS[$i]}" | wc -l > "$cf" &
    COUNT_PIDS+=($!)
    COUNT_FILES+=("$cf")
done
wait "${COUNT_PIDS[@]}"
TOTAL_LINES=0
for cf in "${COUNT_FILES[@]}"; do
    TOTAL_LINES=$(( TOTAL_LINES + $(cat "$cf") ))
done
echo "  Total alignment lines: $TOTAL_LINES"

# ---- create FIFOs to decompress shards in parallel ----
echo "[$(date)] Starting $N decompressor processes..."
for SHARD in "${SORTED_SHARDS[@]}"; do
    FIFO="$TMPDIR/$(basename "$SHARD" .sam.gz).fifo"
    mkfifo "$FIFO"
    zcat "$SHARD" > "$FIFO" &
    BGPIDS+=($!)
    FIFOS+=("$FIFO")
done

# ---- merge + split into N_CHUNKS ----
echo "[$(date)] Merging and splitting into $N_CHUNKS chunks -> ${CHUNK_PREFIX}.chunk*.sorted.sam.gz"
{
    samtools view -H "${BAMS[0]}"
    sort -m -k1,1 "${FIFOS[@]}"
} | python3 "$SPLIT_PY" \
    --n_chunks        "$N_CHUNKS" \
    --out_prefix      "$CHUNK_PREFIX" \
    --total_aln_lines "$TOTAL_LINES"

wait

# ---- validate chunks ----
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
echo "Pre-sorted shards in $PRESORT_DIR can now be deleted to free space."
