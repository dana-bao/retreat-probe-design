#!/usr/bin/env bash
#SBATCH --job-name=bamsort_comp_mod
#SBATCH --partition=icelake
#SBATCH --cpus-per-task=20
#SBATCH --mem=240G
#SBATCH --time=36:00:00
#SBATCH --array=0-6%7
#SBATCH --output=/rds/project/rds-FSWe0O6MFwc/users/zh384/logs/%x_%A_%a.out
#SBATCH --error=/rds/project/rds-FSWe0O6MFwc/users/zh384/logs/%x_%A_%a.err

# Merges *.compressed.bam shards for 7 modern species, sorts by read name,
# and splits into N_CHUNKS queryname-sorted BAM chunks for downstream ngsLCA/bamdam.
# Output: {OUTROOT}/{SPECIES_ID}/{SPECIES_ID}.chunk{00-03}.sorted.bam
# Skips species where all chunks already exist and pass quickcheck.
# Originals are never modified.

set -euo pipefail
export LC_ALL=C.UTF-8
export LANG=C.UTF-8

BAMROOT="/rds/project/rds-FSWe0O6MFwc/users/hb676/modern_bowtie2_bam"
OUTROOT="/rds/project/rds-FSWe0O6MFwc/users/hb676/modern_bowtie2_sorted_compressed"
N_CHUNKS=4

SPECIES_LIST=(
    111071_modern
    1250597_modern
    344590_modern
    357475_modern
    39862_modern
    408822_modern
    42846_modern
)

THREADS="${SLURM_CPUS_PER_TASK:-20}"
THREADS="${THREADS//[^0-9]/}"
[[ -n "$THREADS" ]] || { echo "ERROR: THREADS empty"; exit 10; }
SORT_MEM_PER_THREAD="10G"   # 20 threads × 10G = 200G

mkdir -p "$OUTROOT" /rds/project/rds-FSWe0O6MFwc/users/zh384/logs

module purge
module use /usr/local/software/spack/spack-modules/icelake-20220710/linux-rocky8-icelake
module load samtools/1.14/gcc/v46lwk2d

echo "HOST=$(hostname)"
samtools --version | head -n 2
echo "THREADS=$THREADS  SORT_MEM_PER_THREAD=$SORT_MEM_PER_THREAD"

N="${#SPECIES_LIST[@]}"
TASK="${SLURM_ARRAY_TASK_ID:-0}"
if [[ "$TASK" -ge "$N" ]]; then
    echo "TASK=$TASK out of range (N=$N) -> exiting cleanly."
    exit 0
fi

SPECIES_ID="${SPECIES_LIST[$TASK]}"
SPECIES_PATH="${BAMROOT}/${SPECIES_ID}"
echo "[$(date)] TASK=$TASK SPECIES_ID=$SPECIES_ID"

# ---- gather compressed BAMs ----
shopt -s nullglob
BAMS=( "$SPECIES_PATH"/*.compressed.bam )
shopt -u nullglob

if [[ "${#BAMS[@]}" -eq 0 ]]; then
    echo "No compressed BAMs found in $SPECIES_PATH; skipping."
    exit 0
fi
echo "Input compressed BAMs: ${#BAMS[@]}"

OUTDIR="$OUTROOT/$SPECIES_ID"
mkdir -p "$OUTDIR"
CHUNK_PREFIX="$OUTDIR/${SPECIES_ID}"

# ---- skip if all chunks already done ----
all_done=true
for i in $(seq 0 $((N_CHUNKS - 1))); do
    cf="${CHUNK_PREFIX}.chunk$(printf '%02d' "$i").sorted.bam"
    if ! [[ -s "$cf" ]] || ! samtools quickcheck "$cf" 2>/dev/null; then
        all_done=false; break
    fi
done
if $all_done; then
    echo "All $N_CHUNKS chunks already exist and pass quickcheck. Skipping."
    exit 0
fi

TMPDIR_WORK="$OUTDIR/tmp_${SLURM_JOB_ID:-local}_${SLURM_ARRAY_TASK_ID:-0}"
mkdir -p "$TMPDIR_WORK"
trap 'rm -rf "$TMPDIR_WORK"' EXIT

MERGED_TMP="$TMPDIR_WORK/${SPECIES_ID}.merged.tmp.bam"
SORTED_TMP="$TMPDIR_WORK/${SPECIES_ID}.sorted.tmp.bam"

# ---- merge all compressed shards ----
echo "[$(date)] Merging ${#BAMS[@]} compressed BAMs -> $MERGED_TMP"
samtools merge -@ "$THREADS" -f -n "$MERGED_TMP" "${BAMS[@]}"

# ---- sort by read name ----
echo "[$(date)] Sorting by name -> $SORTED_TMP"
samtools sort -@ "$THREADS" -m "$SORT_MEM_PER_THREAD" -n \
    -T "$TMPDIR_WORK/sorttmp" -O bam -o "$SORTED_TMP" "$MERGED_TMP"
rm -f "$MERGED_TMP"

# ---- count reads for splitter ----
echo "[$(date)] Counting reads..."
TOTAL=$(samtools view -c -@ "$THREADS" "$SORTED_TMP")
echo "  Total reads: $TOTAL"

# ---- write Python BAM splitter ----
SPLIT_PY="$TMPDIR_WORK/bam_split.py"
cat > "$SPLIT_PY" << 'PYEOF'
"""Split a SAM stream (samtools view -h) into N queryname-sorted BAM chunks."""
import argparse, subprocess, sys

ap = argparse.ArgumentParser()
ap.add_argument("--n_chunks",        type=int, required=True)
ap.add_argument("--out_prefix",               required=True)
ap.add_argument("--total_aln_lines", type=int, required=True)
args = ap.parse_args()

n      = args.n_chunks
target = max(1, args.total_aln_lines // n)

header_lines, first_aln = [], None
for raw in sys.stdin.buffer:
    if raw.startswith(b"@"):
        header_lines.append(raw)
    else:
        first_aln = raw
        break

header_bytes = b"".join(header_lines)

def open_chunk(i):
    path = f"{args.out_prefix}.chunk{i:02d}.sorted.bam"
    proc = subprocess.Popen(["samtools", "view", "-b", "-o", path, "-"],
                            stdin=subprocess.PIPE)
    proc.stdin.write(header_bytes)
    return proc, path

if first_aln is None:
    print("[bam_split] WARNING: no alignment lines; writing empty chunks", file=sys.stderr)
    for i in range(n):
        proc, path = open_chunk(i)
        proc.stdin.close(); proc.wait()
        print(f"[bam_split] chunk {i:02d} -> {path} (empty)", file=sys.stderr, flush=True)
    sys.exit(0)

chunk_idx = 0
proc, path = open_chunk(0)
print(f"[bam_split] chunk {chunk_idx:02d} -> {path}", file=sys.stderr, flush=True)

lines_in_chunk = 0
prev_name      = None
line           = first_aln

while line:
    name = line.split(b"\t", 1)[0]
    if chunk_idx < n - 1 and lines_in_chunk >= target and prev_name and name != prev_name:
        proc.stdin.close(); proc.wait()
        chunk_idx += 1
        lines_in_chunk = 0
        proc, path = open_chunk(chunk_idx)
        print(f"[bam_split] chunk {chunk_idx:02d} -> {path}", file=sys.stderr, flush=True)
    proc.stdin.write(line)
    lines_in_chunk += 1
    prev_name = name
    line = sys.stdin.buffer.readline()

proc.stdin.close(); proc.wait()

# fill any remaining chunks (e.g. species with very few reads)
for i in range(chunk_idx + 1, n):
    proc, path = open_chunk(i)
    proc.stdin.close(); proc.wait()
    print(f"[bam_split] chunk {i:02d} -> {path} (empty)", file=sys.stderr, flush=True)

print(f"[bam_split] done — {n} chunks", file=sys.stderr, flush=True)
PYEOF

# ---- split into N_CHUNKS BAM files ----
echo "[$(date)] Splitting into $N_CHUNKS chunks -> ${CHUNK_PREFIX}.chunk*.sorted.bam"
samtools view -h -@ "$THREADS" "$SORTED_TMP" | python3 "$SPLIT_PY" \
    --n_chunks        "$N_CHUNKS" \
    --out_prefix      "$CHUNK_PREFIX" \
    --total_aln_lines "$TOTAL"
rm -f "$SORTED_TMP"

# ---- validate ----
echo "[$(date)] Validating chunks..."
for i in $(seq 0 $((N_CHUNKS - 1))); do
    cf="${CHUNK_PREFIX}.chunk$(printf '%02d' "$i").sorted.bam"
    if [[ -s "$cf" ]]; then
        samtools quickcheck "$cf" && echo "  chunk $i ok: $cf"
    else
        echo "  WARNING: chunk $i missing or empty: $cf"
    fi
done

echo "[$(date)] Done: $SPECIES_ID"
