#!/usr/bin/env bash
set -euo pipefail

# Optional: make sure your personal bin (where seqkit lives) is available
export PATH="$HOME/bin:$PATH"

# Safety: if no files match the glob, expand to nothing instead of literal text
shopt -s nullglob

for f in ncbi_dataset/data/*/*_genomic.fna; do
  total=$(seqkit stats -T "$f" | awk 'NR==2{print $5}')
  asm=$(basename "$(dirname "$f")")
  echo -e "$total\t$asm"
done | sort -nr | head -10
