"""
Extract eProbe input FASTA for species where overlap is too low.
Uses bamdam_filtered read names only (no Kraken2 intersection),
then fishes the corresponding undamaged sequences from modern_tiled.

Output: eprobe_input/{sid}_bamdam.fasta

Usage:
    python3 eprobe_input_bamdam_only.py \
        --species_id      56490 \
        --bamdam_filt_dir /rds/.../bamdam_filtered_compressed \
        --modern_tiled_dir /rds/.../modern_tiled \
        --out_dir         /rds/.../eprobe_input

variant of eprobe_input_extract.py, for annotation see eprobe_input_extract.py
"""

import argparse
import os
import sys


def find_bamdam_filtered(bamdam_filt_dir, species_id):
    path = os.path.join(bamdam_filt_dir, f'{species_id}.shrunk.correct_genus.lca')
    return path if os.path.isfile(path) else None


def find_modern_tiled(modern_tiled_dir, species_id):
    for name in (f'{species_id}.fasta', f'{species_id}_modern.fasta'):
        path = os.path.join(modern_tiled_dir, name)
        if os.path.isfile(path):
            return path
    return None


def read_names_bamdam(path):
    """Column 0 of the LCA file, stripped at the first ':'."""
    names = set()
    with open(path) as f:
        for line in f:
            if line.strip():
                parts = line.split('\t')
                if parts:
                    names.add(parts[0].split(':')[0])
    return names


def write_subset_from_fasta(fasta_path, wanted, out_fh):
    written = 0
    write_this = False
    seen = set()
    with open(fasta_path) as f:
        for line in f:
            if line.startswith('>'):
                name = line[1:].split()[0]
                write_this = name in wanted
                if write_this:
                    out_fh.write(line)
                    seen.add(name)
                    written += 1
            elif write_this:
                out_fh.write(line)
    not_found = len(wanted - seen)
    return written, not_found


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--species_id',        required=True)
    ap.add_argument('--bamdam_filt_dir',   required=True)
    ap.add_argument('--modern_tiled_dir',  required=True)
    ap.add_argument('--out_dir',           required=True)
    args = ap.parse_args()

    sid = args.species_id

    bamdam_path = find_bamdam_filtered(args.bamdam_filt_dir, sid)
    tiled_path  = find_modern_tiled(args.modern_tiled_dir, sid)

    missing = []
    if bamdam_path is None:
        missing.append('bamdam_filtered')
    if tiled_path is None:
        missing.append('modern_tiled')
    if missing:
        print(f"SKIP {sid}: missing {', '.join(missing)}")
        sys.exit(1)

    os.makedirs(args.out_dir, exist_ok=True)
    out_path = os.path.join(args.out_dir, f'{sid}_bamdam.fasta')

    if os.path.isfile(out_path) and os.path.getsize(out_path) > 0:
        print(f"Already done: {out_path} — skipping.")
        sys.exit(0)

    b_names = read_names_bamdam(bamdam_path)
    print(f"{sid}: bamdam={len(b_names):,} read names", flush=True)

    with open(out_path, 'w') as fout:
        written, not_found = write_subset_from_fasta(tiled_path, b_names, fout)

    if not_found:
        print(f"  WARNING: {not_found:,} bamdam names not found in modern_tiled FASTA")

    print(f"  -> {written:,} sequences written to {out_path}")


if __name__ == '__main__':
    main()
