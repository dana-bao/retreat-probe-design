"""
Extract eProbe input FASTAs by:
  1. Computing the overlap of Kraken2-correct and bamdam-correct read names
  2. Fishing the corresponding undamaged sequences from modern_tiled FASTAs

For each species, the three required inputs are:
  - kraken2_filtered/{sid}_task*.k2.0.2.core_nt.correct_genus.out
  - bamdam_filtered/{sid}.shrunk.correct_genus.lca
  - modern_tiled/{sid}.fasta

Species missing any of these required files are reported and skipped.
Output: eprobe_input/{sid}.fasta

Usage:
    python3 scripts/eprobe_input_extract.py \\
        --kraken_filt_dir /rds/.../kraken2_filtered \\
        --bamdam_filt_dir /rds/.../bamdam_filtered \\
        --modern_tiled_dir /rds/.../modern_tiled \\
        --out_dir /rds/.../eprobe_input
"""

import argparse
import glob
import os
import sys

# file finders  (same logic as pipeline_visualize.py)
def find_kraken_filtered(kraken_filt_dir, species_id):
    matches = glob.glob(
        os.path.join(kraken_filt_dir,
                     f'{species_id}_task*.k2.0.2.core_nt.correct_genus.out')
    )
    return matches[0] if matches else None

def find_bamdam_filtered(bamdam_filt_dir, species_id):
    path = os.path.join(bamdam_filt_dir, f'{species_id}.shrunk.correct_genus.lca')
    return path if os.path.isfile(path) else None

def find_modern_tiled(modern_tiled_dir, species_id):
    for name in (f'{species_id}.fasta', f'{species_id}_modern.fasta'):
        path = os.path.join(modern_tiled_dir, name)
        if os.path.isfile(path):
            return path
    return None

# read-name extractors for kraken, which is index 1 in kraken out files
def read_names_kraken(path):
    names = set()
    with open(path) as f:
        for line in f:
            if line.strip():
                parts = line.split('\t')
                if len(parts) >= 2:
                    names.add(parts[1].strip())
    return names

# read-name extractor for bamdam, which is index 0 in the LCA file, stripped at the first ':'
def read_names_bamdam(path):
    names = set()
    with open(path) as f:
        for line in f:
            if line.strip():
                parts = line.split('\t')
                if parts:
                    names.add(parts[0].split(':')[0])
    return names

# streaming FASTA writer
def write_overlap_from_fasta(fasta_path, wanted, out_fh):
    """
    Stream through fasta_path, writing records whose name is in `wanted`.
    Returns (written_count, not_found_count).
    """
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


# main processing function for each species
def process_species(sid, args):
    kraken_path = find_kraken_filtered(args.kraken_filt_dir, sid)
    bamdam_path = find_bamdam_filtered(args.bamdam_filt_dir, sid)
    tiled_path  = find_modern_tiled(args.modern_tiled_dir, sid)

    missing = []
    if kraken_path is None:
        missing.append('kraken2_filtered')
    if bamdam_path is None:
        missing.append('bamdam_filtered')
    if tiled_path is None:
        missing.append('modern_tiled')

    if missing:
        print(f"SKIP {sid}: missing {', '.join(missing)}")
        sys.exit(0)

    out_path = os.path.join(args.out_dir, f'{sid}.fasta')
    if os.path.isfile(out_path) and os.path.getsize(out_path) > 0:
        print(f"Already done: {out_path} — skipping.")
        sys.exit(0)

    k_names = read_names_kraken(kraken_path)
    b_names = read_names_bamdam(bamdam_path)
    overlap = k_names & b_names

    print(f"{sid}: kraken={len(k_names):,}  bamdam={len(b_names):,}  "
          f"overlap={len(overlap):,}", flush=True)

    with open(out_path, 'w') as fout:
        written, not_found = write_overlap_from_fasta(tiled_path, overlap, fout)

    if not_found:
        print(f"  WARNING: {not_found:,} overlap names not found in modern_tiled FASTA")

    print(f"  -> {written:,} sequences written to {out_path}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--species_id',       required=True,
                    help='Species taxid to process')
    ap.add_argument('--kraken_filt_dir',  required=True,
                    help='kraken2_filtered directory')
    ap.add_argument('--bamdam_filt_dir',  required=True,
                    help='bamdam_filtered directory (aeDNA)')
    ap.add_argument('--modern_tiled_dir', required=True,
                    help='modern_tiled directory')
    ap.add_argument('--out_dir',          required=True,
                    help='Output directory for eprobe_input FASTAs')
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    process_species(args.species_id, args)


if __name__ == '__main__':
    main()
