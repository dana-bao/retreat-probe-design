#!/usr/bin/env python3
'''
Counting statistics per species across the compressed filtering pipeline.

Columns:
  species_id           - species taxid
  total_tiles          - reads in original tiled FASTA
  ngslca_out           - reads with any bowtie2 alignment (summed across chunks)
  bamdam_out           - reads assigned at genus level or below (summed across chunks)
  bamdam_correct_genus - reads assigned to the correct genus (bamdam_filter output)
  kraken_total         - all reads processed by K-mer-based approach
  kraken_genus_level   - reads assigned to any genus level or below
  kraken_correct_genus - reads assigned to the correct genus

Prints "missing" for any file not available.

Variant of count_stats.py, adapted for the compressed pipeline where some outputs are split into chunks and/or merged. 
'''

import csv
import glob
import os
import re
import subprocess

TILED_READS  = "/rds/project/rds-FSWe0O6MFwc/users/hb676/tiled_reads"
NGSLCA_DIR   = "/rds/project/rds-FSWe0O6MFwc/users/hb676/ngslca_out_compressed"
BAMDAM_DIR   = "/rds/project/rds-FSWe0O6MFwc/users/hb676/bamdam_out_compressed"
BAMDAM_FILT  = "/rds/project/rds-FSWe0O6MFwc/users/hb676/bamdam_filtered_compressed"
KRAKEN_RAW    = "/rds/project/rds-FSWe0O6MFwc/users/hb676/kraken2"
KRAKEN_FILT   = "/rds/project/rds-FSWe0O6MFwc/users/hb676/kraken2_filtered"
EPROBE_INPUT  = "/rds/project/rds-FSWe0O6MFwc/users/hb676/eprobe_input"

OUT_CSV = "/rds/project/rds-FSWe0O6MFwc/users/hb676/pipeline_stats_compressed.csv"


def fasta_count(path):
    if not os.path.isfile(path):
        return None
    try:
        out = subprocess.check_output(['grep', '-c', '^>', path], universal_newlines=True)
        return int(out.strip())
    except subprocess.CalledProcessError:
        return 0


def line_count(path):
    if path is None or not os.path.isfile(path):
        return None
    with open(path) as f:
        return sum(1 for _ in f)


def sum_ngslca_chunks(ngslca_dir, species_id):
    '''Count assigned reads (taxid != 0) across all chunk .lca files.'''
    chunks = sorted(glob.glob(os.path.join(ngslca_dir, f'{species_id}_chunk*.lca')))
    if not chunks:
        return None
    total = 0
    for path in chunks:
        with open(path) as f:
            for line in f:
                if not line.strip():
                    continue
                parts = line.split('\t')
                if not parts[0] or parts[0].startswith('#'):
                    continue
                if len(parts) >= 2 and parts[1].strip().split(':')[0] != '0':
                    total += 1
    return total


def sum_bamdam_chunks(bamdam_dir, species_id):
    '''Sum line counts across all chunk .shrunk.lca files.'''
    pattern = os.path.join(bamdam_dir, str(species_id), f'{species_id}_chunk*.shrunk.lca')
    chunks = sorted(glob.glob(pattern))
    if not chunks:
        return None
    return sum(line_count(p) or 0 for p in chunks)


def sum_bamdam_filt_chunks(bamdam_filt_dir, species_id):
    '''Sum line counts across correct_genus filtered files (single or chunks).'''
    # Try single merged file first
    path = os.path.join(bamdam_filt_dir, f'{species_id}.shrunk.correct_genus.lca')
    if os.path.isfile(path):
        return line_count(path)
    # Fall back to chunk files
    chunks = sorted(glob.glob(os.path.join(bamdam_filt_dir, f'{species_id}_chunk*.shrunk.correct_genus.lca')))
    if not chunks:
        return None
    return sum(line_count(p) or 0 for p in chunks)


def sum_kraken(kraken_dir, species_id, suffix):
    pat = re.compile(
        rf"^{re.escape(str(species_id))}_task\d+\.k2\.0\.05\.core_nt\.{re.escape(suffix)}$"
    )
    total = 0
    found = False
    for fname in os.listdir(kraken_dir):
        if pat.match(fname):
            found = True
            with open(os.path.join(kraken_dir, fname)) as f:
                total += sum(1 for _ in f)
    return total if found else None


def fmt(v):
    return str(v) if v is not None else 'missing'


def main():
    species = sorted(
        f[:-6] for f in os.listdir(TILED_READS) if f.endswith('.fasta')
    )

    headers = ['species_id', 'total_tiles', 'ngslca_out', 'bamdam_out',
               'bamdam_correct_genus', 'kraken_total',
               'kraken_genus_level', 'kraken_correct_genus', 'overlap']

    with open(OUT_CSV, 'w', newline='') as fout:
        writer = csv.writer(fout)
        writer.writerow(headers)

        for sp in species:
            print(f"Processing {sp} ...", flush=True)
            tiles       = fasta_count(os.path.join(TILED_READS, f"{sp}.fasta"))
            ngslca      = sum_ngslca_chunks(NGSLCA_DIR, sp)
            bamdam      = sum_bamdam_chunks(BAMDAM_DIR, sp)
            bamdam_corr = sum_bamdam_filt_chunks(BAMDAM_FILT, sp)
            k_total     = sum_kraken(KRAKEN_RAW,  sp, 'out')
            k_genus     = sum_kraken(KRAKEN_FILT, sp, 'genus_level.out')
            k_correct   = sum_kraken(KRAKEN_FILT, sp, 'correct_genus.out')
            overlap     = fasta_count(os.path.join(EPROBE_INPUT, f"{sp}.fasta"))

            writer.writerow([sp, fmt(tiles), fmt(ngslca), fmt(bamdam),
                             fmt(bamdam_corr), fmt(k_total),
                             fmt(k_genus), fmt(k_correct), fmt(overlap)])

    print(f"Done. Output written to {OUT_CSV}")


if __name__ == '__main__':
    main()
