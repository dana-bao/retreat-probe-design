#!/usr/bin/env python3
'''
Counting statistics per species across the filtering pipeline.

Columns:
  species_id           - species taxid
  total_tiles          - reads in original tiled FASTA
  ngslca_out           - reads with any bowtie2 alignment (all ranks, including root-assigned;
                         one line per read in .lca)
  bamdam_out           - reads assigned at genus level or below (global alignment-based path, shrunk LCA)
  bamdam_correct_genus - reads assigned to the correct genus (global alignment-based path, bamdam_filter output)
  kraken_total         - all reads processed by K-mer-based approach (classified + unclassified)
  kraken_genus_level   - reads assigned to any genus level or below (K-mer-based path, relaxed)
  kraken_correct_genus - reads assigned to the correct genus (K-mer-based path, strict)
  overlap              - reads agreed upon by both K-mer-based and global alignment-based approach (eprobe_input)

Prints "missing" for any file not available.
'''

import csv
import os
import re
import subprocess

TILED_READS  = "/rds/project/rds-FSWe0O6MFwc/users/hb676/tiled_reads"
NGSLCA_DIR   = "/rds/project/rds-FSWe0O6MFwc/users/hb676/ngslca_out_new"
BAMDAM_DIR   = "/rds/project/rds-FSWe0O6MFwc/users/hb676/bamdam_out_new"
BAMDAM_FILT  = "/rds/project/rds-FSWe0O6MFwc/users/hb676/bamdam_filtered_new"
KRAKEN_RAW   = "/rds/project/rds-FSWe0O6MFwc/users/hb676/kraken2"
KRAKEN_FILT  = "/rds/project/rds-FSWe0O6MFwc/users/hb676/kraken2_filtered"
EPROBE_INPUT = "/rds/project/rds-FSWe0O6MFwc/users/hb676/eprobe_input_new"


def fasta_count(path):
    """Count FASTA records (lines starting with '>') using grep."""
    if not os.path.isfile(path):
        return None
    try:
        out = subprocess.check_output(['grep', '-c', '^>', path], universal_newlines=True)
        return int(out.strip())
    except subprocess.CalledProcessError:
        return 0


def line_count(path):
    """Count lines in a file."""
    if not os.path.isfile(path):
        return None
    with open(path) as f:
        return sum(1 for _ in f)


def sum_kraken(kraken_dir, species_id, suffix):
    """Sum lines across all task files for this species with the given suffix."""
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


OUT_CSV = "/rds/project/rds-FSWe0O6MFwc/users/hb676/pipeline_stats.csv"


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
            tiles        = fasta_count(os.path.join(TILED_READS,  f"{sp}.fasta"))
            ngslca       = line_count( os.path.join(NGSLCA_DIR,   f"{sp}.lca"))
            bamdam       = line_count( os.path.join(BAMDAM_DIR,   sp, f"{sp}.shrunk.lca"))
            bamdam_corr  = line_count( os.path.join(BAMDAM_FILT,  f"{sp}.shrunk.correct_genus.lca"))
            k_total      = sum_kraken(KRAKEN_RAW,  sp, 'out')
            k_genus      = sum_kraken(KRAKEN_FILT, sp, 'genus_level.out')
            k_correct    = sum_kraken(KRAKEN_FILT, sp, 'correct_genus.out')
            overlap      = fasta_count(os.path.join(EPROBE_INPUT, f"{sp}.fasta"))

            writer.writerow([sp, fmt(tiles), fmt(ngslca), fmt(bamdam),
                             fmt(bamdam_corr), fmt(k_total),
                             fmt(k_genus), fmt(k_correct), fmt(overlap)])

    print(f"Done. Output written to {OUT_CSV}")


if __name__ == '__main__':
    main()
