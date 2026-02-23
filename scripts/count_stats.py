#!/usr/bin/env python3
'''
Counting statistics per species across the filtering pipeline.

Columns:
  species_id          - species taxid
  total_tiles         - reads in original tiled FASTA
  ngslca_processed    - reads with any bowtie2 alignment (all ranks, one line per read in .lca)
  bamdam_shrunk       - reads assigned at genus level or below (bowtie2 path)
  kraken_correct_genus - reads assigned to the correct genus (Kraken2 path, strict)
  kraken_genus_level  - reads assigned to any genus level or below (Kraken2 path, relaxed)

Prints "missing" for any file not yet available.
'''

import os
import re
import subprocess

TILED_READS  = "/rds/project/rds-FSWe0O6MFwc/users/hb676/tiled_reads"
NGSLCA_DIR   = "/rds/project/rds-FSWe0O6MFwc/users/hb676/ngslca_out"
BAMDAM_FASTA = "/rds/project/rds-FSWe0O6MFwc/users/hb676/bamdam_fasta"
KRAKEN_FILT  = "/rds/project/rds-FSWe0O6MFwc/users/hb676/kraken2_filtered"


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


def sum_kraken(filtered_dir, species_id, suffix):
    """Sum lines across all task files for this species with the given suffix."""
    pat = re.compile(
        rf"^{re.escape(str(species_id))}_task\d+\.k2\.[\d.]+\.core_nt\.{re.escape(suffix)}$"
    )
    total = 0
    found = False
    for fname in os.listdir(filtered_dir):
        if pat.match(fname):
            found = True
            with open(os.path.join(filtered_dir, fname)) as f:
                total += sum(1 for _ in f)
    return total if found else None


def fmt(v):
    return str(v) if v is not None else 'missing'


def main():
    species = sorted(
        f[:-6] for f in os.listdir(TILED_READS) if f.endswith('.fasta')
    )

    headers = ['species_id', 'total_tiles', 'ngslca_processed', 'bamdam_shrunk',
               'kraken_correct_genus', 'kraken_genus_level']
    print('\t'.join(headers))

    for sp in species:
        tiles     = fasta_count(os.path.join(TILED_READS,  f"{sp}.fasta"))
        ngslca    = line_count( os.path.join(NGSLCA_DIR,   f"{sp}.lca"))
        bamdam    = fasta_count(os.path.join(BAMDAM_FASTA, f"{sp}.shrunk.fasta"))
        k_correct = sum_kraken(KRAKEN_FILT, sp, 'correct_genus.out')
        k_genus   = sum_kraken(KRAKEN_FILT, sp, 'genus_level.out')

        print('\t'.join([sp, fmt(tiles), fmt(ngslca), fmt(bamdam),
                         fmt(k_correct), fmt(k_genus)]))


if __name__ == '__main__':
    main()
