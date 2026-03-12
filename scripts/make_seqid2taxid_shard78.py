"""
Generate seqid2taxid entries for the 82 new assemblies (shard 78, 79, 80).

Reads assembly_data_report.jsonl to get GCA accession → taxid, then scans
each *_genomic.fna for sequence headers.  Produces two output files ready
to append to the existing kraken2_taxonomy files:

  shard78_seqid2taxid.map     → append to seqid2taxid.map
  shard78_seqid2taxid.acc2tax → append to seqid2taxid.acc2tax
                                 (same format as the awk-generated file:
                                  seqid \\t seqid \\t taxid \\t 0)

Usage:
    python3 scripts/make_seqid2taxid_shard78.py \\
        --assembly_dir /home/hb676/rds/hpc-work/probe_design/data/assemblies/ncbi_dataset/data \\
        --out_dir      /home/hb676/rds/hpc-work/probe_design/data/assemblies

Then on the cluster append:
    cat shard78_seqid2taxid.map     >> kraken2_taxonomy/seqid2taxid.map
    cat shard78_seqid2taxid.acc2tax >> kraken2_taxonomy/seqid2taxid.acc2tax
"""

import argparse
import gzip
import json
import os
import sys


def load_taxids(jsonl_path):
    """
    Parse assembly_data_report.jsonl.
    Returns dict: GCA_accession (no version) -> taxid str.
    e.g. 'GCA_000239015' -> '12345'
    """
    taxid_map = {}
    n_records = 0
    with open(jsonl_path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            rec = json.loads(line)
            acc = rec.get('accession', '')           # e.g. 'GCA_000239015.2'
            taxid = str(rec.get('organism', {}).get('taxId', ''))
            if acc and taxid:
                taxid_map[acc] = taxid
                taxid_map[acc.rsplit('.', 1)[0]] = taxid   # bare, no version
                n_records += 1
    print(f"[jsonl] {n_records} records in {os.path.basename(jsonl_path)} "
          f"(may include assemblies not present as directories)", file=sys.stderr)
    return taxid_map


def iter_seqids(fna_path):
    """
    Yield sequence IDs (text before first space on '>' lines) from a FASTA.
    Handles plain or .gz files.
    """
    opener = gzip.open if fna_path.endswith('.gz') else open
    with opener(fna_path, 'rt') as f:
        for line in f:
            if line.startswith('>'):
                seqid = line[1:].split()[0]
                yield seqid


def main():
    ap = argparse.ArgumentParser(
        description='Generate seqid2taxid entries for 82 new assemblies'
    )
    ap.add_argument('--assembly_dir', required=True,
                    help='Path to ncbi_dataset/data/ containing GCA_* subdirs and assembly_data_report.jsonl')
    ap.add_argument('--out_dir', default='.',
                    help='Directory for output files (default: current dir)')
    args = ap.parse_args()

    assembly_dir = args.assembly_dir
    jsonl = os.path.join(assembly_dir, 'assembly_data_report.jsonl')

    if not os.path.isfile(jsonl):
        print(f"ERROR: {jsonl} not found", file=sys.stderr)
        sys.exit(1)

    taxid_map = load_taxids(jsonl)

    os.makedirs(args.out_dir, exist_ok=True)
    out_map     = os.path.join(args.out_dir, 'shard78_seqid2taxid.map')
    out_acc2tax = os.path.join(args.out_dir, 'shard78_seqid2taxid.acc2tax')

    n_seqs   = 0
    n_done   = 0
    n_miss   = 0
    missing_gcas = []

    with open(out_map, 'w') as fm, open(out_acc2tax, 'w') as fa:
        for entry in sorted(os.scandir(assembly_dir), key=lambda e: e.name):
            if not entry.is_dir():
                continue
            gca = entry.name          # e.g. GCA_000239015.2 or GCF_000239015.2
            if not (gca.startswith('GCA_') or gca.startswith('GCF_')):
                continue

            # find taxid
            taxid = taxid_map.get(gca) or taxid_map.get(gca.rsplit('.', 1)[0])
            if not taxid:
                print(f"  WARNING: no taxid for {gca} — skipping", file=sys.stderr)
                missing_gcas.append(gca)
                n_miss += 1
                continue

            # find .fna file
            fnas = [f for f in os.listdir(entry.path) if f.endswith('_genomic.fna')]
            if not fnas:
                print(f"  WARNING: no *_genomic.fna in {entry.path} — skipping", file=sys.stderr)
                missing_gcas.append(gca)
                n_miss += 1
                continue

            fna_path = os.path.join(entry.path, fnas[0])
            for seqid in iter_seqids(fna_path):
                fm.write(f"{seqid}\t{taxid}\n")
                fa.write(f"{seqid}\t{seqid}\t{taxid}\t0\n")
                n_seqs += 1
            n_done += 1

    print(f"\nWrote {n_seqs} sequence entries across {n_done} assemblies.",
          file=sys.stderr)
    if missing_gcas:
        print(f"Skipped {n_miss} assemblies (no taxid or no .fna): {missing_gcas}", file=sys.stderr)

    print(f"\nOutput files:")
    print(f"  {out_map}")
    print(f"  {out_acc2tax}")
    print(f"\nTo update the taxonomy, run on the cluster:")
    print(f"  cat {out_acc2tax} >> /rds/project/rds-FSWe0O6MFwc/users/hb676/kraken2_taxonomy/seqid2taxid.acc2tax")
    print(f"  cat {out_map}     >> /rds/project/rds-FSWe0O6MFwc/users/hb676/kraken2_taxonomy/seqid2taxid.map")


if __name__ == '__main__':
    main()
