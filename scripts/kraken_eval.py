'''
Evaluate Kraken2 genus-level classification accuracy across confidence thresholds.
Read names are expected in the format: {true_species_taxid}|{contig}|{tile}|{pos}

For each read:
  - Unclassified (U)
  - Classified (C):
      - Walk up taxonomy from assigned taxid to find genus
      - Walk up taxonomy from true species taxid to find genus
      - above_genus: assigned taxid has no genus ancestor (classified too broadly)
      - correct_genus: assigned genus is the true genus
      - wrong_genus: assigned genus is not the true genus

Metrics:
  sensitivity = correct_genus / total_reads
  precision   = correct_genus / (correct_genus + wrong_genus)
  f1          = harmonic mean of sensitivity and precision
'''

import argparse
import csv
import os
import re
import sys
from collections import defaultdict

# import taxonomy information from kraken2 nodes.dmp file, which includes taxid, parent taxid, and rank for each node in the taxonomy tree.
def load_taxonomy(nodes_dmp):
    print(f"Loading taxonomy from {nodes_dmp} ...", flush=True)
    tax = {}
    with open(nodes_dmp) as f:
        for line in f:
            parts = line.split("\t|\t")
            taxid  = int(parts[0].strip())
            parent = int(parts[1].strip())
            rank   = parts[2].strip()
            tax[taxid] = (parent, rank)
    print(f"  Loaded {len(tax):,} nodes.", flush=True)
    return tax

# cache taxid to avoid redundant tree walks
_genus_cache = {}

# find genus-level taxid using current taxid, returns none if no genus found or above genus level
def get_genus_taxid(taxid, taxonomy):
    """Walk up from taxid to find genus-level ancestor. Returns genus taxid or None."""
    if taxid in _genus_cache:
        return _genus_cache[taxid]

    path = []
    current = taxid
    result = None

    ABOVE_GENUS = {"family", "order", "class", "phylum", "kingdom", "superkingdom"}

    while True:
        if current not in taxonomy:
            result = None
            break
        parent, rank = taxonomy[current]
        if rank == "genus":
            result = current
            break
        if rank in ABOVE_GENUS:
            result = None
            break
        if parent == current:
            result = None
            break
        path.append(current)
        current = parent

    for node in path:
        _genus_cache[node] = result
    _genus_cache[taxid] = result
    return result

# evaluate kraken2 output file using total reads, unclassified reads, and classification outcomes at genus level
def evaluate_file(out_file, taxonomy, true_species_taxid):
    true_genus = get_genus_taxid(true_species_taxid, taxonomy)
    if true_genus is None:
        print(f"  WARNING: no genus found for true species {true_species_taxid}", flush=True)

    counts = defaultdict(int)

    with open(out_file) as f:
        for line in f:
            counts["total"] += 1
            if line[0] == "U":
                counts["unclassified"] += 1
                continue

            # C, read_name, assigned_taxid, read_len, kmer_hits
            tab1 = line.index("\t")
            tab2 = line.index("\t", tab1 + 1)
            tab3 = line.index("\t", tab2 + 1)
            assigned_taxid = int(line[tab2 + 1:tab3])

            assigned_genus = get_genus_taxid(assigned_taxid, taxonomy)

            if assigned_genus is None:
                counts["above_genus"] += 1
            elif assigned_genus == true_genus:
                counts["correct_genus"] += 1
            else:
                counts["wrong_genus"] += 1

    return counts

# compute sensitivity, precision, and f1 score
def compute_metrics(counts):
    total      = counts["total"]
    correct    = counts["correct_genus"]
    wrong      = counts["wrong_genus"]
    above      = counts["above_genus"]
    unclassified = counts["unclassified"]
    at_genus   = correct + wrong

    sensitivity = correct / total    if total    > 0 else 0.0
    precision   = correct / at_genus if at_genus > 0 else 0.0
    f1 = (2 * sensitivity * precision / (sensitivity + precision)
          if (sensitivity + precision) > 0 else 0.0)

    return {
        "total_reads":       total,
        "unclassified":      unclassified,
        "above_genus":       above,
        "correct_genus":     correct,
        "wrong_genus":       wrong,
        "classified_at_genus": at_genus,
        "sensitivity":       round(sensitivity, 6),
        "precision":         round(precision,   6),
        "f1":                round(f1,           6),
    }

# extract confidence level and species taxid from kraken2 output file names 
def discover_files(results_dir):
    """Match {species_id}_task{N}.k2.{confidence}.core_nt.out"""
    pattern = re.compile(r"^(\d+)_task\d+\.k2\.([\d.]+)\.core_nt\.out$")
    files = []
    for fname in sorted(os.listdir(results_dir)):
        m = pattern.match(fname)
        if m:
            files.append((int(m.group(1)), m.group(2),
                          os.path.join(results_dir, fname)))
    return files

# main body of running
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--results_dir", required=True,
                        help="Directory containing Kraken2 .out files")
    parser.add_argument("--nodes", required=True,
                        help="Path to nodes.dmp from Kraken2 DB taxonomy")
    parser.add_argument("--out", default="kraken_eval.csv",
                        help="Output CSV (default: kraken_eval.csv)")
    args = parser.parse_args()

    taxonomy = load_taxonomy(args.nodes)

    files = discover_files(args.results_dir)
    if not files:
        print(f"ERROR: no matching .out files in {args.results_dir}", file=sys.stderr)
        sys.exit(1)

    print(f"Found {len(files)} files to process.", flush=True)

    rows = []
    for species_id, confidence, filepath in files:
        print(f"  [conf={confidence}] species {species_id} ...", flush=True)
        counts  = evaluate_file(filepath, taxonomy, species_id)
        metrics = compute_metrics(counts)
        rows.append({"species_id": species_id, "confidence": confidence, **metrics})

    fieldnames = ["species_id", "confidence", "total_reads", "unclassified",
                  "above_genus", "correct_genus", "wrong_genus",
                  "classified_at_genus", "sensitivity", "precision", "f1"]

    with open(args.out, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"\nResults written to: {args.out}", flush=True)

    # Quick summary table
    print(f"\n{'species_id':>12} {'conf':>6} {'sensitivity':>12} {'precision':>10} {'f1':>8}")
    print("-" * 54)
    for r in sorted(rows, key=lambda x: (x["species_id"], x["confidence"])):
        print(f"{r['species_id']:>12} {r['confidence']:>6} "
              f"{r['sensitivity']:>12.4f} {r['precision']:>10.4f} {r['f1']:>8.4f}")


if __name__ == "__main__":
    main()
