'''
Filter bamdam .shrunk.lca files to correct genus.

Read names expected in format: {true_species_taxid}|{contig}|{tile}|{pos}
True species taxid is acquired from the filename:
  {species_id}.shrunk.lca  or  {species_id}_modern.shrunk.lca

.lca file format (tab-separated):
  Column 1: read name
  Column 2: assigned taxid (0 = unclassified)
  Column 3+: additional fields (ignored)

Output correct_genus.lca: reads where assigned genus matches true genus
'''

import argparse
import os
import re
import sys

# import taxonomy information from nodes.dmp file, which includes taxid, parent taxid, and rank for each node in the taxonomy tree.
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


_genus_cache = {}

# find genus-level taxid using current taxid, returns none if no genus found or above genus level
def get_genus_taxid(taxid, taxonomy):
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

# filter based on the criteria described on top, writing correct genus reads only
def filter_file(in_file, out_correct, taxonomy, true_species_taxid):
    true_genus = get_genus_taxid(true_species_taxid, taxonomy)
    if true_genus is None:
        print(f"  WARNING: no genus found for true species {true_species_taxid}", flush=True)

    counts = {"total": 0, "unclassified": 0, "above_genus": 0,
              "correct_genus": 0, "wrong_genus": 0}

    with open(in_file) as fin, open(out_correct, "w") as fout_correct:
        for line in fin:
            if not line.strip():
                continue
            counts["total"] += 1

            parts = line.split("\t")
            assigned_taxid = int(parts[1].split(':')[0].strip())

            if assigned_taxid == 0:
                counts["unclassified"] += 1
                continue

            assigned_genus = get_genus_taxid(assigned_taxid, taxonomy)

            if assigned_genus is None:
                counts["above_genus"] += 1
            elif assigned_genus == true_genus:
                counts["correct_genus"] += 1
                fout_correct.write(line)
            else:
                counts["wrong_genus"] += 1

    return counts

# find chunk .shrunk.lca files 
def discover_files(lca_dir):
    pattern = re.compile(r"^(\d+)(_modern)?\.shrunk\.lca$")
    files = []
    for root, dirs, fnames in os.walk(lca_dir):
        for fname in sorted(fnames):
            m = pattern.match(fname)
            if m:
                files.append((int(m.group(1)), fname, os.path.join(root, fname)))
    return sorted(files, key=lambda x: x[1])

# main functions
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--lca_dir", required=True,
                        help="Directory containing bamdam .shrunk.lca files (searched recursively)")
    parser.add_argument("--nodes", required=True,
                        help="Path to nodes.dmp from taxonomy")
    parser.add_argument("--out_dir", required=True,
                        help="Directory to write filtered output files")
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    taxonomy = load_taxonomy(args.nodes)

    files = discover_files(args.lca_dir)
    if not files:
        print(f"ERROR: no .shrunk.lca files found under {args.lca_dir}", file=sys.stderr)
        sys.exit(1)

    print(f"Found {len(files)} files to process.", flush=True)

    for species_id, fname, filepath in files:
        base = fname.replace(".shrunk.lca", "")
        out_correct_path = os.path.join(args.out_dir, f"{base}.shrunk.correct_genus.lca")

        if os.path.isfile(out_correct_path):
            print(f"  species {species_id} ({fname}) -> already done, skipping.", flush=True)
            continue

        print(f"  species {species_id} ({fname}) -> {base}.shrunk.correct_genus.lca", flush=True)
        counts = filter_file(filepath, out_correct_path, taxonomy, species_id)

        total   = counts["total"]
        correct = counts["correct_genus"]
        pct_correct = 100 * correct / total if total > 0 else 0
        print(f"    correct_genus: {correct:,} / {total:,} reads ({pct_correct:.3f}%)", flush=True)


if __name__ == "__main__":
    main()
