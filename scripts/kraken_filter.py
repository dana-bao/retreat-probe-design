'''
Filter Kraken2 .out file to keep reads classified at genus level.

Produces two output files per input:
  1. correct_genus.out, which contain reads assigned to the correct genus (strict filter):
       1. Classified
       2. Assigned taxid is at genus level or below
       3. Assigned genus matches the true genus encoded in the read name
  2. genus_level.out, which contain reads assigned to any genus (relaxed filter, analogous to bamdam --upto genus):
       1. Classified
       2. Assigned taxid is at genus level or below
       (genus does not need to match true genus)

Read names expected in format: {true_species_taxid}|{contig}|{tile}|{pos}
True species taxid is acquired from the filename: {species_id}_task{N}.k2.{conf}.core_nt.out
'''

import argparse
import os
import re
import sys

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

# filter based on the criteria described on top, writing two output files
def filter_file(in_file, out_correct, out_genus, taxonomy, true_species_taxid):
    true_genus = get_genus_taxid(true_species_taxid, taxonomy)
    if true_genus is None:
        print(f"  WARNING: no genus found for true species {true_species_taxid}", flush=True)

    counts = {"total": 0, "unclassified": 0, "above_genus": 0,
              "correct_genus": 0, "wrong_genus": 0}

    with open(in_file) as fin, open(out_correct, "w") as fout_correct, open(out_genus, "w") as fout_genus:
        for line in fin:
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
                fout_correct.write(line)
                fout_genus.write(line)
            else:
                counts["wrong_genus"] += 1
                fout_genus.write(line)

    return counts

# only process files with confidence level 0.2 and extract species taxid
def discover_files(results_dir, confidence="0.2"):
    pattern = re.compile(r"^(\d+)(_modern)?_task\d+\.k2\.([\d.]+)\.core_nt\.out$")
    files = []
    for fname in sorted(os.listdir(results_dir)):
        m = pattern.match(fname)
        if m and m.group(3) == confidence:
            files.append((int(m.group(1)), m.group(3), fname,
                          os.path.join(results_dir, fname)))
    return files

# main body of running
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--results_dir", required=True,
                        help="Directory containing Kraken2 .out files")
    parser.add_argument("--nodes", required=True,
                        help="Path to nodes.dmp from Kraken2 DB taxonomy")
    parser.add_argument("--out_dir", required=True,
                        help="Directory to write filtered output files")
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    taxonomy = load_taxonomy(args.nodes)

    files = discover_files(args.results_dir)
    if not files:
        print(f"ERROR: no matching .out files in {args.results_dir}", file=sys.stderr)
        sys.exit(1)

    print(f"Found {len(files)} files to process.", flush=True)

    for species_id, confidence, fname, filepath in files:
        out_correct_fname = fname.replace(".core_nt.out", ".core_nt.correct_genus.out")
        out_correct_path  = os.path.join(args.out_dir, out_correct_fname)
        out_genus_fname   = fname.replace(".core_nt.out", ".core_nt.genus_level.out")
        out_genus_path    = os.path.join(args.out_dir, out_genus_fname)

        if os.path.isfile(out_correct_path) and os.path.isfile(out_genus_path):
            print(f"  [conf={confidence}] species {species_id} -> already done, skipping.", flush=True)
            continue

        print(f"  [conf={confidence}] species {species_id} -> {out_correct_fname}, {out_genus_fname}", flush=True)
        counts = filter_file(filepath, out_correct_path, out_genus_path, taxonomy, species_id)

        total   = counts["total"]
        correct = counts["correct_genus"]
        genus   = counts["correct_genus"] + counts["wrong_genus"]
        pct_correct = 100 * correct / total if total > 0 else 0
        pct_genus   = 100 * genus   / total if total > 0 else 0
        print(f"    correct_genus: {correct:,} / {total:,} reads ({pct_correct:.3f}%)", flush=True)
        print(f"    genus_level:   {genus:,} / {total:,} reads ({pct_genus:.3f}%)", flush=True)


if __name__ == "__main__":
    main()
