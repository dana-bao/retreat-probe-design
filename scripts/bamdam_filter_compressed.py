'''
Filter bamdam chunk .shrunk.lca files to reads belonging to the correct genus

Discovers {species_id}[_modern]_chunk*.shrunk.lca files under lca_dir,
groups them by species, reads all chunks in memory sequentially, and
writes filtered output in full. 

Read names expected in format: {true_species_taxid}|{contig}|{tile}|{pos}

Produces per species {base}.shrunk.correct_genus.lca

Variant of bamdam_filter.py
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

# filter based on the criteria described on top, writing only one output file with correct genus reads
def filter_chunks(chunk_paths, out_correct, taxonomy, true_species_taxid):
    '''Read multiple chunk .shrunk.lca files sequentially, write filtered output.'''
    true_genus = get_genus_taxid(true_species_taxid, taxonomy)
    if true_genus is None:
        print(f"  WARNING: no genus found for true species {true_species_taxid}", flush=True)

    counts = {"total": 0, "unclassified": 0, "above_genus": 0,
              "correct_genus": 0, "wrong_genus": 0}

    with open(out_correct, "w") as fout_correct:
        for path in chunk_paths:
            with open(path) as fin:
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

# find chunk .shrunk.lca files and group by species
def discover_chunk_groups(lca_dir):

    pattern = re.compile(r'^(\d+(_modern)?)_chunk\d+\.shrunk\.lca$')
    groups = {}
    for root, dirs, fnames in os.walk(lca_dir):
        for fname in sorted(fnames):
            m = pattern.match(fname)
            if m:
                base = m.group(1)
                groups.setdefault(base, []).append(os.path.join(root, fname))
    return {k: sorted(v) for k, v in sorted(groups.items())}

# main functions
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--lca_dir", required=True,
                        help="Directory containing bamdam chunk .shrunk.lca files (searched recursively)")
    parser.add_argument("--nodes", required=True,
                        help="Path to nodes.dmp from taxonomy")
    parser.add_argument("--out_dir", required=True,
                        help="Directory to write filtered output files")
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    taxonomy = load_taxonomy(args.nodes)

    groups = discover_chunk_groups(args.lca_dir)
    if not groups:
        print(f"ERROR: no chunk .shrunk.lca files found under {args.lca_dir}", file=sys.stderr)
        sys.exit(1)

    print(f"Found {len(groups)} species to process.", flush=True)

    for base, chunk_paths in groups.items():
        out_correct = os.path.join(args.out_dir, f"{base}.shrunk.correct_genus.lca")

        if os.path.isfile(out_correct):
            print(f"  {base} -> already done, skipping.", flush=True)
            continue

        species_id = int(base.replace("_modern", ""))
        print(f"  {base} ({len(chunk_paths)} chunks) -> {base}.shrunk.correct_genus.lca", flush=True)
        counts = filter_chunks(chunk_paths, out_correct, taxonomy, species_id)

        total   = counts["total"]
        correct = counts["correct_genus"]
        pct_correct = 100 * correct / total if total > 0 else 0
        print(f"    correct_genus: {correct:,} / {total:,} reads ({pct_correct:.3f}%)", flush=True)


if __name__ == "__main__":
    main()
