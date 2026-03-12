"""
Off-target alignment analysis.

For each read in a SAM file:
  1. Collect all mapped alignments; compute similarity = (aligned_len - NM) / aligned_len
  2. Map each reference to its ancestor at --level (default: genus) via acc2taxid + nodes.dmp
  3. Keep only the best-similarity alignment per taxon at that level
  4. Remove the group belonging to the target taxid's ancestor at that level
  5. Report how many reads have off-target hits above the similarity threshold,
     and which off-target taxa are the most frequent best alternative hits

Usage:
    python3 primary_secondary_analysis.py \\
        --sam           input.sam[.gz] | - \\
        --target_taxid  128015 \\
        --acc2tax       seqid2taxid.acc2tax \\
        --nodes         nodes.dmp \\
        --names         names.dmp \\
        --out_summary   output.txt \\
        --level         genus \\
        --simscorelow   0.95 \\
        --top           20
"""

import argparse
import contextlib
import gzip
import itertools
import re
import subprocess
import sys
from collections import Counter, defaultdict


# stream BAM files as SAM text via samtools view -h, or open SAM/SAM.gz files directly
class _SamtoolsView:
    def __init__(self, path):
        self.path = path
        self.proc = None

    def __enter__(self):
        self.proc = subprocess.Popen(
            ['samtools', 'view', '-h', self.path],
            stdout=subprocess.PIPE, text=True
        )
        return self.proc.stdout

    def __exit__(self, *args):
        if self.proc:
            self.proc.stdout.close()
            self.proc.wait()


def open_sam(path):
    if path == '-':
        return contextlib.nullcontext(sys.stdin)
    if path.endswith('.bam'):
        return _SamtoolsView(path)
    if path.endswith('.gz'):
        return gzip.open(path, 'rt')
    return open(path, 'rt')


# compute similarity = (aligned_length - NM) / aligned_length
def cigar_aln_len(cigar):
    if cigar == '*':
        return 0
    return sum(int(n) for n, op in re.findall(r'(\d+)([MIDNSHP=X])', cigar)
               if op in 'MI=X')


def get_nm(tag_fields):
    for tag in tag_fields:
        if tag.startswith('NM:i:'):
            return int(tag[5:])
    return None


def compute_similarity(cigar, nm):
    aln_len = cigar_aln_len(cigar)
    if aln_len == 0 or nm is None:
        return None
    return (aln_len - nm) / aln_len


# accession extraction from RNAME
def extract_accession(rname):
    if '|' not in rname:
        return rname
    parts = [p for p in rname.split('|') if p]
    for p in reversed(parts):
        if re.match(r'^[A-Z]{1,8}_?\d+(\.\d+)?$', p):
            return p
    return parts[-1]


# SAM scan: best similarity per (qname, rname) pair
def scan_sam(sam_path):
    """
    Single pass through the SAM file.
    For each (read, reference) pair, keep the best similarity seen.

    Returns:
        best_by_rname : dict  qname -> {rname -> sim (float or None)}
        rnames        : set   all unique reference names hit
    """
    best_by_rname = defaultdict(dict)
    rnames = set()
    n_lines = 0

    print("[Scan] Reading SAM ...", file=sys.stderr)
    with open_sam(sam_path) as f:
        for line in f:
            if line.startswith('@'):
                continue
            cols = line.split('\t')
            if len(cols) < 6:
                continue
            flag = int(cols[1])
            if flag & 4:      # unmapped
                continue
            if flag & 2048:   # supplementary
                continue

            qname = cols[0]
            rname = cols[2]
            if rname == '*':
                continue

            cigar = cols[5]
            nm    = get_nm(cols[11:]) if len(cols) > 11 else None
            sim   = compute_similarity(cigar, nm)

            rnames.add(rname)
            prev = best_by_rname[qname].get(rname)
            if prev is None or (sim is not None and (prev is None or sim > prev)):
                best_by_rname[qname][rname] = sim

            n_lines += 1
            if n_lines % 5_000_000 == 0:
                print(f"  ... {n_lines:,} lines, {len(best_by_rname):,} reads",
                      file=sys.stderr)

    print(f"  Scan done: {n_lines:,} alignment lines, "
          f"{len(best_by_rname):,} reads, "
          f"{len(rnames):,} unique references", file=sys.stderr)
    return dict(best_by_rname), rnames


# acc2taxid lookup
def load_acc2taxid(acc2tax_path, rnames):
    acc_to_rnames = defaultdict(list)
    for rname in rnames:
        acc = extract_accession(rname)
        acc_to_rnames[acc].append(rname)
        bare = acc.rsplit('.', 1)[0]
        if bare != acc:
            acc_to_rnames[bare].append(rname)

    needed = set(acc_to_rnames)
    rname_to_taxid = {}

    print(f"[acc2taxid] looking up {len(rnames):,} references "
          f"({len(needed):,} accession strings) ...", file=sys.stderr)

    opener = gzip.open if acc2tax_path.endswith('.gz') else open
    with opener(acc2tax_path, 'rt') as f:
        first = next(f)
        lines = f if first.startswith('accession') else itertools.chain([first], f)
        for i, line in enumerate(lines, 1):
            if not needed:
                break
            if i % 50_000_000 == 0:
                print(f"  ... {i:,} lines, {len(needed):,} remaining", file=sys.stderr)
            cols = line.split('\t')
            if len(cols) < 3:
                continue
            for acc in (cols[1], cols[0]):   # versioned first, then bare
                if acc in needed:
                    for rn in acc_to_rnames[acc]:
                        rname_to_taxid[rn] = cols[2].strip()
                    needed.discard(acc)

    if needed:
        print(f"[acc2taxid] WARNING: {len(needed):,} accession(s) not resolved "
              f"(first 5: {list(needed)[:5]})", file=sys.stderr)
    print(f"[acc2taxid] matched {len(rname_to_taxid):,} / {len(rnames):,} references",
          file=sys.stderr)
    return rname_to_taxid



# loading taxonomy information
def load_nodes(nodes_path):
    nodes = {}
    with open(nodes_path) as f:
        for line in f:
            parts = line.split('\t|\t')
            if len(parts) < 3:
                continue
            nodes[parts[0].strip()] = (parts[1].strip(), parts[2].strip())
    print(f"[nodes] loaded {len(nodes):,} nodes", file=sys.stderr)
    return nodes


def load_names(names_path, taxids):
    names = {}
    needed = {str(t) for t in taxids}
    with open(names_path) as f:
        for line in f:
            parts = line.split('\t|\t')
            if len(parts) < 4:
                continue
            tid = parts[0].strip()
            if tid in needed and 'scientific name' in parts[3]:
                names[tid] = parts[1].strip()
    return names


def ancestor_at_level(taxid, level, nodes, cache):
    tid = str(taxid)
    if tid in cache:
        return cache[tid]

    path = []
    cur  = tid
    while True:
        if cur in cache:
            result = cache[cur]
            break
        entry = nodes.get(cur)
        if entry is None:
            result = None
            break
        parent, rank = entry
        if rank == level:
            result = cur
            break
        if cur == parent:
            result = None
            break
        path.append(cur)
        cur = parent

    for t in path:
        cache[t] = result
    cache[tid] = result
    return result


# analysis logic
def analyse(best_by_rname, rname_to_taxid, target_anc, nodes, level, simscorelow):
    """
    For each read:
      - rname -> taxid -> level ancestor
      - keep best sim per level ancestor
      - skip target ancestor
      - count off-target ancestors whose best sim >= simscorelow

    Returns:
        off_target_counter : Counter  ancestor_taxid -> n_reads
        n_reads_any        : int      reads with >= 1 qualifying off-target hit
        n_reads_total      : int
    """
    cache            = {}
    off_target_ctr   = Counter()
    n_reads_any      = 0
    n_reads_total    = len(best_by_rname)

    for qname, rname_sims in best_by_rname.items():
        # collapse to best sim per level ancestor, excluding target
        best_per_anc = {}
        for rname, sim in rname_sims.items():
            taxid = rname_to_taxid.get(rname)
            if taxid is None:
                continue
            anc = ancestor_at_level(taxid, level, nodes, cache)
            if anc is None or anc == target_anc:
                continue
            prev = best_per_anc.get(anc)
            if prev is None or (sim is not None and (prev is None or sim > prev)):
                best_per_anc[anc] = sim

        # count off-target ancestors above threshold
        found = False
        for anc, sim in best_per_anc.items():
            if sim is not None and sim >= simscorelow:
                off_target_ctr[anc] += 1
                found = True
        if found:
            n_reads_any += 1

    return off_target_ctr, n_reads_any, n_reads_total


# main functions
def main():
    ap = argparse.ArgumentParser(
        description='Off-target alignment analysis: best hit per taxonomic level'
    )
    ap.add_argument('--sam',          required=True,
                    help='SAM / SAM.gz / BAM / - (stdin, text SAM)')
    ap.add_argument('--target_taxid', required=True,
                    help='Target taxid; its level-ancestor will be excluded')
    ap.add_argument('--acc2tax',      required=True,
                    help='acc2taxid file (nucl_gb.accession2taxid or seqid2taxid.acc2tax)')
    ap.add_argument('--nodes',        required=True, help='nodes.dmp')
    ap.add_argument('--names',        required=True, help='names.dmp')
    ap.add_argument('--out_summary',  required=True, help='Output summary text file')
    ap.add_argument('--level',        default='genus',
                    help='Taxonomic level to collapse hits to (default: genus)')
    ap.add_argument('--simscorelow',  type=float, default=0.95,
                    help='Similarity threshold: (aln_len - NM) / aln_len (default: 0.95)')
    ap.add_argument('--top',          type=int, default=20,
                    help='Number of top off-target taxa to report (default: 20)')
    args = ap.parse_args()

    target_taxid = str(args.target_taxid)

    # Step 1: scan SAM
    best_by_rname, rnames = scan_sam(args.sam)

    # Step 2: acc2taxid
    rname_to_taxid = load_acc2taxid(args.acc2tax, rnames)

    # Step 3: taxonomy — find the level ancestor of the target
    nodes = load_nodes(args.nodes)
    anc_cache    = {}
    target_anc   = ancestor_at_level(target_taxid, args.level, nodes, anc_cache)
    if target_anc is None:
        print(f"WARNING: no {args.level} ancestor found for {target_taxid}; "
              f"using exact taxid match", file=sys.stderr)
        target_anc = target_taxid
    print(f"[taxonomy] target {target_taxid} -> {args.level} {target_anc}",
          file=sys.stderr)

    # Step 4: analyse
    off_target_ctr, n_any, n_total = analyse(
        best_by_rname, rname_to_taxid, target_anc, nodes, args.level, args.simscorelow
    )

    # Step 5: names
    top_taxa = off_target_ctr.most_common(args.top)
    names    = load_names(args.names, [t for t, _ in top_taxa] + [target_anc, target_taxid])

    # Step 6: write summary
    lines = []
    lines.append("Off-target Alignment Summary")
    lines.append(f"  Input           : {args.sam}")
    lines.append(f"  Target taxid    : {target_taxid}  "
                 f"({names.get(target_taxid, 'unknown')})")
    lines.append(f"  Target {args.level:<10}: {target_anc}  "
                 f"({names.get(target_anc, 'unknown')})")
    lines.append(f"  Level           : {args.level}")
    lines.append(f"  Similarity      : (aligned_len - NM) / aligned_len  "
                 f"threshold = {args.simscorelow}")
    lines.append("")
    lines.append(f"  Total reads with at least one mapped alignment : {n_total:>10,}")
    lines.append(f"  Reads with >= 1 off-target {args.level} hit")
    lines.append(f"    with sim >= {args.simscorelow}                  : {n_any:>10,}  "
                 f"({100 * n_any / max(n_total, 1):.1f}%)")
    lines.append("")
    lines.append(f"Top {args.top} off-target {args.level}-level taxa "
                 f"(best hit per {args.level} per read, sim >= {args.simscorelow})")
    lines.append(f"  {'Rank':>4}  {'Taxid':>10}  {'Reads':>10}  {'%Total':>7}  Name")
    lines.append("  " + "-" * 72)
    for rank, (taxid, count) in enumerate(top_taxa, 1):
        pct  = 100 * count / max(n_total, 1)
        name = names.get(str(taxid), 'unknown')
        lines.append(f"  {rank:>4}  {taxid:>10}  {count:>10,}  {pct:>6.1f}%  {name}")

    summary = "\n".join(lines) + "\n"
    with open(args.out_summary, 'w') as f:
        f.write(summary)
    print(summary)
    print(f"[Done] {args.out_summary}", file=sys.stderr)


if __name__ == '__main__':
    main()
