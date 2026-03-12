"""
Randomly subsample N reads from a FASTA file.

Usage:
    python3 subsample_fasta.py --input in.fasta --output out.fasta --n 300000
    python3 subsample_fasta.py --input in.fasta --output out.fasta --n 300000 --seed 42
"""

import argparse
import random
import sys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--input',  required=True, help='Input FASTA file')
    ap.add_argument('--output', required=True, help='Output FASTA file')
    ap.add_argument('--n',      required=True, type=int, help='Number of reads to sample')
    ap.add_argument('--seed',   type=int, default=None, help='Random seed (optional)')
    args = ap.parse_args()

    reads = []
    with open(args.input) as f:
        header = None
        seq = []
        for line in f:
            line = line.rstrip()
            if line.startswith('>'):
                if header is not None:
                    reads.append((header, ''.join(seq)))
                header = line
                seq = []
            else:
                seq.append(line)
        if header is not None:
            reads.append((header, ''.join(seq)))

    total = len(reads)
    print(f"Total reads: {total:,}", file=sys.stderr)

    if args.n >= total:
        print(f"WARNING: requested N={args.n:,} >= total {total:,}; writing all reads.",
              file=sys.stderr)
        sampled = reads
    else:
        random.seed(args.seed)
        sampled = random.sample(reads, args.n)

    with open(args.output, 'w') as f:
        for header, seq in sampled:
            f.write(header + '\n')
            f.write(seq + '\n')

    print(f"Written: {len(sampled):,} reads -> {args.output}", file=sys.stderr)


if __name__ == '__main__':
    main()
