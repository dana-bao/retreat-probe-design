from Bio import SeqIO
import argparse
import csv
import os
import random
import re
import sys


# Extract assembly accession from filename
def extract_assembly_accession(path: str) -> str:
    filename = os.path.basename(path)
    m = re.search(r"(GCA|GCF)_[0-9]{9}\.[0-9]+", filename)
    if not m:
        sys.exit(f"ERROR: Cannot extract assembly accession from filename: {filename}")
    return m.group(0)

# Search for species taxid using assembly accession
def load_taxid_from_metadata(metadata_tsv: str, accession: str) -> str:
    with open(metadata_tsv, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")

        if "Assembly Accession" not in reader.fieldnames or \
           "Species Taxonomic ID" not in reader.fieldnames:
            sys.exit("ERROR: Metadata missing required columns")

        for row in reader:
            if row["Assembly Accession"] == accession:
                return str(int(float(row["Species Taxonomic ID"])))

    sys.exit(f"ERROR: No taxid found for accession {accession}")

# Apply changes of deamination/mutation if required
DNA_bases = ("A", "C", "G", "T")
def apply_change(tile: list[str], tile_len: int, deamination: bool, mutation: bool, rng: random.Random) -> list[str]:
    n = len(tile)
    if n == 0:
        return tile

    # Deamination in first 3 bp
    if deamination:
        first3 = min(3, n)
        c_pos = [i for i in range(first3) if tile[i].upper() == "C"]
        if c_pos:
            pos = rng.choice(c_pos)
            tile[pos] = "T" if tile[pos].isupper() else "t"

    # One random mutation at the middle of the read (26 for full length 52bp reads, middle for shorter tile)
    if mutation:
        target_idx = (tile_len - 1) // 2
        actual_mid = (n - 1) // 2
        center_idx = min(target_idx, actual_mid)
        old = tile[center_idx]
        old_u = old.upper()
        choices = [b for b in DNA_bases if b != old_u] if old_u in DNA_bases else list(DNA_bases)
        new = rng.choice(choices)
        tile[center_idx] = new if old.isupper() else new.lower()

    return tile


# Main program
def main():
    # Taking in arguments
    ap = argparse.ArgumentParser()
    ap.add_argument("--assembly", required=True)
    ap.add_argument("--metadata", required=True)
    ap.add_argument("--output_dir", default=".")
    ap.add_argument("--tile_len", type=int, default=52)
    ap.add_argument("--step", type=int, default=5)
    ap.add_argument("--deamination", action="store_true",
                    help="Enable C to T deamination in first 3 bp")
    ap.add_argument("--mutation", action="store_true",
                    help="Enable one random mutation in positions 4–52")
    ap.add_argument("--seed", type=int, default=None)
    ap.add_argument("--wrap", type=int, default=80)

    args = ap.parse_args()

    # randomization
    rng = random.Random(args.seed)

    # taxid search
    accession = extract_assembly_accession(args.assembly)
    taxid = load_taxid_from_metadata(args.metadata, accession)

    # creating output directory and file
    os.makedirs(args.output_dir, exist_ok=True)
    out_path = os.path.join(
        args.output_dir,
        f"{taxid}.fasta"
    )

    tile_count = 0

    # Process assembly and write tiles to output
    with open(out_path, "w", encoding="utf-8") as out_fh:
        for record in SeqIO.parse(args.assembly, "fasta"):
            read_code = record.id
            seq = str(record.seq)
            L = len(seq)

            if L < args.tile_len:
                tile = list(seq)
                tile = apply_change(tile, args.tile_len, args.deamination, args.mutation, rng)

                tile_idx = 1
                tile_count += 1
                start_1 = 1
                end_1 = L
                header = f">{taxid}|{read_code}|tile{tile_idx}|pos{start_1}-{end_1}"
                out_fh.write(header + "\n")
                tile_seq = "".join(tile)

                if args.wrap > 0:
                    for i in range(0, len(tile_seq), args.wrap):
                        out_fh.write(tile_seq[i:i + args.wrap] + "\n")
                else:
                    out_fh.write(tile_seq + "\n")
                continue

            tile_idx = 0
            for start in range(0, L - args.tile_len + 1, args.step):
                end = start + args.tile_len
                tile = list(seq[start:end])

                tile = apply_change(tile, args.tile_len, args.deamination, args.mutation, rng)

                tile_idx += 1
                tile_count += 1
                start_1 = start + 1
                end_1 = end
                header = f">{taxid}|{read_code}|tile{tile_idx}|pos{start_1}-{end_1}"
                out_fh.write(header + "\n")
                tile_seq = "".join(tile)

                if args.wrap > 0:
                    for i in range(0, len(tile_seq), args.wrap):
                        out_fh.write(tile_seq[i:i + args.wrap] + "\n")
                else:
                    out_fh.write(tile_seq + "\n")

    print(
        f"Done\n"
        f"  Assembly accession: {accession}\n"
        f"  Species taxid:      {taxid}\n"
        f"  Tiles written:      {tile_count}\n"
        f"  Output:             {out_path}",
        file=sys.stderr
    )


if __name__ == "__main__":
    main()
