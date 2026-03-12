'''
Check which target species have representation in both ngsLCA and Kraken2 databases.

For species with a species-level assembly, the assembly organism taxid is the species taxid.
For species using a same-genus representative assembly, the assembly organism taxid is looked up
from genus_metadata.tsv by matching the assembly accession.

If --in_ngslca and --in_kraken are not provided, produces organism_taxid_mapping.tsv
  and organism_taxids_to_check.txt for use in cluster-side database checks
  (see pipeline_overview.md for database check commands).

If --in_ngslca and --in_kraken are provided, produces db_coverage.tsv with per-species
  database presence flags, and species_in_both_dbs.txt listing target taxids found
  in both databases.

Usage:
  # Without database results:
  python3 check_db_coverage.py \
    --all_avail data/all_availability.tsv \
    --genus_meta data/genus_metadata.tsv \
    --out_dir data/

  # With database results:
  python3 check_db_coverage.py \
    --all_avail data/all_availability.tsv \
    --genus_meta data/genus_metadata.tsv \
    --out_dir data/ \
    --in_ngslca data/in_ngslca.txt \
    --in_kraken data/in_kraken2.txt
'''

import argparse
import os
import pandas as pd

# building genus-level taxid mapping from all_availability.tsv and genus_metadata.tsv
def build_organism_taxid_map(all_avail_path, genus_meta_path):
    all_avail = pd.read_csv(all_avail_path, sep='\t', dtype=str)
    genus_meta = pd.read_csv(genus_meta_path, sep='\t', dtype=str)

    # map assembly accession -> organism taxid from genus metadata
    acc_to_org_taxid = dict(zip(
        genus_meta['Assembly Accession'],
        genus_meta['Organism Taxonomic ID']
    ))

    rows = []
    for _, row in all_avail.iterrows():
        target_taxid = str(row['Species Taxonomic ID']).strip()
        has_species_assembly = str(row.get('Species level Assembly', '')).strip().lower() == 'yes'
        accession = str(row.get('Assembly Accession', '')).strip()
        merge = str(row.get('_merge', '')).strip()

        if accession in ('', 'nan'):
            # no assembly available
            rows.append({
                'target_taxid': target_taxid,
                'organism_taxid': None,
                'assembly_accession': None,
                'note': 'no_assembly'
            })
            continue

        if has_species_assembly or merge == 'right_only':
            # species-level assembly: organism is the target species itself
            organism_taxid = target_taxid
        else:
            # genus-representative assembly: look up from genus metadata
            organism_taxid = acc_to_org_taxid.get(accession)

        rows.append({
            'target_taxid': target_taxid,
            'organism_taxid': organism_taxid,
            'assembly_accession': accession,
            'note': 'species_level' if (has_species_assembly or merge == 'right_only') else 'genus_representative'
        })

    return pd.DataFrame(rows)

# main functions
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--all_avail',   required=True, help='Path to all_availability.tsv')
    parser.add_argument('--genus_meta',  required=True, help='Path to genus_metadata.tsv')
    parser.add_argument('--out_dir',     required=True, help='Output directory')
    parser.add_argument('--in_ngslca',   default=None,  help='File of taxids found in nucl_gb.accession2taxid (step 3)')
    parser.add_argument('--in_kraken',   default=None,  help='File of taxids found in Kraken2 db (step 3)')
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    df = build_organism_taxid_map(args.all_avail, args.genus_meta)

    # Write full mapping table
    mapping_path = os.path.join(args.out_dir, 'organism_taxid_mapping.tsv')
    df.to_csv(mapping_path, sep='\t', index=False)
    print(f'Wrote organism taxid mapping -> {mapping_path}')

    # Write unique organism taxids to check (excluding no_assembly)
    to_check = df[df['organism_taxid'].notna()]['organism_taxid'].unique()
    check_path = os.path.join(args.out_dir, 'organism_taxids_to_check.txt')
    with open(check_path, 'w') as f:
        for t in sorted(to_check):
            f.write(t + '\n')
    print(f'Wrote {len(to_check)} unique organism taxids -> {check_path}')
    print(f'  ({df["note"].value_counts().to_dict()})')

    print('Copy organism_taxids_to_check.txt to the cluster, run the database checks (see pipeline_overview.md),')
    print('then rerun with --in_ngslca and --in_kraken')

    # Step 3: if results provided, generate final coverage table
    if args.in_ngslca and args.in_kraken:
        with open(args.in_ngslca) as f:
            in_ngslca = set(l.strip() for l in f if l.strip())
        with open(args.in_kraken) as f:
            in_kraken = set(l.strip() for l in f if l.strip())

        df['in_ngslca'] = df['organism_taxid'].apply(
            lambda t: t in in_ngslca if pd.notna(t) else False)
        df['in_kraken2'] = df['organism_taxid'].apply(
            lambda t: t in in_kraken if pd.notna(t) else False)
        df['in_both'] = df['in_ngslca'] & df['in_kraken2']

        coverage_path = os.path.join(args.out_dir, 'db_coverage.tsv')
        df.to_csv(coverage_path, sep='\t', index=False)
        print(f'\nWrote coverage table -> {coverage_path}')

        in_both = df[df['in_both']]['target_taxid'].tolist()
        in_both_path = os.path.join(args.out_dir, 'species_in_both_dbs.txt')
        with open(in_both_path, 'w') as f:
            for t in in_both:
                f.write(t + '\n')
        print(f'Species in BOTH databases: {len(in_both)} -> {in_both_path}')
        print(f'  in ngsLCA only: {df["in_ngslca"].sum() - df["in_both"].sum()}')
        print(f'  in Kraken2 only: {df["in_kraken2"].sum() - df["in_both"].sum()}')
        print(f'  in neither: {(~df["in_ngslca"] & ~df["in_kraken2"] & df["organism_taxid"].notna()).sum()}')
        print(f'  no assembly: {(df["note"] == "no_assembly").sum()}')


if __name__ == '__main__':
    main()
