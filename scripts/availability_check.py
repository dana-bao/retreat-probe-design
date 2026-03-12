import pandas as pd

# File loadings
species_file_name = input("Enter species genome metadata TSV file name: ")
genus_file_name = input("Enter genus genome metadata TSV file name: ")
list_name = input("Enter species taxaID list file name (one taxid per line): ")
species_to_genus_name = input("Enter species-to-genus taxaID table TSV file name: ")
species_rep_to_genus_name = input("Enter species-rep-to-genus taxaID table TSV file name: ")
output_name = input("Enter output TSV file name: ")

species_df = pd.read_csv(species_file_name, sep="\t")
genus_df = pd.read_csv(genus_file_name, sep="\t")
species_list = pd.read_csv(list_name, sep="\t", header=None, names=["Species Taxonomic ID"])
species_to_genus = pd.read_csv(species_to_genus_name, sep="\t")
species_rep_to_genus = pd.read_csv(species_rep_to_genus_name, sep="\t")

# Use classification logic to assign a priority and status describing the assembly we have
def classify_row(row):
    # refseq and reference genome, by assembly level
    if row["Source Database"] == "SOURCE_DATABASE_REFSEQ" and row["Assembly Refseq Category"] == "reference genome" and row["Assembly Level"] == "Complete Genome":
        return 1, "RefSeq reference genome, complete genome"
    if row["Source Database"] == "SOURCE_DATABASE_REFSEQ" and row["Assembly Refseq Category"] == "reference genome" and row["Assembly Level"] == "Chromosome":
        return 2, "RefSeq reference genome, chromosome"
    if row["Source Database"] == "SOURCE_DATABASE_REFSEQ" and row["Assembly Refseq Category"] == "reference genome" and row["Assembly Level"] == "Scaffold":
        return 3, "RefSeq reference genome, scaffold"
    if row["Source Database"] == "SOURCE_DATABASE_REFSEQ" and row["Assembly Refseq Category"] == "reference genome" and row["Assembly Level"] == "Contig":
        return 4, "RefSeq reference genome, contig"
    # genbank and reference genome, by assembly level
    if row["Source Database"] == "SOURCE_DATABASE_GENBANK" and row["Assembly Refseq Category"] == "reference genome" and row["Assembly Level"] == "Complete Genome":
        return 5, "Genbank reference genome, complete genome"
    if row["Source Database"] == "SOURCE_DATABASE_GENBANK" and row["Assembly Refseq Category"] == "reference genome" and row["Assembly Level"] == "Chromosome":
        return 6, "Genbank reference genome, chromosome"
    if row["Source Database"] == "SOURCE_DATABASE_GENBANK" and row["Assembly Refseq Category"] == "reference genome" and row["Assembly Level"] == "Scaffold":
        return 7, "Genbank reference genome, scaffold"
    if row["Source Database"] == "SOURCE_DATABASE_GENBANK" and row["Assembly Refseq Category"] == "reference genome" and row["Assembly Level"] == "Contig":
        return 8, "Genbank reference genome, contig"

    # representative genome (any assembly level, any database) - have not appear in any of the metadata for species involved in this project
    if row["Assembly Refseq Category"] == "representative genome":
        return 9, "representative genome"

    # refseq by assembly level
    if row["Source Database"] == "SOURCE_DATABASE_REFSEQ" and row["Assembly Level"] == "Complete Genome":
        return 10, "RefSeq complete genome"
    if row["Source Database"] == "SOURCE_DATABASE_REFSEQ" and row["Assembly Level"] == "Chromosome":
        return 11, "RefSeq chromosome"
    if row["Source Database"] == "SOURCE_DATABASE_REFSEQ" and row["Assembly Level"] == "Scaffold":
        return 12, "RefSeq scaffold"
    if row["Source Database"] == "SOURCE_DATABASE_REFSEQ" and row["Assembly Level"] == "Contig":
        return 13, "RefSeq contig"

    # Otherwise, compare largest total sequence length later
    return 99, "Largest total sequence length"

# Applying classification logic and generate df containing them
def add_priority_and_status(df):
    df[["Priority", "Status"]] = df.apply(classify_row, axis=1, result_type="expand")
    return df

# Selection of best assembly/assemblies based on assigned priority or total sequence length
def select_best_per_taxid_keep_ties(df, taxid_col):
    df = df.copy()

    # Determine for each taxid the min priority
    min_pri = df.groupby(taxid_col)["Priority"].transform("min")
    has_preferred = min_pri < 99

    # Priority preference winners
    keep_preferred = has_preferred & (df["Priority"] == min_pri)

    # Total length winners, only when no priority preference match exists
    max_len = df.groupby(taxid_col)["Assembly Stats Total Sequence Length"].transform("max")
    keep_by_len = (~has_preferred) & (df["Assembly Stats Total Sequence Length"] == max_len)

    best = df[keep_preferred | keep_by_len].copy()

    # Fill in Status for total length winners
    best.loc[best["Priority"] == 99, "Status"] = "Largest total sequence length"

    return best



# Compute Species-level best assemblies
species_df = add_priority_and_status(species_df)

species_best = select_best_per_taxid_keep_ties(species_df, taxid_col="Organism Taxonomic ID")

species_best = species_best.rename(columns={
    "Organism Taxonomic ID": "Species Taxonomic ID",
    "Organism Name": "Species Name",
})

# Merge to species taxid list, keeping right_only subspecies/strain assembly
base = species_list.merge(
    species_best,
    on="Species Taxonomic ID",
    how="outer",
    indicator=True
)

base["Species level Assembly"] = base["_merge"].map({
    "both": "yes",
    "right_only": "yes",
    "left_only": "no",
})

# Genus fallback for species_list taxids with no species-level assembly (left_only)
need_genus = base["_merge"] == "left_only"

# Attach genus_taxid, genus_name and species_name from mapping
mapping = species_to_genus.rename(columns={
    "species_taxid": "Species Taxonomic ID",
    "species_name": "Species Name (from mapping)",
    "genus_taxid": "Genus Taxonomic ID",
    "genus_name": "Genus Name",
})

base = base.merge(
    mapping[["Species Taxonomic ID", "Species Name (from mapping)", "Genus Taxonomic ID", "Genus Name"]],
    on="Species Taxonomic ID",
    how="left",
    suffixes=("", "_map")
)

# Fill Species Name from mapping
base.loc[base["Species Name"].isna(), "Species Name"] = base.loc[base["Species Name"].isna(), "Species Name (from mapping)"]
base = base.drop(columns=["Species Name (from mapping)"])

# Compute Genus-level best assemblies
genus_df = add_priority_and_status(genus_df)

rep_map = species_rep_to_genus.rename(columns={
    "species_taxid": "Organism Taxonomic ID",
    "genus_taxid": "Genus Taxonomic ID",
    "genus_name": "Genus Name",
})[["Organism Taxonomic ID", "Genus Taxonomic ID", "Genus Name"]]

genus_df = genus_df.merge(
    rep_map,
    on="Organism Taxonomic ID",
    how="left"
)

genus_best = select_best_per_taxid_keep_ties(genus_df, taxid_col="Genus Taxonomic ID")

# Merge genus_best for rows needing genus fallback
left_only_rows = base[need_genus].copy()
other_rows = base[~need_genus].copy()
left_only_rows["Genus Taxonomic ID"] = left_only_rows["Genus Taxonomic ID"].astype("Int64") # Used because subspecies/strains do not have genus level and would cause datatype issue
genus_best["Genus Taxonomic ID"] = genus_best["Genus Taxonomic ID"].astype("Int64")

overlap = [
    "Assembly Accession",
    "Source Database",
    "Assembly Level",
    "Assembly Refseq Category",
    "Assembly Status",
    "Assembly Release Date",
    "Status",
    "Priority",
    "Assembly Stats Total Sequence Length",
    "Organism Name",
]

left_only_rows = left_only_rows.drop(columns=[c for c in overlap if c in left_only_rows.columns])

left_only_rows = left_only_rows.merge(
    genus_best[[
        "Genus Taxonomic ID",
        "Genus Name",
        "Organism Name",
        "Assembly Accession",
        "Source Database",
        "Assembly Level",
        "Assembly Refseq Category",
        "Assembly Status",
        "Assembly Release Date",
        "Status",
        "Priority",
        "Assembly Stats Total Sequence Length",
    ]],
    on="Genus Taxonomic ID",
    how="left",
)

final = pd.concat([other_rows, left_only_rows], ignore_index=True)

# if a species still has more than 1 assemblies, keep the one with the longest total length
final["Assembly Stats Total Sequence Length"] = pd.to_numeric(
    final["Assembly Stats Total Sequence Length"],
    errors="coerce"
)

final = (
    final.sort_values(
        by=["Species Taxonomic ID", "Assembly Stats Total Sequence Length"],
        ascending=[True, False]
    )
    .drop_duplicates(subset=["Species Taxonomic ID"], keep="first")
    .reset_index(drop=True)
)


final = final[[
    "Species Taxonomic ID",
    "Species Name",
    "Species level Assembly",
    "Genus Taxonomic ID",
    "Genus Name",
    "Assembly Accession",
    "Source Database",
    "Assembly Level",
    "Assembly Refseq Category",
    "Assembly Status",
    "Assembly Release Date",
    "Assembly Stats Total Sequence Length",
    "Status",
    "_merge",
]]

final.to_csv(output_name, sep="\t", index=False)
print(f"Done writing {output_name}")
