import pandas as pd

# File loadings
file_name = input("Enter genome metadata file name: ")
df = pd.read_csv(file_name, sep="\t")
list_name = input("Enter taxaID list file name: ")
taxa_list = pd.read_csv(list_name, sep="\t", header=None, names=["Organism Taxonomic ID"])

# Classification logic
def classify_row(row):
    if row["Assembly Refseq Category"] == "reference genome":
        return 1, "reference genome"
    if row["Assembly Refseq Category"] == "representative genome":
        return 2, "representative genome"
    if row["Source Database"] == "SOURCE_DATABASE_REFSEQ" and row["Assembly Level"] == "Complete Genome":
        return 3, "Refseq complete genome"
    if row["Source Database"] == "SOURCE_DATABASE_REFSEQ" and row["Assembly Level"] == "Chromosome":
        return 4, "Refseq complete chromosome"
    if row["Source Database"] == "SOURCE_DATABASE_GENBANK" and row["Assembly Level"] == "Complete Genome":
        return 5, "Genbank complete genome"
    if row["Source Database"] == "SOURCE_DATABASE_GENBANK" and row["Assembly Level"] == "Chromosome":
        return 6, "Genbank complete chromosome"
    if row["Source Database"] == "SOURCE_DATABASE_REFSEQ" and row["Assembly Level"] == "Contig":
        return 7, "Refseq contig"
    if row["Source Database"] == "SOURCE_DATABASE_GENBANK" and row["Assembly Level"] == "Contig":
        return 8, "Genbank contig"
    return 9, "Other assembly"

# Apply classification
df[["Priority", "Status"]] = df.apply(classify_row, axis=1, result_type="expand")

# For each taxonomy ID, pick the best assemblies
best_assemblies = df[df["Priority"] == df.groupby("Organism Taxonomic ID")["Priority"].transform("min")].copy()

# Merge with original taxa list
total_assemblies = taxa_list.merge(
    best_assemblies,
    on="Organism Taxonomic ID",
    how="outer",
    indicator=True)

# Export the best_assemblies
total_assemblies.to_csv("taxid_genome_availability.tsv", sep="\t", index=False)

print("Done writing taxid_genome_availability.tsv")
