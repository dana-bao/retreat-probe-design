import json
import pandas as pd

file_name = input("Enter taxonomy jsonl file name: ")
output_name = input("Enter output tsv file name: ")
rows = []

with open(file_name) as f:
    for line in f:
        rec = json.loads(line)

        tax = rec["taxonomy"]
        classification = tax.get("classification", {})

        species_taxid = tax["tax_id"]
        species_name = tax["current_scientific_name"]["name"]

        genus = classification.get("genus", {})
        genus_taxid = genus.get("id")
        genus_name = genus.get("name")

        rows.append({
            "species_taxid": species_taxid,
            "species_name": species_name,
            "genus_taxid": genus_taxid,
            "genus_name": genus_name
        })

df = pd.DataFrame(rows)
df.to_csv(output_name, sep="\t", index=False)
