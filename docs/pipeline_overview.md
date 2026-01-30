## Pipeline overview

This document contains step-by-step computational workflow used in this project, with relevant parameter choices.

**Taxa List Synthesis**
List of taxa that may respond to glacial retreat is synthesized from recent literature (see `taxids.txt`). Taxonomic ID at species level are extracted from NCBI.

**Genome Assemblies Acquisition**
Species taxonomic IDs are used to search for taxonomy information:
```
datasets summary taxonomy taxon --inputfile species_id.txt --as-json-lines > taxonomy.jsonl
```
This is then used to extract genus level Taxonomic ID using `taxa_extraction.py`. Output is shown in `species_to_genus.tsv`.  

Using NCBI CLI, a search on the availability of the genomes for these taxa at both species and genus levels is conducted:
```
datasets summary genome taxon --inputfile species_id.txt --as-json-lines --limit all > species_genomes.jsonl
datasets summary genome taxon --inputfile genus_id.txt --as-json-lines --limit all > genus_genomes.jsonl
```
These are then converted to metadata tsv file:
```
dataformat tsv genome --inputfile species_genomes.jsonl --fields organism-tax-id,organism-name,accession,source_database,assminfo-level,assminfo-refseq-category,assminfo-status,assminfo-release-date > species_metadata.tsv
dataformat tsv genome --inputfile genus_genomes.jsonl --fields organism-tax-id,organism-name,accession,source_database,assminfo-level,assminfo-refseq-category,assminfo-status,assminfo-release-date > genus_metadata.tsv
```
As genus metadata contains, in fact, the metadata of representative species within the genus, species of interest need to be mapped to the corresponding representative species based on genus. 
This is done using `taxa_extraction.py` again. Output is shown in `species_rep_to_genus.tsv`

Using these following files:
*species genome metadata table
*genus genome metadata table
*species taxaID list
*species of interest to genus taxaID table
*species representative to genus taxaID table

`availability_check.py` try identify the best available assemblies for each species of interest
The script first identifies species level assemblies for species of interest. For those that do not have any assembly available, it falls back to genus level, then identifies representatives species assemblies within the genus. 
From all identified assemblies, the ones with best quality is selected based on the following hierarchical priority:
```
refseq database and reference genome, by assembly level (complete genome -> chromosome -> scaffold -> contig)
genbank database and reference genome, by assembly level
refseq database, not reference genome, by assembly level
otherwise, choose the assembly with largest total sequence length
``` 
Those with no genome are left blank. Certain species have their subspecies/strains identified under taxonomic IDs different from the original ones, and are also included. 
Outputs are shown in `all_availability.tsv`

Genomes are downloaded according to assembly accession using `download_assemblies.sh`. 

**Ancient DNA Read Simulation**
To simulate ancient DNA reads from modern assemblies, gargammel is used. 

First, reorganize the downloaded assemblies using:
``` 
find ncbi_dataset/data -type f -name "*_genomic.fna" \
  -exec cat {} + > all_assemblies.fna
``` 


