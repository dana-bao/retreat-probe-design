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
This is done using `taxa_extraction.py` again. Output is shown in `species_rep_to_genus.tsv`.  

Using these following files:
*species genome metadata table
*genus genome metadata table
*species taxaID list
*species of interest to genus taxaID table
*species representative to genus taxaID table

`availability_check.py` try identify the best available assemblies for each species of interest.  
The script first identifies species level assemblies for species of interest. For those that do not have any assembly available, it falls back to genus level, then identifies representatives species assemblies within the genus.  
From all identified assemblies, the one with best quality is selected based on the following hierarchical priority:
```
refseq database and reference genome, by assembly level (complete genome -> chromosome -> scaffold -> contig)
genbank database and reference genome, by assembly level
refseq database, not reference genome, by assembly level
``` 
If none of the standard met or multiple assemblies with the same hierarchy appear for one species, choose the assembly with largest total sequence length.  
Those with no genome are left blank. Certain species have their subspecies/strains identified under taxonomic IDs different from the original ones, and are also included.  
Outputs are shown in `all_availability.tsv`.  

Genomes are downloaded according to assembly accession using `download_assemblies.sh`.  

**Ancient DNA Read Simulation**
To simulate ancient DNA reads from modern assemblies, `aeDNA_simulation.py` is developed and used. The program allows the following arguments:
* `--assembly`: the name of genome assembly file
* `--metadata`: the name of metadata table file
* `--output_dir`: the desired output directory, default current directory
* `--tile_len`: the length of each tile, default 52bp
* `--step`: the window sliding step between each tile, default 5bp
* `--deamination`: adding deamination damage that turns C to T at the firt 3bp of the 5' end of each tiled read
* `--mutation`: adding 1bp of mutation at the center of the read (e.g. bp 26 for the default 52bp length)
* `--seed`: setting randomize seed, only when testing
* `--wrap`: output fasta file wrapping, suitable for long tiled reads, default 80bp

for the purpose of this project, parameters used are as follow:
```
python3 aeDNA_simulation.py \
  --assembly assembly_accession_genomic.fna \
  --metadata all_availability2.tsv \
  --deamination \
  --mutation
```

To allow for parallel processing, a list of assemblies to process is created via:
``` 
find ../data/assemblies/ncbi_dataset/data -type f -name "*_genomic.fna" | sort > assemblies.list
wc -l assemblies.list
head assemblies.list
``` 

As there is limited storage, the assemblies are split into different batches for tiling. Assemblies selected for each batch is based on total_seq_length and test runs.  
This is then run with `slurm_tile.sbatch` using:
```
BATCH=batches/batch_XXX.list
N=$(wc -l < "$BATCH")
sbatch --array=1-"$N"%6 slurm_tile.sbatch "$BATCH"
```
where `batch_XXX.list` contains the list of assemblies directories for this particular batch.  
Sample tiling output is in `1118155.fasta`






cd hit
kraken2

The running process is then done using `run_tiling_array.sbatch`.  