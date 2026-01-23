# retreat-probe-design
Github for Part II Project

List of taxa that may respond to glacier retreat is synthesized from recent literature (see `taxids.txt`). Taxanomic ID were extracted from NCBI. 
Using NCBI CLI, a search on the availability of the full genomes for these taxa is conducted:
```
datasets summary genome taxon --inputfile taxids.txt --as-json-lines --limit all > genomes.jsonl
```
This is then converted to tsv file:
```
dataformat tsv genome --inputfile genomes.jsonl --fields organism-tax-id,organism-name,accession,source_database,assminfo-level,assminfo-refseq-category,assminfo-status,assminfo-release-date > genomes.tsv
```
Sample outputs are shown in `genomes.tsv`. 

From all identified genome, the ones with best quality is selected based on the following hierarchical priority using `availability_check.py`:
```
reference genome (any assembly level, any database)
representative genome (any assembly level, any database)
complete genome, refseq 
chromosome, refseq 
complete genome, genbank
chromosome, genbank
contig, refseq
contig, genbank
other
``` 
Those with no genome are left blank. Certain species have their subspecies/strains identified under taxonomic IDs different from the original ones, and are also included. 
Sample outputs are shown in `taxid_genome_availability.tsv`

