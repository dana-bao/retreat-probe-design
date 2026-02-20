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
* `--wrap`: output FASTA file wrapping, suitable for long tiled reads, default 80bp

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

All 84 assemblies were tiled in parallel using `run_tiling_array.sbatch`.  
Sample tiling output is in `1118155.fasta`.  

**Reassigning Reads with Taxonomic Identification Algorithms**
An existing, wildly adopted taxonomic identification tools, Kraken2, is first used. 
Kraken2 core_nt database is chosen as the reference database for the purpose of this project.  
Prebuild core_nt database is acquired from Kraken 2 index zone. It is subsequently run using `sbatch_mom2.sh` with the following relevant parameters: 
``` 
kraken2 \
  --db cort_nt \
  --threads 48 \
  --confidence 0.05 \
  --report-minimizer-data \
  --report species_id.c0.2.report.txt" \
  --output species_id.c0.2.kraken.out" \
``` 
* `--db`: reference database
* `--threads`: number of cpus exploited
* `--confidence`: threshold for fraction of k-mers supporting the classification, otherwise considered unclassified
* `--report-minimizer-data`: report minimizer and distinct minimizer count information
* `--report`: return summary report
* `--output`: return output

Sample Kraken2 output can be found in `kraken2_sample_output`.  

To evaluate the best confidence threshold, 4 different values (0, 0.05, 0.1, 0.2) are tested on 10 sample species selected to represent the diversity of glacial retreat taxa composition.  The resulting kraken2 outputs are compared on their sensitivity (percentage of correctly assigned reads out of all reads) and precision (percentage of correctly assigned reads out of all classified reads) using f1 score calculated in `kraken_eval.py`. The 10 sample species can be found in `sample_10_species.txt`, with evaluation results in `kraken_eval.csv`. For this speicfic project, the most ideal confidence threshold have been identified as 0.05 and is adopted for further processing.  

The Kraken2 outputs are filtered to only include reads with a correct genus level assignment using `kraken_filter.py`. Sample output can be found in `129212_task1.k2.0.05.core_nt.correct_genus.out`.  

In addition, an in-house competitive mapping pipeline is used to compare the quality of resulting probes.  
The pipeline requires bowtie2. To ensure that bowtie2 uses the same database as Kraken2, NCBI BLAST core_nt is acquired and converted to FASTA files via BLAST: 
``` 
update_blastdb.pl --decompress core_nt
blastdbcmd -db core_nt -entry all -out core_nt.fasta
``` 
Bowtie2 database is then built and indexed using `bowtie2_index.sbatch`.  
Alignment is done using `bowtie2_run.sbatch`, where each array task aligns one FASTA file against one database shard. An example command with relevant parameters for a single task is:
```
bowtie2 -p 16 -k 100 -x core_nt.00 -f -U species.fasta --no-unal 2> logfile | samtools view -@ 16 -b -o tmp_bam -
```
* `-p`: number of threads
* `-k`: maximum number of alignments to report per read
* `-x`: bowtie2 index prefix (one shard of core_nt)
* `-f`: input in FASTA format
* `-U`: unpaired input reads
* `--no-unal`: suppress unaligned reads in output
* `-@` (samtools): number of threads for BAM compression
* `-b` (samtools): output in BAM format
* `-o` (samtools): output file path

After bowtie2 alignment, the output BAM files are merged and name-sorted with `bamsort.sbatch`. For species with more aligned reads, the merged BAM header can exceed the BAM format limit; `bamsort_merge_sam.sbatch` is used for these species instead, streaming reads from each BAM individually with `samtools view` rather than `samtools merge`, and sorting with GNU sort, outputting gzip-compressed SAM to bypass the limit.  

ngsLCA is then run on the sorted files with `ngslca.sbatch`. An example command with relevant parameters for a single species is:
```
ngsLCA \
    -simscorelow 0.95 \
    -simscorehigh 1.0 \
    -fix-ncbi 0 \
    -names names.dmp \
    -nodes nodes.dmp \
    -acc2tax nucl_gb.accession2taxid \
    -bam species_id.merged.sorted.bam \
    -outnames species_id
```
* `-simscorelow`: minimum alignment similarity score to consider a read
* `-simscorehigh`: maximum alignment similarity score to consider a read
* `-fix-ncbi`: whether to apply NCBI-specific accession fixes (0 = off)
* `-names`: NCBI taxonomy names file
* `-nodes`: NCBI taxonomy nodes file
* `-acc2tax`: accession-to-taxid mapping file
* `-bam`: input name-sorted BAM or SAM file
* `-outnames`: prefix for output files

Bamdam is then run on the ngsLCA outputs with `bamdam.sbatch` to filter the BAM and LCA files down to reads assigned at genus level or below. An example command with relevant parameters is:
```
bamdam shrink \
    --in_bam species_id.merged.sorted.bam \
    --in_lca species_id.lca \
    --out_bam species_id.shrunk.bam \
    --out_lca species_id.shrunk.lca \
    --stranded ss \
    --upto genus
```
* `--in_bam`: input name-sorted BAM file
* `--in_lca`: input LCA file from ngsLCA
* `--out_bam`: output filtered BAM
* `--out_lca`: output filtered LCA
* `--stranded`: library strandedness (`ss` = single-stranded)
* `--upto`: taxonomic rank to filter reads up to

For species with sorted file in `.sam.gz` format, the script first converts it to BAM with `samtools view`, as bamdam requires BAM input.  

