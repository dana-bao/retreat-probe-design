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
From all identified assemblies, the ones with best quality is selected based on the following hierarchical priority:
```
refseq database and reference genome, by assembly level (complete genome -> chromosome -> scaffold -> contig)
genbank database and reference genome, by assembly level
refseq database, not reference genome, by assembly level
otherwise, choose the assembly with largest total sequence length
``` 
Those with no genome are left blank. Certain species have their subspecies/strains identified under taxonomic IDs different from the original ones, and are also included.  
Outputs are shown in `all_availability.tsv`.  

Genomes are downloaded according to assembly accession using `download_assemblies.sh`.  

**Ancient DNA Read Simulation**
To simulate ancient DNA reads from modern assemblies, gargammel is used.  

Reorganize the downloaded assemblies using:
``` 
find ncbi_dataset/data -type f -name "*_genomic.fna" \
  -exec cat {} + > all_assemblies.fna
``` 
This will be the endogenous input that will be processed by gargammel to simulate aeDNA reads.  
However, gargammel is designed to also consider modern species contamination and modern microbial contamination. Therefore, the script require additional input for these. Dummy scripts are thus created as a filler for both: 
```
>dummy
TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
```
Gargammel is subsequently ran with the following parameters:
```
perl ./gargammel.pl -c 1 -f ./src/sizefreq.size.gz -matfile ./src/matrices/single- -o ./output/output ./ 
```
* `-c`: coverage
* `-f`: file containing list of frequency of different fragment sizes
* `-matfile`: matrix file containing substitution misincorporation due to deamination
* `-o`: specify output directory and file name
As empirical evidence of for frequency and deamination of aeDNA in glacial retreat context is lacking, precalculated values provided by gargammel based on chosen studies are used (`sizefreq.size.gz` and `src/matrices/`).  

Gargammel read simulation adopts a structure similar to sequencing for generating reads from reference genome, instead of using deterministic, tiling approach. Therefore, to ensure the reads to cover as much of the full assemblies as possible without compromising memory space and efficiency, test runs with different coverage is done with one specific assembly, `GCA_000239015.2`.  
The file sizes of input and respective output files are as follow:
``` 
input file:
12M Jan 30 16:02 GCA_000239015.2_ASM23901v2_genomic.fna

output files:
8.7M Jan 30 16:09 c1_test_output_s1.fq.gz
9.0M Jan 30 16:09 c1_test_output_s2.fq.gz
18M Jan 30 16:29 c2_test_output_s1.fq.gz
18M Jan 30 16:29 c2_test_output_s2.fq.gz
26M Jan 30 16:54 c3_test_output_s1.fq.gz
27M Jan 30 16:54 c3_test_output_s2.fq.gz
44M Jan 30 17:03 c5_test_output_s1.fq.gz
45M Jan 30 17:03 c5_test_output_s2.fq.gz
87M Jan 30 17:09 c10_test_output_s1.fq.gz
90M Jan 30 17:09 c10_test_output_s2.fq.gz
``` 
From the test runs, it can be inferred that the output file size shares a relatively linear relationship with coverages. This is later used as a reference to infer the appropriate coverage to avoid storage issue.  
To acquire the optimal coverage, we assume that read assignment follows a poisson distribution and calculate as follow: 
$$
C = - \ln\!\left(1 - P^{L/G}\right)
$$
Where L is average length of reads, G is genome assembly length, and P is the probability of genome covered entirely.  
Longest assemblies are sorted using `long_assemblies_sort.sh`. Outputs are as follow:
``` 
8030382961	GCA_030686995.1
3593360174	GCA_037465115.1
2738132926	GCA_052724335.1
2715530335	GCF_003573695.1
2416637544	GCA_032207245.1
2307496354	GCA_033807815.1
1815982469	GCA_965199645.1
1791703142	GCA_026214975.1
1779762036	GCA_026929855.1
1631492059	GCA_046244935.1
``` 
