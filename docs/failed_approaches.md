This document contains approaches attempted but failed or deemed inappropriate for the purpose of this project. 

**Genome Assemblies Acquisition**

Acquiring only species level assembly without falling back to genus level was attempted, but later discarded due to the low number of assemblies acquired. 

**Ancient DNA Read Simulation**

To simulate ancient DNA reads from modern assemblies, gargammel was attempted.  
Reorganize the downloaded assemblies using:
``` 
find ncbi_dataset/data -type f -name "*_genomic.fna" \
  -exec cat {} + > all_assemblies.fna
``` 
This will be the endogenous input that will be processed by gargammel to simulate aeDNA reads.  
However, gargammel is designed to also consider modern species contamination and modern microbial contamination. Therefore, the script requires additional input for these. Dummy scripts are thus created as a filler for both: 
```
>dummy
TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
```
Gargammel is subsequently run with the following parameters:
```
perl ./gargammel.pl -c 1 -f ./src/sizefreq.size.gz -matfile ./src/matrices/single- -o ./output/output ./ 
```
* `-c`: coverage
* `-f`: file containing list of frequency of different fragment sizes
* `-matfile`: matrix file containing substitution misincorporation due to deamination
* `-o`: specify output directory and file name
As empirical evidence for frequency and deamination of aeDNA in glacial retreat context is lacking, precalculated values provided by gargammel based on chosen studies are used (`sizefreq.size.gz` and `src/matrices/`).  
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
$C = - \ln\left(1 - P^{L/G}\right)$
Where L is the average length of reads, G is genome assembly length, and P is the probability of genome covered entirely.  
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
Considering the storage limit of 1T for this project, risk of overflowing the storage space exists. Additionally, even with a high coverage, there is no guarantee that the final read set will cover the entire genome. Therefore, a deterministic tiling approach was adopted instead. 
As there is limited storage, an effort to split the assemblies into different batches for tiling was made. Assemblies selected for each batch were based on total_seq_length and test runs.  
This is then run with `slurm_tile.sbatch` using:
```
BATCH=batches/batch_XXX.list
N=$(wc -l < "$BATCH")
sbatch --array=1-"$N"%6 slurm_tile.sbatch "$BATCH"
```
where `batch_XXX.list` contains the list of assemblies directories for this particular batch.  
However, this greatly compromised downstream analysis as the newest release of Kraken2 core_nt database also requires 316.2G storage, and splitting into multiple batches would require much more rounds of processing.  
Extra research storage space was thus used, allowing for parallel processing of all 81 assemblies.  

**Reassigning Reads with Taxonomic Identification Algorithms**

When running Kraken2 with core_nt database on Cambridge HPC using SL-3 account, as core_nt database exceed the RAM limit, parameter `--memory-mapping` was used, which allows Kraken2 to search for matches without loading the entire database into RAM.  
However, this requires Kraken2 to frequently load and reload pages during processing, which greatly reduced its efficiency.  
As an effort to overcome this, running it in parallel on smaller files split from the original read files was attempted using `kraken.sbatch`.  
There were little to no improvements in processing efficiency. A node with higher RAM limit was thus borrowed to run without using `--memory-mapping`.  

When merging and sorting BAM files produced by bowtie2, the bowtie2 index was built from the full core_nt database. Each shard alignment BAM file therefore carries 1,562,196 @SQ header lines covering every reference sequence in the database. The failure mode differed depending on the samtools version used.  

The first attempted approach in `bamsort.sbatch` used the standard `samtools merge` followed by `samtools sort -n -O bam` pipeline, outputting a `.merged.sorted.bam` file. Under samtools 1.14, this failed immediately at the merge step with `[E::bam_hdr_write] Header too long for BAM format`, as version 1.14 enforces a strict internal limit on header size. Under samtools 1.19.2, this limit is relaxed, the merge step succeeded, but `samtools sort` ran out of memory for species with large read volumes (e.g. species 2996, 16G of aligned reads). This is because `samtools sort` pre-allocates its full sort buffer (10G per thread × 20 threads = 200G) before any data is processed, leaving only 40G of the 240G job allocation for the concurrent `samtools merge` process, OS overhead, and process memory. Species with smaller read volumes (e.g. species 1085736, 8G) fit within this headroom and succeeded, making the failure appear data-volume-dependent rather than systematic.  

The second attempted approach in `bamsort_sam.sbatch` avoided BAM output by extracting the header from one BAM with `samtools view -H`, streaming reads from each input BAM with `samtools view`, and piping the combined stream through `samtools sort -n -O SAM` before gzipping to `.merged.sorted.sam.gz`. This bypassed writing a BAM-format output file, but `samtools sort` still pre-allocates its full sort buffer and creates intermediate temporary BAM files internally during its merge-sort passes, making it subject to the same memory constraints.  

Thus, `bamsort_merge_sam.sbatch` was used to resolve this issue by replacing `samtools sort` with GNU `sort`, which operates on plain text without SAM header structure, creates no BAM temp files, and does not pre-allocate its sort buffer, but rather fills incrementally and spills to disk freely. Header issue is bypassed with `samtools view -H`, which writes the header directly to the output stream. This makes the approach robust across all species regardless of read volume or samtools version.  

**Reassigning Reads with K-mer Based Taxonomic Identification Algorithms - Kraken2**

To evaluate the best confidence threshold, 4 different values (0, 0.05, 0.1, 0.2) were tested on 10 sample species selected to represent the diversity of glacial retreat taxa composition.  The resulting kraken2 outputs were compared on their sensitivity (percentage of correctly assigned reads out of all reads) and precision (percentage of correctly assigned reads out of all classified reads) using f1 score calculated in `kraken_eval.py`. It accepts the following arguments:
* `--results_dir`: directory containing Kraken2 `.out` files
* `--nodes`: path to `nodes.dmp` from the Kraken2 database taxonomy
* `--out`: output CSV file path, default to `kraken_eval.csv`
The 10 sample species can be found in `sample_10_species.txt`, with evaluation results in `kraken_eval.csv`. For this specific project, the most ideal confidence threshold has been identified as 0.05 based on f1 score, but literature search for appropriate parameters has demonstrated 0.05 to 0.2 all having tests supporting their performance. Further consideration has deemed that to avoid including conserved tiles in the final candidate, a higher confidence threshold is necessary, and 0.2 is thus adopted.  

**References**

- Wood et al. (2019). Kraken2. https://doi.org/10.1186/s13059-019-1891-0  
- Wang et al. (2022). ngsLCA. https://doi.org/10.1111/2041-210X.14006
- Danecek et al. (2021). SAMtools. https://doi.org/10.1093/gigascience/giab008
- Lu et al. (2022). From defaults to databases: parameter and database choice dramatically impact the performance of metagenomic taxonomic classification tools. https://doi.org/10.1099/mgen.0.000949  
- Renaud et al. (2017). gargammel: a sequence simulator for ancient DNA. https://doi.org/10.1093/bioinformatics/btx512  