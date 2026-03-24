## Pipeline overview

This document contains step-by-step computational workflow used in this project, with relevant parameter choices.  

**Taxa List Synthesis**

List of taxa that may respond to glacial retreat is synthesized from past literature (see `species_id.txt`). Taxonomic ID at species level are extracted from NCBI.  

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
As genus metadata contains the metadata of representative species within the genus, species of interest need to be mapped to the corresponding representative species based on genus. This is done using `taxa_extraction.py` again. Output is shown in `species_rep_to_genus.tsv`.  

Using these following files:
* species genome metadata table
* genus genome metadata table
* species taxaID list
* species of interest to genus taxaID table
* species representative to genus taxaID table
`availability_check.py` try identify the best available assemblies for each species of interest. The script first identifies species level assemblies for species of interest. For those that do not have any assembly available, it falls back to genus level, then identifies representatives species assemblies within the genus. From all identified assemblies, the one with best quality is selected based on the following hierarchical priority:

| Priority | Database | Reference Genome | Assembly Level Order                             |
| -------- | -------- | ---------------- | ------------------------------------------------ |
| 1        | RefSeq   | Yes              | Complete genome > Chromosome > Scaffold > Contig |
| 2        | GenBank  | Yes              | Complete genome > Chromosome > Scaffold > Contig |
| 3        | RefSeq   | No               | Complete genome > Chromosome > Scaffold > Contig |

If none of the standard met or multiple assemblies with the same hierarchy appear for one species, choose the assembly with largest total sequence length. Those with no genome are left blank. Certain species have their subspecies/strains identified under taxonomic IDs different from the original ones, and are also included.  
Outputs are shown in `all_availability.tsv`.  

Genomes are downloaded according to assembly accession using `download_assemblies.sbatch`.  

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

**Reassigning Reads with K-mer Based Taxonomic Identification Algorithms - Kraken2**

An existing, wildly adopted taxonomic identification tools, Kraken2, is first used. 
Kraken2 core_nt database is chosen as the reference database for the purpose of this project. Prebuild core_nt database is acquired from Kraken 2 index zone. It is subsequently run using `sbatch_mom2.sh` with the following relevant parameters: 
``` 
kraken2 \
  --db cort_nt \
  --threads 48 \
  --confidence 0.2 \
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

To evaluate the best confidence threshold, 4 different values (0, 0.05, 0.1, 0.2) are tested on 10 sample species selected to represent the diversity of glacial retreat taxa composition.  The resulting kraken2 outputs are compared on their sensitivity (percentage of correctly assigned reads out of all reads) and precision (percentage of correctly assigned reads out of all classified reads) using f1 score calculated in `kraken_eval.py`. It accepts the following arguments:
* `--results_dir`: directory containing Kraken2 `.out` files
* `--nodes`: path to `nodes.dmp` from the Kraken2 database taxonomy
* `--out`: output CSV file path, default to `kraken_eval.csv`
The 10 sample species can be found in `sample_10_species.txt`, with evaluation results in `kraken_eval.csv`. For this speicfic project, the most ideal confidence threshold have been identified as 0.05 and is adopted for further processing.  

The Kraken2 outputs are then filtered using `kraken_filter.py`, which produces two output files per input:
1. `correct_genus.out`: reads classified at genus level or below where the assigned genus matches the true genus encoded in the read name
2. `genus_level.out`: reads classified at genus level or below regardless of whether the genus is correct
`kraken_filter.py` accepts the following arguments:
* `--results_dir`: directory containing Kraken2 `.out` files 
* `--nodes`: path to `nodes.dmp` from the Kraken2 database taxonomy 
* `--out_dir`: directory to write filtered output files 
Sample output can be found in `129212_task1.k2.0.05.core_nt.correct_genus.out` and `129212_task1.k2.0.05.core_nt.genus_level.out`.  

**Reassigning Reads with Local Alignment Based Taxonomic Identification Algorithms - Competitive Mapping**

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

After bowtie2 alignment, the output BAM files are merged and name-sorted with `bamsort.sbatch`. For species with larger read volumes, `bamsort.sbatch` fails due to hard header size limit or out of memory issue. `bamsort_merge_sam.sbatch` is used for these species instead, replacing `samtools sort` with GNU `sort`, which operates on plain text, does not pre-allocate a fixed buffer, and spills to disk freely, and avoid header issue by writing it directly to the output stream, passing only alignment lines through `sort`.  

For species with exceptionally large intermediate tmp files, a two-phase approach is further used. `bamsort_presort.sbatch` first namesorts each bowtie2 shard BAM independently into SAM.gz, `bamsort_merge_presorted.sbatch` then performs a streaming merge with `sort -m` and splits the merged output into multiple chunks as these typically produce large outputs and would risk hitting downstream wall time if output as a single file.  

For species that completed `bamsort_merge_sam.sbatch` successfully but hit downstream ngsLCA wall time limit, `bamsort_split.sbatch` is used to retroactively split the existing merged sorted SAM.gz into `*.chunkNN.sorted.sam.gz` files at read-name boundaries.  

As an alternative to the GNU sort approach, the BAM headers of all alignment outputs can first be compressed using `compressbam.sbatch`, which remove the header lines that are not referred to in the alignment. The compressed shards are then merged, queryname-sorted, and split into 4 chunk BAMs using `bamsort_compressed.sbatch`.  

ngsLCA is then run on the sorted files with `ngslca.sbatch`, which prioritises chunk files over any full-length file for the same species. Because the bowtie2 reference database was converted from the Kraken2 core_nt database, taxonomy files from Kraken2 are used directly: `nodes.dmp`, `names.dmp`, and `seqid2taxid.map`. `seqid2taxid.map` uses a Kraken2-specific key format (`kraken:taxid|TAXID|ACCESSION`) that does not match the bare accession names in the SAM RNAME field. It is therefore converted to a 4-column NCBI acc2tax format by stripping the prefix:
```
awk 'BEGIN{OFS="\t"}{n=split($1,a,"|"); acc=a[n]; print acc, acc, $2, 0}' \
    seqid2taxid.map > seqid2taxid.acc2tax
```
An example ngsLCA command with relevant parameters for a single species is:
```
ngsLCA \
    -simscorelow 0.95 \
    -simscorehigh 1.0 \
    -fix-ncbi 0 \
    -names names.dmp \
    -nodes nodes.dmp \
    -acc2tax seqid2taxid.acc2tax \
    -bam species_id.merged.sorted.bam \
    -outnames species_id
```
* `-simscorelow`: minimum alignment similarity score to consider a read
* `-simscorehigh`: maximum alignment similarity score to consider a read
* `-fix-ncbi`: whether to apply NCBI-specific accession fixes (0 = off)
* `-names`: Kraken2 taxonomy names file
* `-nodes`: Kraken2 taxonomy nodes file
* `-acc2tax`: accession-to-taxid mapping file
* `-bam`: input name-sorted BAM or SAM file
* `-outnames`: prefix for output files

Each chunk file is processed as a separate SLURM array task, producing a corresponding `species_id_chunkNN.lca` file. For species that completed ngsLCA as a single full file but whose SAM was later split into chunks with `bamsort_split.sbatch`, `lca_split.sbatch` can retroactively split the full LCA to match the SAM chunk boundaries by scanning chunk SAM to find its first read name, building a boundary list for dividing LCA sequentially. This ensures the split is robust even when a read is absent from the LCA due to ngsLCA filtering.  

Bamdam is then run on the ngsLCA outputs with `bamdam.sbatch` to filter the BAM and LCA files down to reads assigned at genus level or below. When both a full LCA and chunk LCAs exist for the same species. only the chunk LCAs are processed. An example command with relevant parameters is:
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
For species with sorted file in `.sam.gz` format, the script first converts it to BAM with `samtools view`, as bamdam requires BAM input. For chunk items, bamdam produces per-chunk shrunk outputs. Once all chunks for a species are complete, `bamdam_merge.sbatch` merges the chunk BAMs with `samtools merge` and concatenates the chunk LCA files to produce a single output per species.  

The shrunk LCA files are then further filtered using `bamdam_filter.py`, which produces one output file per species:
1. `correct_genus.lca`: reads assigned at genus level or below where the assigned genus matches the true genus encoded in the read name
`bamdam_filter.py` accepts the following arguments:
* `--lca_dir`: directory containing bamdam `.shrunk.lca` files
* `--nodes`: path to `nodes.dmp` from the NCBI taxonomy
* `--out_dir`: directory to write filtered output files

**Database Coverage Check**

As the presence or absence of targeted or representative species assmebly in the database may affect the performance of both pipelines, to assess their corresponding performance under different level of database coverage, species that have any genetic assembly in the database is identified. This is done using `check_db_coverage.py`, which first checks for the list of representative species taxid by referring to `all_availability.tsv` and `genus_metadata.tsv`: 
```
python3 check_db_coverage.py \
  --all_avail data/all_availability.tsv \
  --genus_meta data/genus_metadata.tsv \
  --out_dir data/
```
It is important to point out that this only checks for the presence of any genetic assembly of these targeted or representative species, as search for the exact assembly used for tiling have been done and revealed that no species have their exact assembly accession representation in the database.  

Outputs `organism_taxid_mapping.tsv` and `organism_taxids_to_check.txt` are then compared to ngsLCA and Kraken2 accession-taxid metadata: 
```
awk 'NR==FNR{t[$1]=1;next} $3 in t {print $3}' \
  organism_taxids_to_check.txt \
  nucl_gb.accession2taxid \
  | sort -u > in_ngslca.txt

awk '$1>0 {print $5}' inspect.txt | sort -u > kraken2_all_taxids.txt
grep -Fwf organism_taxids_to_check.txt kraken2_all_taxids.txt | sort -u > in_kraken2.txt
```
`check_db_coverage.py` are then run again with the ngslca and kraken2 comparison results: 
```
python3 check_db_coverage.py \
  --all_avail data/all_availability.tsv \
  --genus_meta data/genus_metadata.tsv \
  --out_dir data/ \
  --in_ngslca data/in_ngslca.txt \
  --in_kraken data/in_kraken2.txt
```
Outputs `db_coverage.tsv` with per-species database presence flags, and `species_in_both_dbs.txt` listing target taxids found in both databases. Only species in `species_in_both_dbs.txt` are used for downstream pipeline comparison.  

**Extending the Local Alignment Pipeline with New Assemblies**

A key advantage of the local alignment pipeline over the k-mer based pipeline is its modularity. While Kraken2's k-mer index must be reconstructed from scratch whenever new sequences are added, with local alignment pipeline, new assemblies can be appended to the core_nt database and indexed as additional bowtie2 shards, allowing reads to be mapped against the extended database incrementally. The local alignment database is thus extended to include the assemblies of target species and representative speices that were absent from the original core_nt database.  

New assemblies are added to core_nt as shards 78, 79, and 80. The bowtie2 index for these new shards is built using `bowtie2_index_shard78.sbatch`. Reads are then aligned against these shards using `bowtie2_run_shard78.sbatch`. The compressed shard BAMs are then merged, queryname-sorted, and split into chunks with the alignments of previously existed shards. The resulting chunk BAMs feed into the same ngsLCA → bamdam → bamdam_filter pipeline, with the taxonomy files updated to include the new assemblies using `make_seqid2taxid_shard78.py`. This is done using corresponding `_compressed` variants of each script (`ngslca_compressed.sbatch`, `bamdam_compressed.sbatch`, `bamdam_filter_compressed.py`).  

**Modern DNA Comparison**

To enable a direct comparison between simulated damaged aeDNA reads and undamaged modern reads, the same kraken2 → kraken-filter pipeline and bowtie2 → bamsort → ngsLCA → bamdam → bamdam-filter pipeline are applied to modern reads for the 7 species found in both databases. The corresponding `_modern` variants of each script (`bowtie2_modern.sbatch`, `bowtie2_run_shard78_modern.sbatch`, `compressbam_modern.sbatch`, `bamsort_modern.sbatch`, `bamsort_compressed_modern.sbatch`, `ngslca_modern.sbatch`, `ngslca_compressed_modern.sbatch`, `bamdam_modern.sbatch`, `bamdam_compressed_modern.sbatch`, `bamdam_filter_compressed_modern.sbatch`) are used with the same parameters as the aeDNA steps. Modern reads are generated again with `aeDNA_simulation.py`:
```
python3 aeDNA_simulation.py \
  --assembly assembly_accession_genomic.fna \
  --metadata all_availability2.tsv \
```

**Biophysical Filter**

To prepare eprobe input, the overlap of reads correctly assigned at genus level by both k-mer and the local alignment pipeline is computed for each species, then extracted from the undamaged modern tiled FASTAs using `eprobe_input_extract.py`. The script accepts the following arguments:
* `--species_id`: species taxid to process
* `--kraken_filt_dir`: directory containing Kraken2 `correct_genus.out` files
* `--bamdam_filt_dir`: directory containing bamdam `correct_genus.lca` files
* `--modern_tiled_dir`: directory containing modern tiled FASTA files
* `--out_dir`: output directory for eprobe input FASTAs
Using its sbatch wrapper `eprobe_input_extract.sbatch`, output files are named as `{sid}.fasta`.  

For species where the overlap is too low to yield sufficient input reads, `eprobe_input_bamdam_only.py` is used instead, extracting sequences based on bamdam-filtered read names only without requiring Kraken2 agreement. These species are run with `eprobe_input_bamdam_only.sbatch` and output as `{sid}_bamdam.fasta`.  

The eProbe biophysical filter is then applied to each input FASTA using `eprobe_filter.sbatch`, which runs as a SLURM array with one task per species. An example command with relevant parameters for a single task is:
```
eprobe util filter \
    --input  species_id.fasta \
    --output species_id \
    --step   biophysical \
    --threads 4
```
Output filtered FASTAs are written as `{species_id}.filtered.fa`.  

For species with large input FASTAs that cause memory errors in the eprobe dimer filter, `subsample_fasta.py` is used to randomly subsample to 300,000 reads before proceeding, producing `{sid}_sub300k.fasta`. It accepts the following arguments:
* `--input`: input FASTA file
* `--output`: output FASTA file
* `--n`: number of reads to sample
* `--seed`: random seed, optional

Filter results per species are summarized using `eprobe_summary.py`, which counts input, passed, and rejected sequences for each species and records the input type (plain, bamdam, sub300k, or bamdam_sub300k). It accepts the following arguments:
* `--eprobe_input_dir`: eprobe_input directory
* `--eprobe_filtered_dir`: eprobe_filtered directory
* `--out`: output CSV file path
Output is shown in `eprobe_summary.csv`.   


**Shuffle, Subsampling and Deduplication**

Before deduplication, each species' eprobe-filtered output is subsampled according to a pre-computed plan in `subsample_plan.csv` using `subsample_for_cdhit.py`. `subsample_plan.csv` contains three columns per species: `species_id`, `eprobe_passed`, and `ratio` (eprobe_passed / 625). The target number of reads is fixed at 625 for all species. The script compares each species' ratio against a set threshold; species with ratio above the threshold are subsampled down to 625 × threshold reads, while species at or below the threshold retain all available reads. The script accepts the following arguments:
* `--plan`: subsample plan CSV
* `--filtered_dir`: eprobe_filtered directory
* `--out_dir`: output directory for subsampled FASTAs
* `--threshold`: ratio threshold above which subsampling is applied, default 5.0; 4.0 is used for this project
* `--seed`: random seed for reproducibility (optional)
Output FASTAs are written to `deduplicate_input/` as `{species_id}.fasta`.  

All subsampled FASTAs are then pooled and deduplicated across species using `cdhit.sbatch`, which concatenates all inputs into a single `pooled.fasta` and runs CD-HIT-EST with the following parameters:
```
cd-hit-est \
    -i  pooled.fasta \
    -o  pooled_dedup \
    -c  0.95 \
    -n  8 \
    -aS 1.0 \
    -r  1 \
    -d  0
```
* `-c`: sequence identity threshold
* `-n`: word length
* `-aS`: minimum alignment coverage of the shorter sequence
* `-r`: align both strands
* `-d`: full FASTA header in cluster file
Outputs are `pooled_dedup`, which contains one representative sequence per cluster, and `pooled_dedup.fasta.clstr`, which shows the cluster assignments. The final probe panel contains 51220 reads.  

**Summary, Analysis and Visualization**

Off-target alignment analysis is performed using `primary_secondary_analysis.py`, which reads the merged bowtie2 BAM for each species, computes alignment similarity using `(aligned_length − NM) / aligned_length` for each alignment, collapses hits to the best similarity per taxonomic group at genus level, excludes the target genus, and reports what fraction of reads have off-target hits above the similarity threshold, and the top off-target taxa. It accepts the following arguments:
* `--sam`: input SAM / SAM.gz / BAM file or stdin
* `--target_taxid`: target species taxid whose genus ancestor will be excluded
* `--acc2tax`: accession-to-taxid mapping file
* `--nodes`: path to `nodes.dmp`
* `--names`: path to `names.dmp`
* `--out_summary`: output summary text file
* `--level`: taxonomic level to collapse hits to, default genus
* `--simscorelow`: similarity threshold, default 0.95
* `--top`: number of reported top off-target taxa, default 20
This is ran via `primary_secondary_lowsignal.sbatch` for 20 low-signal simulated ancient species, 1 low-signal simulated modern species, and 1 high-signal simulated ancient species as comparison.  

Read counts at each pipeline stage per species are acquired using `count_stats.py` and `count_stats_compressed.py`. These scripts collect counts across the tiled reads, ngsLCA output, bamdam output, bamdam correct_genus, Kraken2 total, Kraken2 genus level, Kraken2 correct_genus, and eprobe input overlap columns. `count_stats_compressed.py` handles chunk-based outputs from the compressed pipeline by summing across chunk files. Output can be seen at `pipeline_stats.csv` and `pipeline_stats_compressed.csv` respectively.  

Pipeline results for the 7 species found in both databases are visualized using `pipeline_visualize.py` and `pipeline_visualize_compressed.py`. Both scripts produce three figures comparing the k-mer based and local alignment based pipelines for simulated aeDNA-damaged and modern reads side by side:
1. `kraken_summary.png`: grouped bar chart of read counts at each k-mer stage per species
2. `competitive_summary.png`: grouped bar chart of read counts at each local alignment stage per species
3. `intersection.png`: stacked bar chart of correct-genus read overlap between k-mer and local alignment per species
A `summary_table.csv` is also written with all counts and overlaps for reference. `pipeline_visualize_compressed.py` uses the compressed pipeline outputs.  

Read count comparisons across pipeline stages for the 7 species are produced by `pipeline_compare.py` and `pipeline_compare_compressed.py`, outputting one row per (species, mode) combination to `pipeline_comparison.csv`. All directory arguments are optional, with missing files reported as empty cells.  

`taxonomy_availability_plot.py` produces a stacked bar chart grouping retreat species by taxonomic group, default at phylum, coloured by assembly availability (species-level assembly, genus-level assembly, or none). Higher-level taxonomy is fetched from NCBI Entrez with results cached locally for subsequent runs. It accepts the following arguments:
* `--input`: path to `all_availability.tsv`
* `--out`: output directory for the figure
* `--email`: email address for NCBI Entrez
* `--rank`: taxonomic rank for grouping, default phylum
* `--cache`: TSV cache file for NCBI taxonomy lookups, default `data/taxonomy_cache.tsv`
Output figure can be seen at `taxonomy_availability_{rank}.png`.  

`probe_funnel_plot.py` produces a Sankey-style funnel plot showing read count changes across the four probe design stages: eprobe input, epobre biophysical filter, shuffle with subsampling, and CD-HIT deduplication. Each bar represents 100% of reads entering that stage, with the coloured portion showing what fraction is retained and the grey portion showing what is dropped. Percentages of read kept are computed relative to the previous stage, with absolute read counts shown above each bar. It accepts the following arguments:
* `--summary`: path to `eprobe_summary.csv`
* `--plan`: path to `subsample_plan.csv`
* `--post_dedup`: read count after CD-HIT deduplication (integer)
* `--threshold`: subsample ratio threshold used, default 5.0; 4.0 is used for this project
* `--out`: output figure path, default `figures/probe_funnel.png`
Output figure can be seen at `probe_funnel.png`.  

`pooled_dedup_piechart.py` produces two pie charts for comparing species composition by taxonomic rank in the final probe set and the input species distribution, using the same NCBI taxonomy cache as `taxonomy_availability_plot.py`. Both charts share an identical color scheme so the same phylum receives the same color across both figures. It accepts the following arguments:
* `--pooled_dedup`: path to the `pooled_dedup` FASTA
* `--summary`: path to `eprobe_summary.csv`
* `--cache`: TSV taxonomy cache file, default `data/taxonomy_cache.tsv`
* `--email`: email for NCBI Entrez, required only if taxids are missing from cache
* `--rank`: taxonomic rank for grouping, default phylum
* `--out_reads`: output path for reads pie chart, default `figures/pooled_dedup_piechart.png`
* `--out_species`: output path for species-count pie chart, default `figures/species_piechart.png`
Output figures can be seen at `pooled_dedup_piechart.png` and `species_piechart.png`.

`probe_count_lollipop.py` produces a Manhattan-style lollipop plot showing the number of probes in `pooled_dedup` per species, grouped and colour-coded by taxonomic rank. It also uses the same NCBI taxonomy cache and `names.dmp` for species name lookup. It accepts the following arguments:
* `--pooled_dedup`: path to the `pooled_dedup` FASTA from CD-HIT
* `--names`: path to NCBI `names.dmp` for species name lookup
* `--cache`: TSV taxonomy cache file, default `data/taxonomy_cache.tsv`
* `--email`: email for NCBI Entrez, required only if taxids are missing from cache
* `--rank`: taxonomic rank for grouping and colouring, default phylum
* `--out`: output figure path, default `figures/probe_count_lollipop.png`
Output figure can be seen at `probe_count_lollipop.png`.

`pooled_dedup_counts.py` counts the number of probes per species in the `pooled_dedup` FASTA, attaches species names and phylum from the taxonomy cache, and outputs a sorted TSV table. It accepts the following arguments:
* `--pooled_dedup`: path to the `pooled_dedup` FASTA from CD-HIT
* `--names`: path to NCBI `names.dmp` for species name lookup
* `--cache`: TSV taxonomy cache file, default `data/taxonomy_cache.tsv`
* `--out`: output TSV path, default `data/pooled_dedup_counts.tsv`
Output is written to `pooled_dedup_counts.tsv` and also printed to stdout.

`correct_rate_correlation.py` produces a scatter plot comparing the local alignment correct rate against the k-mer correct rate per species. It accepts the following arguments:
* `--input`: path to `pipeline_stats_compressed.csv`, default `data/pipeline_stats_compressed.csv`
* `--out`: output figure path, default `figures/correct_rate_correlation_compressed.png`
Output figure can be seen at `correct_rate_correlation_compressed.png`.  

**Declaration of AI Usage**

All sbatch scripts are written with Claude Code Sonnet 4.6, then proofread and annotated by me, except `sbatch_mom2.sh`, which is provided by Prof. Anton Enright, CGS, Department of Pathology, University of Cambridge.  

Python scripts written and annotated by me, debug with ChatGPT or Claude Code Sonnet 4.6:
```
taxa_extraction.py
availability_check.py
aeDNA_simulation.py
kraken_filter.py
bamdam_filter_compressed.py (based on bamdam_filter.py)
eprobe_input_bamdam_only.py (based on eprobe_input_extract.py)
eprobe_summary.py
subsample_fasta.py
subsample_for_cdhit.py
count_stats_compressed.py (based on count_stats.py)
pipeline_visualize_compressed.py (based on pipeline_visualize.py)
pipeline_compare_compressed.py (based on pipeline_compare.py)
```
The rest of the python scripts are written with Claude Code Sonnet 4.6, then proofread and annotated by me. 