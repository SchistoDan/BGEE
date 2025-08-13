# Barcode Gene Extractor & Evaluator (BGEE) Snakemake workflow #
Snakemake workflow for recovering high-quality barcode sequences from genome skim data, built around MitoGeneExtractor and adapted for genome skims of museum specimens.

# Contents # 
 - [Requirements](#Requirements)
 - [Workflow](#Workflow)
 - [Running](#running)
 - [Cluster configuration](#Cluster-configuration-using-profiles/slurm/config.yaml)
 - [Cluster submission](#Cluster-submission)
 - [Results structure](#Results-structure)
 - [Contributing](#Contributing)

# Requirements #
- [MitoGeneExtractor](https://github.com/cmayer/MitoGeneExtractor) version 1.9.6 installed. Clone repository and follow installation instructons. 
- Paired-end reads in .fastq.gz or .fastq format.
- samples_file.csv (generated manually, or as outlined below if working from BOLD sample metadata).
- sequence_references_file.csv (generated manually, or using [Gene Fetch](https://github.com/bge-barcoding/gene_fetch?tab=readme-ov-file) within the workflow).
- Activated conda env (see bgee_env.yaml).

# Workflow #
1. Preprocessing mode (both pre-processing modes are run in parallel):
  - 'concat':
    - Initial raw read quality control - Adapter, trimming, quality trimming, poly-g trimming, and deduplication of paired-end reads using [fastp](https://github.com/OpenGene/fastp) (fastp_pe_concat).
    - Concantenation of trimmed PE reads (fastq_concat) and associated log files (aggregate_concat_logs).
    - Compress trimmed reads (gziped_trimmed).
    - Secondary read quality control - Additional quality trimming of concatenated reads using [trim galore](https://github.com/FelixKrueger/TrimGalore) (quality_trim) and concatenation of associated log files (Aggregate_trim_galore_logs).
  - 'merge:
    - Raw read quality control and merging - Adapter, trimming, quality trimming, poly-g trimming, deduplication, and merging of paired-end reads using [fastp](https://github.com/OpenGene/fastp) (fastp_pe_merge).
    - 'Clean' headers of input files, as required by MitoGeneExtractor (clean_headers_merge), and concentation of associated log files (Aggregate_clean_headers_logs).

<p align="center">
  <img src="https://github.com/user-attachments/assets/139b8c7c-b0dc-465c-8c95-e3a58ea1ab96" width="500"/>
</p>

2. Retrieval of sample-specific pseudo-references from GenBank using [Gene-Fetch](https://github.com/bge-barcoding/gene_fetch). (gene_fetch).
3. Protein reference-guided barcode recovery using [MitoGeneExtractor](https://github.com/cmayer/MitoGeneExtractor) (MitoGeneExtractor_concat & MitoGeneExtractor_merge).
4. Rename raw consensus sequence headers and concatenation into multi-FASTA (rename_and_combine_cons_concat & rename_and_combine_cons_merge) (uses supplementary [rename_headers.py](https://github.com/bge-barcoding/MitoGeneExtractor-BGE/blob/main/workflow/scripts/rename_headers.py).
5. Remove any remaining exonerate intermediate (Concatenated_exonerate_input_*) files leftover after MGE (remove_exonerate_intermediates).
6. Create list of MGE (FASTA) alignment files for downstream processing (create_alignment_log).
7. Filter MGE alignment files to remove low-quality, contaminant, or outlier sequences before repeat consensus sequence generation.
   - Remove aligned reads with similarity to human COI (a common contaminant of museum specimens) (uses supplementary [human_cox1_filter.py](https://github.com/bge-barcoding/MitoGeneExtractor-BGE/blob/main/workflow/scripts/human_cox1_filter.py).
   - Remove aligned reads with AT content over and/or below a specified threshold (high AT content can be indicative of contamination (e.g. fungal)) (uses supplementary [at_content_filter.py](https://github.com/bge-barcoding/MitoGeneExtractor-BGE/blob/main/workflow/scripts/at_content_filter.py).
   - Remove reads that are statistical outliers compared to the original consensus (uses supplementary [statistical_outliers.py](https://github.com/bge-barcoding/MitoGeneExtractor-BGE/blob/main/workflow/scripts/statistical_outlier_filter.py).
   - [optional] Remove reads with similarity to supplied reference sequence(s) (uses supplementary [reference_filter.py](https://github.com/bge-barcoding/MitoGeneExtractor-BGE/blob/main/workflow/scripts/reference_filter.py).
   - Generation of 'cleaned' consensus sequence (uses supplementary [consensus_generator.py](https://github.com/bge-barcoding/MitoGeneExtractor-BGE/blob/main/workflow/scripts/consensus_generator.py).
   - Aggregate metrics from each stage of filtering (uses supplementary [aggregate_filter_metrics.py](https://github.com/bge-barcoding/MitoGeneExtractor-BGE/blob/main/workflow/scripts/aggregate_filter_metrics.py).
   - Remove intermediate files and unecessary logs generated during the consensus cleaning process (remove_fasta_cleaner_files).
8. Evaluate barcode consensus sequence quality based on various metrics (length, ambiguous base content, etc.), and select the 'best' sequences according to specific ranking criteria (either 'normal' or 'relaxed'). (see docstring of supplementary [fasta_compare](https://github.com/bge-barcoding/MitoGeneExtractor-BGE/blob/main/workflow/scripts/fasta_compare.py).
9. Compile statistics from read QC, MGE, and consensus cleaning metrics into a CSV report for both 'concat' and 'merge' modes (uses supplementary [mge_stats.py](https://github.com/bge-barcoding/MitoGeneExtractor-BGE/blob/main/workflow/scripts/mge_stats.py) and combine_stats_files).
10. Remove temporary files, sample-specific logs once aggregated, etc. (cleanup_files).



# Running: #
## Clone BGEE github repository and set up conda environment ##
- [Install miniconda](https://www.anaconda.com/docs/getting-started/miniconda/install#quickstart-install-instructions).
```bash
git clone https://github.com/bge-barcoding/BGEE.git [path/to/desired/install/location/]
cd installation/dir/BGEE
conda env create -f bgee_env.yaml
git status
```

## Generate sample input file ###
- Can be created manually, or via [sample-processing](https://github.com/bge-barcoding/sample-processing) workflow.
- Must contain `ID`, `forward` (read paths), `reverse` (read paths), and `taxid` **or** `hierarchical taxonomy` columns (see below for examples).
- Column 1 can be named 'ID', 'process_id', 'Process ID', 'process id', 'Process id', 'PROCESS ID', 'sample', 'SAMPLE', or 'Sample'.
- Due to regex matching and statistics aggregation, the sample ID will be considered as the string before the first underscore. It is therefore recommended that sample names do not use '_' characters. E.g. BSNHM002-24 instead of BSNHM002_24, or P3-1-A10-2-G1 instead of P3_1_A10_2_G1.
- Taxid's can be found manually by searching the expected species/genus/family of each sample in the [NCBI taxonomy database](https://www.ncbi.nlm.nih.gov/taxonomy).
- If hierarchical taxonomic information for each sample is provided, `Gene Fetch` will, starting from species, find the closest valid taxid on NCBI for the provided taxonomy before proceeding with pseudo-reference fetching.
  
**samples.csv example (taxids)**
| ID | forward | reverse | taxid |
| --- | --- | --- | --- |
| BSNHM002-24  | abs/path/to/R1.fq.gz | abs/path/to/R2.fq.gz | 177658 |
| BSNHM038-24 | abs/path/to/R1.fq.gz | abs/path/to/R2.fq.gz | 177627 |
| BSNHM046-24 | abs/path/to/R1.fq.gz | abs/path/to/R2.fq.gz | 3084599 |

**samples.csv example (hierarchical taxonomy)**
| ID | forward | reverse | phylum | class | order | family | genus | species |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | 
| BSNHM002-24  | abs/path/to/R1.fq.gz | abs/path/to/R2.fq.gz | phylum | class | order | family | genus | species |
| BSNHM038-24 | abs/path/to/R1.fq.gz | abs/path/to/R2.fq.gz | phylum | class | order | family | genus | species |
| BSNHM046-24 | abs/path/to/R1.fq.gz | abs/path/to/R2.fq.gz | phylum | class | order | family | genus | species |

## Gathering sample-specific pseudo-references ##
- This can be created manually, or using [Gene-fetch](https://github.com/bge-barcoding/gene_fetch) integrated into the workflow. If enabled (in the config.yaml by setting `run_gene_fetch` to 'true'), gene-fetch will retrieve the necessary protein pseudo-references for each sample from NCBI GenBank using the samples taxonomic identifier (taxid)/taxonomic lineages for each sample, a sequence target (e.g. COI or rbcL), and NCBI API credentials (email address & API key - see [guidance](https://support.nlm.nih.gov/kbArticle/?pn=KA-05317) on getting a key). 
- Must contain 'process_id', 'reference'name' and 'protein_reference_path' at a minimum.

**sample_references.csv example**
| process_id | reference_name | protein_reference_path | 
| --- | --- | --- |
| BSNHM002-24  | BSNHM002-24 | path/to/BSNHM002-24.fasta |
| BSNHM038-24 | BSNHM038-24 | path/to/BSNHM038-24.fasta |
| BSNHM046-24 | BSNHM046-24 | path/toBSNHM046-24.fasta |

## Customising snakemake configuration file ##
- Update `config/config.yaml` with neccessary paths and variables - See [MitoGeneExtractor README.md](https://github.com/bge-barcoding/MitoGeneExtractor-BGE/blob/main/README.md) for a more detailed explanation of Exonernate run paramters.
- Rule-specific resources in the [config.yaml](https://github.com/bge-barcoding/MitoGeneExtractor-BGE/blob/main/config/config.yaml) - Each rule can specify the necessary number of threads and memory resources (in Mb) for every job (e.g. specifying 4 threads and 4G memory for fastp_pe_merge would allocate those resources for every 'fastp_pe_merge' job).
```
# config/config.yaml

## General BGEE pipeline parameters and paths
# BGEE run name identifier
run_name: "BGEE_RUN_1"
# Path to MGE installation (MitoGeneExtractor-vX.X.X file)
mge_path: "path/to/MitoGeneExtractor-1.9.6beta3/MitoGeneExtractor-v1.9.6beta3"
# Path to samples.csv
samples_file: "path/to/samples.csv"
# Path to references.csv (comment out 'sequence_reference_file' line below or remove path if 'run_gene_fetch' is true)
sequence_reference_file: ""
# Path to output directory. Will make final dir if it does not exist already
output_dir: "path/to/BGEE-RUN_1-out"

## Gene Fetch parameters (https://github.com/bge-barcoding/gene_fetch)
run_gene_fetch: true  # Set true to use gene-fetch to generate reference sequences
gene_fetch:
  email: "example@example.ac.uk"      # Required: Email for NCBI API
  api_key: "api_key1234567890"        # Required: NCBI API key
  gene: "cox1"                        # Target gene (e.g., cox1, rbcl, matk)
  input_type: "taxid"                                     # Does the 'samples.csv' contain a 'taxid' column or 'hierarchical' taxonomic information column per sample (default: taxid)? - see Gene Fetch docs for more detail
  genbank: true                       # Download GenBank records for corresponding pseudo-references
  # output_dir will be set automatically to {output_dir}/references/

## MGE running parameters (see https://github.com/cmayer/MitoGeneExtractor/tree/main?tab=readme-ov-file#command-line-options)
# Exonerate relative score threshold parameter
r:
  - 1
  - 1.3
  - 1.5
# Exonerate minimum score threshold parameter
s:
  - 50
  - 100
# Number of bp to extend beyond the Exonerate alignment
n: 0
# Genetic code (see https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi)
C: 5
# Consensus threshold (e.g. 0.5 = 50%)
t: 0.5
  
## Post-processing of aligned reads for cleaning (using five supplementary scripts)
# human coi filtering -> at content filtering -> statistical outlier filtering -> (optional) reference-base filtering -> 'cleaned' consensus generation
# Default parameters
fasta_cleaner:
  consensus_threshold: 0.5      # Threshold at which bases at each position must 'agree' to be incldued in the consensus (0.5 = ≥50% of bases at each position must agree)
  human_threshold: 0.95         # Threshold at which reads are removed due to similarity with human COI (0.95 = reads with ≥95% similarity are removed)
  at_difference: 0.1            # Threshold at which reads are removed due to AT content variation (0.1 = reads with AT% differing by >10% from the consensus are removed)
  at_mode: "absolute"           # At content filtering mode: Absolute (reads removed when AT% differs from consensus in either direction by `at_difference` (e.g. >10%)), Higher (reads removed when AT% is higher than `at_difference` threshold), Lower (reads removed when AT% is lower than `at_difference` threshold).
  outlier_percentile: 90.0      # Threshold at which reads are flagged as statistical outliers to the consensus and removed (90.0 = reads <90% 'similar' to the consensus are removed)
  disable_human: false          # If false, 'human_cox1_filter.py' will run
  disable_at: false             # If false, 'at_content_filter.py' will run
  disable_outliers: false       # If false, 'statistical_outliers.py' will run
  reference_dir: null           # To skip reference filtering enter "null"/"None". To enable reference filtering, provide path to a parent directory (flat structure) containing '{sample}_reference.fasta' files.

## Comparison and selection of the 'best' consensus sequence produced by MGE and fasta_cleaner for each sample (using fasta_compare.py)
run_fasta_compare: true  # Set false to skip fasta_compare step
fasta_compare:
  target: "cox1"    # Options: cox1, rbcl, matk
  verbose: false    # Enable verbose logging (default: false)
  rank: 3           # Barcode rank and below considered a 'pass' (default: 3)
  relaxed: false    # Use 'relaxed' ranking system (default: false)

## Resource allocation for each rule. Memory shown in Mb (mem_mb), e.g. 2048 = 2G memory.
# Rules have dynamic memory scaling upon retry (mem_mb * retry #).
rules:
  gene_fetch:
    mem_mb: 4096
    threads: 2
  fastp_pe_merge:
    mem_mb: 8192
    threads: 4
  clean_headers_merge:
    mem_mb: 4096
    threads: 2
  aggregate_clean_headers_logs_merge:
    mem_mb: 2048
    threads: 2
  fastp_pe_concat:
    mem_mb: 8192
    threads: 4  
  fastq_concat:
    mem_mb: 4096
    threads: 4
  gzip_trimmed:
    mem_mb: 4096
    threads: 2 
  aggregate_concat_logs:
    mem_mb: 2048
    threads: 2
  quality_trim:
    mem_mb: 8192
    threads: 4
  aggregate_trim_galore_logs:
    mem_mb: 2048
    threads: 2
  MitoGeneExtractor_merge:
    mem_mb: 20480
    threads: 4
    partition: medium        # Change according to your available partitions
  MitoGeneExtractor_concat:
    mem_mb: 20480
    threads: 4
    partition: medium        # Change according to your available partitions
  rename_and_combine_cons_merge:
    mem_mb: 2048
    threads: 2
  rename_and_combine_cons_concat:
    mem_mb: 2048
    threads: 2
  create_alignment_log_merge:
    mem_mb: 4096
    threads: 2
  create_alignment_log_concat:
    mem_mb: 4096
    threads: 2
  human_cox1_filter_merge:
    mem_mb: 10240
    threads: 4
  at_content_filter_merge:
    mem_mb: 15360
    threads: 8
  statistical_outlier_filter_merge:
    mem_mb: 8192
    threads: 4
  reference_filter_merge:
    mem_mb: 8192
    threads: 4
  consensus_generation_merge:
    mem_mb: 8192
    threads: 4
  aggregate_filter_metrics_merge:
    mem_mb: 2048
    threads: 1
  human_cox1_filter_concat:
    mem_mb: 10240
    threads: 4
  at_content_filter_concat:
    mem_mb: 15360
    threads: 8
  statistical_outlier_filter_concat:
    mem_mb: 4096
    threads: 4
  reference_filter_concat:
    mem_mb: 8192
    threads: 4
  consensus_generation_concat:
    mem_mb: 8192
    threads: 4
  aggregate_filter_metrics_concat:
    mem_mb: 2048
    threads: 1
  extract_stats_to_csv_merge:
    mem_mb: 4096
    threads: 4
  extract_stats_to_csv_concat:
    mem_mb: 4096
    threads: 4
  combine_stats_files:
    mem_mb: 1024
    threads: 1
  fasta_compare:
    mem_mb: 4096
    threads: 4
  remove_fasta_cleaner_files:
    mem_mb: 4096
    threads: 4
  remove_exonerate_intermediates:
    mem_mb: 4096
    threads: 2
  cleanup_files:
    mem_mb: 1024
    threads: 1
```

# Cluster configuration using Snakemake profiles #
- See `profiles/` directory for 'slurm' and 'local' (i.e. non-SLURM) cluster submission parameters.
- The profile (`profiles/local` or `profiles/slurm`) will need to be changed depending on your system (see `$PROFILE` variable in `snakemake_run.sh`).

# Cluster submission #
- [snakemake_run.sh](https://github.com/bge-barcoding/MitoGeneExtractor-BGE/blob/main/snakemake_run.sh) handles submission of the snakemake workflow to the HPC cluster. The working directory will initially be unlocked (using `--unlock`) and then the snakemake workflow will be run.  

If using `profiles/slurm`, SLURM will orchestrate submission of each step in the workflow as a separate job.
1. Set up interactive session on compute node using `srun --job-name=snakeflow --partition=[YOUR_PARTITION] --cpus-per-task=2 --mem=4G --time=24:00:00 --pty bash`. Change `[YOUR_PARTITION]` to a suitable partition for your HPC. Remember choose a partition and request enough walltime (with `--time`) for the scheduler to last until all jobs finish, otherwise it will be killed when interactive session ends.
2. Submit snakemake_run.sh to the compute node running the interactive session you just set up with `bash snakemake_run.sh` (when in `BGEE/` directory). You might see the warning "You are running snakemake in a SLURM job context. This is not recommended, as it may lead to unexpected behavior. Please run Snakemake directly on the login node." - this should be fine as the interactive (SLURM) session is merely handling workflow job submission. If you do experience problems, try `./snakemake_run.sh` within a screen session on your cluster's head/login node.  

If using `profiles/local`, all workflow steps will be run as a single job. 
1. Simply run `./snakemake_run.sh` on your desired cluster compute node.



# Results structure #
```
output_dir/
├── references/
│   ├── protein/                                               # Contains FASTA pseudo-references for each sample.
│   ├── genbank/                                               # Contains GenBank records corresponding to fetched pseudo-references (if 'genbank: true' set in config.yaml).
│   └── sequence_references.csv                                # Metadata on fetched pseudo-references.
│
├── merge_mode/
│   ├── consensus/
│   │   ├── {sample}_r_{r}_s_{s}_con_{reference_name}.fas      # Individual consensus files for each sample/parameter combination. 
│   │   └── {run_name}_merge.fasta                             # Combined consensus sequences multi-FASTA for merge mode. sequences will have the fasta header '{sample}_r_{r}_s_{s}_con_{reference_name}_merge' after renaming.
│   ├── alignment/                                             # Alignment files containing trimmed reads aligned by MGE/exonerate (i.e reads that went into generating the consensus sequence).
│   ├── trimmed_data/                                          # Intermediate. Removed at end of run.
│   ├── fasta_cleaner/                                         # Consensus 'cleaning' intermediates and outputs. Each 01-05 directory contains metrics related to each cleaning step.
│   │   ├── 01_human_filtered/
│   │   ├── 02_at_filtered/
│   │   ├── 03_outlier_filtered/
│   │   ├── [optional] 04_reference_filtered/
│   │   ├── 05_cleaned_consensus/
│   │   ├── combined_statistics.csv                            # Combined statistics from all cleaning steps.
│   │   └── all_consensus_sequences.fasta                      # Cleaned consensus sequences multi-FASTA.
│   ├── logs/
│   │   ├── clean_headers.log                                  # Aggregated logs for fasta header cleaning for all samples.
│   │   ├── cleaning_complete.txt
│   │   ├── file_cleanup_complete.txt
│   │   ├── rename_fasta.log
│   │   ├── fasta_cleaner/                                     # Contains logs for each consensus sequence cleaning step.
│   │   ├── fastp/                                             # Contains fastp trimming logs for each sample.
│   │   ├── mge/
│   │   │   └── {sample}_r_{r}_s_{s}_{reference_name}/         # MGE vulgar files for each sample
│   │   │   └── alignment_files.log
│   │   │   └── mge_stats.log
│   ├── out/                                                   # MGE output files for each sample.
│   ├── err/                                                   # MGE error/processing logs for each sample.
│   ├── {run_name}_merge-stats.csv                             # Summary statistics file for merge mode. 
│   └── cleanup_complete.txt
│
├── concat_mode/
│   ├── consensus/
│   │   ├── {sample}_r_{r}_s_{s}_con_{reference_name}.fas       # Individual consensus files for each sample/parameter combination.
│   │   └── {run_name}_concat.fasta                            # Combined consensus sequences multi-FASTA for merge mode. sequences will have the fasta header '{sample}_r_{r}_s_{s}_con_{reference_name}' after renaming.
│   ├── alignment/                                             # Alignment files containing trimmed reads aligned by MGE/exonerate (i.e reads that went into generating the consensus sequence).
│   ├── trimmed_data/
│   │   ├── {sample}/                                          # Trimmed fwd and rev reads (fq.gz), fastp JSON and HTML reports, and Trim Galore QC report.
│   ├── fasta_cleaner/                                         # Consensus 'cleaning' intermediates and outputs. Each 01-05 directory contains metrics related to each cleaning step.
│   │   ├── 01_human_filtered/
│   │   ├── 02_at_filtered/
│   │   ├── 03_outlier_filtered/
│   │   ├── [optional] 04_reference_filtered/
│   │   ├── 05_cleaned_consensus/
│   │   ├── combined_statistics.csv                            # Combined statistics from all cleaning steps.
│   │   └── all_consensus_sequences.fasta                      # Cleaned consensus sequences multi-FASTA.
│   ├── logs/
│   │   ├── alignment_files.log
│   │   ├── concat/
│   │   │   └── {sample}.log
│   │   ├── concat_reads.log
│   │   ├── cleaning_complete.txt
│   │   ├── fasta_cleaner.log
│   │   ├── mge/
│   │   │   └── {sample}_vulgar_r_{r}_s_{s}_{reference_name}.txt
│   │   ├── rename_complete.txt
│   │   ├── trim_galore/
│   │   │   └── {sample}.log
│   ├── logs/
│   │   ├── concat_reads.log                                   # Aggregated logs for trimmed PE read concatentation.
│   │   ├── cleaning_complete.txt
│   │   ├── file_cleanup_complete.txt
│   │   ├── rename_fasta.log
│   │   ├── fasta_cleaner/                                     # Contains logs for each consensus sequence cleaning step.
│   │   ├── fastp/                                             # Contains fastp trimming logs for each sample.
│   │   ├── gzip/                                              # Logs for each trimmed read compression.
│   │   ├── mge/
│   │   │   └── {sample}_r_{r}_s_{s}_{reference_name}/         # MGE vulgar files for each sample
│   │   │   └── alignment_files.log
│   │   │   └── mge_stats.log
│   │   └── trim_galore.log
│   ├── out/                                                   # MGE output files for each sample.
│   ├── err/                                                   # MGE error/processing logs for each sample.
│   ├── {run_name}_concat-stats.csv                             # Summary statistics file for concat mode. 
│   └── cleanup_complete.txt
│
├── fasta_compare/
│   ├── {run_name}_fasta_compare.csv                          # Comparison results between modes
│   ├── {run_name}_full_sequences.fasta                       # Full consensus sequences from both modes
│   ├── {run_name}_barcode_sequences.fasta                    # Barcode sequences from both modes
│   └── fasta_compare_{run_name}.log
└── {run_name}_combined_stats.csv                             # Combined statistics from both modes
```


# Contributing #
- Please feel free to submit issues, fork the repository, and create pull requests for any improvements.
- This snakemake pipeline was produced by Dan Parsons @ NHMUK for BioDiversity Genomics Europe.
- Since this snakemake pipeline uses [MitogeneExtractor](https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.14075) at its core, please cite:
  Brasseur, M.V., Astrin, J.J., Geiger, M.F., Mayer, C., 2023. MitoGeneExtractor: Efficient extraction of mitochondrial genes from next-generation sequencing libraries. Methods in Ecology and Evolution.


  ## To do ##
- Split Snakefile into .smk files
- Add pre-MGE rule to subset/downsample very large input files (Based on file size or sequence number?)
  
