# Barcode Gene Extractor & Evaluator (BGEE) Snakemake workflow #
Snakemake workflow for recovering high-quality barcode sequences from genome skim data, built around MitoGeneExtractor and adapted for genome skims of museum specimens.

# Contents # 
 - [Requirements](#Requirements)
 - [Workflow](#Workflow)
 - [Installation and set up](#Installation-and-set-up)
 - [Cluster configuration](#Cluster-configuration-using-Snakemake-profiles)
 - [Results structure](#Results-structure)
 - [Contributing](#Contributing)

# Requirements #
- [MitoGeneExtractor](https://github.com/cmayer/MitoGeneExtractor) version 1.9.6 installed. Clone repository and follow installation instructons. 
- Paired-end reads in .fastq.gz or .fastq format.
- samples_file.csv (generated manually, or as outlined below if working from BOLD sample metadata).
- sequence_references_file.csv (generated manually, or using [Gene Fetch](https://github.com/bge-barcoding/gene_fetch?tab=readme-ov-file) within the workflow).
- Activated conda env (see bgee_env.yaml).

# Workflow #
1. **Preprocessing mode** (both pre-processing modes are run in parallel):
   - 'concat':
     - Initial raw read quality control - Adapter, trimming, quality trimming, poly-g trimming, and deduplication of paired-end reads using [fastp](https://github.com/OpenGene/fastp) (fastp_pe_concat).
     - Concantenation of trimmed PE reads (fastq_concat) and associated log files (aggregate_concat_logs).
     - Compress trimmed reads (gziped_trimmed).
     - Secondary read quality control - Additional quality trimming of concatenated reads using [trim galore](https://github.com/FelixKrueger/TrimGalore) (quality_trim) and concatenation of associated log files (Aggregate_trim_galore_logs).
   - 'merge:
     - Raw read quality control and merging - Adapter, trimming, quality trimming, poly-g trimming, deduplication, and merging of paired-end reads using [fastp](https://github.com/OpenGene/fastp) (fastp_pe_merge).
     - 'Clean' headers of input files, as required by MitoGeneExtractor (clean_headers_merge), and concentation of associated log files (Aggregate_clean_headers_logs).

<p align="center">
  <img src="https://github.com/user-attachments/assets/139b8c7c-b0dc-465c-8c95-e3a58ea1ab96" width="650"/>
</p>

2. **Sample-specific pseudo-reference retrieval** from GenBank using [Gene-Fetch](https://github.com/bge-barcoding/gene_fetch). (gene_fetch).
3. **Protein reference-guided barcode recovery** using [MitoGeneExtractor](https://github.com/cmayer/MitoGeneExtractor) (MitoGeneExtractor_concat & MitoGeneExtractor_merge).
4. **'Raw' consensus sequence header manipulation and concatenation** into multi-FASTA (rename_and_combine_cons_concat & rename_and_combine_cons_merge) (uses supplementary [rename_headers.py](https://github.com/bge-barcoding/MitoGeneExtractor-BGE/blob/main/workflow/scripts/rename_headers.py).
5. **Remove exonerate intermediates** if remaining after MGE (remove_exonerate_intermediates).
6. **Generate list of MGE alignment files** for downstream processing (create_alignment_log).
7. **Filter MGE alignment files to remove low-quality, contaminant, or outlier sequences before repeat consensus sequence generation.**
   - Remove aligned reads with similarity to human COI (a common contaminant of museum specimens) (uses supplementary [01_human_cox1_filter.py](https://github.com/bge-barcoding/MitoGeneExtractor-BGE/blob/main/workflow/scripts/01_human_cox1_filter.py)).
   - Remove aligned reads with AT content over and/or below a specified threshold (high AT content can be indicative of contamination (e.g. fungal)) (uses supplementary [02_at_content_filter.py](https://github.com/bge-barcoding/MitoGeneExtractor-BGE/blob/main/workflow/scripts/02_at_content_filter.py)).
   - Remove reads that are statistical outliers compared to the original consensus (uses supplementary [03_statistical_outliers.py](https://github.com/bge-barcoding/MitoGeneExtractor-BGE/blob/main/workflow/scripts/03_statistical_outlier_filter.py)).
   - [optional] Remove reads with similarity to supplied reference sequence(s) (uses supplementary [04_reference_filter.py](https://github.com/bge-barcoding/MitoGeneExtractor-BGE/blob/main/workflow/scripts/04_reference_filter.py)).
   - Generation of 'cleaned' consensus sequence (uses supplementary [05_consensus_generator.py](https://github.com/bge-barcoding/MitoGeneExtractor-BGE/blob/main/workflow/scripts/05_consensus_generator.py)).
   - Aggregate metrics from each stage of filtering (uses supplementary [06_aggregate_filter_metrics.py](https://github.com/bge-barcoding/MitoGeneExtractor-BGE/blob/main/workflow/scripts/06_aggregate_filter_metrics.py)).
   - Remove intermediate files and unecessary logs generated during the consensus cleaning process (remove_fasta_cleaner_files).
  
<p align="center">
  <img width="285" height="443" alt="image" src="https://github.com/user-attachments/assets/957d43a7-0c00-40ce-bab8-1827d0e37e1b" />
</p>

8. **Validate and triage barcodes** 
   - Structural validation - Processes all barcode consensus sequences produced by MGE (and fasta_cleaner), combining structural analysis (barcode region HMM-alignment, gaps, ambiguous bases, sequence length) with functional analysis (reading frame determination, stop codon counting, protein translation) to ensure the 'best' consensus sequence is selected for each barcode (uses [supplementary structural_validation.py](https://github.com/bge-barcoding/BGEE/blob/main/workflow/scripts/structural_validation.py))
   - Local BLASTn search - of structurally validated barcode sequences (currently uses BOLDistilled for rapid BLASTn searches, and thus is only suitable for COI) (uses supplementary [tv_local_blast.py](https://github.com/bge-barcoding/BGEE/blob/main/workflow/scripts/tv_local_blast.py)).
   - Taxonomic validation of BLAST results - through parsing output local BLASTn results, and comparison of top BLAST hit taxonomy (observed taxonomy) with input (expected) taxonomy (uses supplementary [tv_blast2taxonomy.py](https://github.com/bge-barcoding/BGEE/blob/main/workflow/scripts/tv_blast2taxonomy.py)).
9. **Compile statistics** from read QC, MGE, and consensus cleaning metrics into a CSV report for both 'concat' and 'merge' modes (uses supplementary [mge_stats.py](https://github.com/bge-barcoding/MitoGeneExtractor-BGE/blob/main/workflow/scripts/mge_stats.py) and combine_stats_files).
10. **Final mergin of metrics** generated during Fastp, TrimGalore, Gene Fetch, MGE, Fasta_cleaner, structural validation and taxononmic validation steps (uses supplementary [val_csv_merger.py](https://github.com/bge-barcoding/BGEE/blob/main/workflow/scripts/val_csv_merger.py)).
11. **Clean up** temporary files, sample-specific logs once aggregated, etc. (cleanup_files).



# Installation and set up: #
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
- **Must contain `ID`, `forward` (read paths), `reverse` (read paths), and `taxid` _OR_ `hierarchical taxonomy` columns (see below for examples).**
- Due to regex matching and statistics aggregation, the sample ID will be considered as the string before the first underscore. **It is therefore recommended that sample names do not use '_' characters.** E.g. BSNHM002-24 instead of BSNHM002_24, or P3-1-A10-2-G1 instead of P3_1_A10_2_G1.
- Taxid's can be found manually by searching the expected species/genus/family of each sample in the [NCBI taxonomy database](https://www.ncbi.nlm.nih.gov/taxonomy).
  
**samples.csv example (taxid)**
| ID | forward | reverse | taxid |
| --- | --- | --- | --- |
| BSNHM002-24  | abs/path/to/R1.fq.gz | abs/path/to/R2.fq.gz | 177658 |
| BSNHM038-24 | abs/path/to/R1.fq.gz | abs/path/to/R2.fq.gz | 177627 |
| BSNHM046-24 | abs/path/to/R1.fq.gz | abs/path/to/R2.fq.gz | 3084599 |

**samples.csv example (hierarchical taxonomy)**
| ID | forward | reverse | phylum | class | order | family | genus | species |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | 
| BSNHM002-24  | abs/path/to/R1.fq.gz | abs/path/to/R2.fq.gz | Arthropoda | Insecta | Hemiptera | Cicadidae | Tibicina | Tibicina tomentosa |
| BSNHM038-24 | abs/path/to/R1.fq.gz | abs/path/to/R2.fq.gz | Tracheophyta | Pinopsida | Pinales | Pinaceae | Abies |  |
| BSNHM046-24 | abs/path/to/R1.fq.gz | abs/path/to/R2.fq.gz | Annelida | Polychaeta | Terebellida | Ampharetidae | Samytha | Samytha sexcirrata |

## Gathering sample-specific pseudo-references ##
- This can be created manually, or using [Gene-fetch](https://github.com/bge-barcoding/gene_fetch) integrated into the workflow (highly recommended). If enabled (in the config.yaml by setting `run_gene_fetch` to 'true'), gene-fetch will retrieve the necessary protein pseudo-references for each sample from NCBI GenBank using the sample's taxonomic identifier (taxid)/taxonomic hierarchy for each sample, a sequence target (e.g. COI or rbcL), and NCBI API credentials (email address & API key - see [guidance](https://support.nlm.nih.gov/kbArticle/?pn=KA-05317) on getting a key).
- If hierarchical taxonomic information is provided (see `samples.csv example (hierarchical taxonomy)` above) for each sample instead of a taxid (see `samples.csv example (taxid)` above), `Gene Fetch` will, starting from the lowest given rank (e.g. species), find the closest valid taxid on NCBI for the provided taxonomy before proceeding with pseudo-reference fetching.
- **Must contain 'process_id', 'reference'name' and 'protein_reference_path' at a minimum.**

**sample_references.csv example**
| process_id | reference_name | protein_reference_path | 
| --- | --- | --- |
| BSNHM002-24  | BSNHM002-24 | path/to/BSNHM002-24.fasta |
| BSNHM038-24 | BSNHM038-24 | path/to/BSNHM038-24.fasta |
| BSNHM046-24 | BSNHM046-24 | path/toBSNHM046-24.fasta |
* **Currently, it is crucial that the sample name (process_id), reference sequence FASTA file, and corresponding reference sequence FASTA header are all identical for correct sample-reference file mapping.**

## Customising snakemake configuration file ##
- Update [config/config.yaml](https://github.com/bge-barcoding/BGEE/blob/main/config/config.yaml) with the neccessary paths and variables.
- **Currently, BGEE's structural_validation only works for COI-5P and rbcL barcodes due to HMM availability (to be updated). BGEE's taxonomic_validation currently only works for COI-5P due to sole use of the BOLDistilled BLAST DB (to be updated).**
- Each rule in the config.yaml can specify the number of requested threads and memory resources (in Mb) for every job (e.g. specifying 4 threads and 4G memory for fastp_pe_merge would allocate those resources for every 'fastp_pe_merge' job).

## Cluster configuration using Snakemake profiles ##
- See `profiles/` directory for 'slurm' and 'local' (i.e. non-SLURM) cluster submission parameters. The `jobs` parameter is likely the most important as it dictates how many workflow jobs can be run concurrently.
- The profile (`profiles/local` or `profiles/slurm`) will need to be changed depending on your system (see `$PROFILE` variable in `snakemake_run.sh`).

### Cluster submission ###
- [snakemake_run.sh](https://github.com/bge-barcoding/MitoGeneExtractor-BGE/blob/main/snakemake_run.sh) handles submission of the snakemake workflow to the HPC cluster. The working directory will initially be unlocked (using `--unlock`) and then the snakemake workflow will be run. 


#### If using `profiles/slurm`, SLURM will orchestrate submission of each step in the workflow as a separate job: ####
- The safer way is to launch Snakemake inside a persistent screen session. This ensures the workflow keeps running even if you disconnect.
```
# Start a persistent screen session
screen -S [SESSION_NAME]

# Allocate resources for Snakemake to use
salloc --job-name=[SESSION_NAME] \
       --partition=[YOUR_PARTITION] \
       --cpus-per-task=16 \
       --mem=8G
```
- Inside this allocation, you can launch Snakemake within an interactive session with srun:
```
srun ./snakemake_run.sh
```
- To detach from the screen session (disconnect but keep it running): `Ctrl + A + D`
- To reconnect again: `screen -r [SESSION_NAME]`
- You may see a warning such as "You are running snakemake in a SLURM job context. This is not recommended..." - This can generally be ignored because the salloc session is only acting as a submission manager. If you do encounter problems, try running `./snakemake_run.sh` within the screen session without running `salloc`.


#### If using `profiles/local`, all workflow steps will be run as a single job: ####
-  Simply run `./snakemake_run.sh` on your desired cluster compute node. This node will handle all job scheduling and job computation.



# Results structure #
```
output_dir/
├── 01_preprocessing/
│   ├── merge_mode/
│   │   ├── trimmed_data/
│   │   │   ├── {sample}_merged.fq                             # Merged paired-end reads
│   │   │   ├── {sample}_merged_clean.fq                       # Header-cleaned merged reads
│   │   │   ├── {sample}_fastp_report.html                     # FastP HTML report
│   │   │   ├── {sample}_fastp_report.json                     # FastP JSON report
│   │   │   └── unpaired/                                      # Unpaired reads from merging
│   │   └── logs/
│   │       ├── clean_headers/
│   │       │   └── clean_headers.log                          # Aggregated header cleaning logs
│   │       ├── fastp/                                         # Individual FastP logs per sample
│   │       └── final_cleanup_complete.txt
│   └── concat_mode/
│       ├── trimmed_data/
│       │   └── {sample}/
│       │       ├── {sample}_R1_trimmed.fastq.gz               # Trimmed forward reads
│       │       ├── {sample}_R2_trimmed.fastq.gz               # Trimmed reverse reads
│       │       ├── {sample}_concat_trimmed.fq                 # Quality-trimmed concatenated reads
│       │       ├── {sample}_fastp_report.html                 # FastP HTML report
│       │       ├── {sample}_fastp_report.json                 # FastP JSON report
│       │       └── {sample}_concat.fastq_trimming_report.txt  # Trim Galore report
│       └── logs/
│           ├── concat/
│           │   └── concat_reads.log                           # Aggregated concatenation logs
│           ├── trim_galore/
│           │   └── trim_galore.log                            # Aggregated Trim Galore logs
│           ├── fastp/                                         # Individual FastP logs per sample
│           ├── gzip/                                          # Compression logs per sample
│           └── final_cleanup_complete.txt
│
├── 02_references/                                             # Only if run_gene_fetch = true
│   ├── protein/
│   │   └── {sample}.fasta                                     # Protein references for each sample
│   ├── genbank/                                               # GenBank records (if genbank: true)
│   └── sequence_references.csv                                # Reference metadata
│
├── 03_barcode_recovery/
│   ├── merge_mode/
│   │   ├── alignment/
│   │   │   └── {sample}_r_{r}_s_{s}_align_{sample}.fas        # MGE alignment files
│   │   ├── consensus/
│   │   │   ├── {sample}_r_{r}_s_{s}_con_{sample}.fas          # Individual consensus files
│   │   │   └── {run_name}_cons_combined-merge.fasta           # Combined consensus sequences
│   │   ├── fasta_cleaner/
│   │   │   ├── 01_human_filtered/
│   │   │   │   ├── human_filtered.txt                         # List of filtered files
│   │   │   │   └── human_filter_metrics.csv                   # Human filtering metrics
│   │   │   ├── 02_at_filtered/
│   │   │   │   ├── at_filtered_sequences/                     # Individual filtered files
│   │   │   │   ├── at_filtered.txt                            # List of filtered files
│   │   │   │   └── at_filter_summary.csv                      # AT filtering summary
│   │   │   ├── 03_outlier_filtered/
│   │   │   │   ├── outlier_filtered.txt                       # List of filtered files
│   │   │   │   ├── outlier_filter_summary_metrics.csv         # Summary metrics
│   │   │   │   └── outlier_filter_individual_metrics.csv      # Individual metrics
│   │   │   ├── 04_reference_filtered/                         # Optional - if reference filtering enabled
│   │   │   │   ├── reference_filtered.txt                     # List of filtered files
│   │   │   │   └── reference_filter_metrics.csv               # Reference filtering metrics
│   │   │   ├── 05_cleaned_consensus/
│   │   │   │   └── cleaned_cons_metrics-merge.csv             # Consensus generation metrics
│   │   │   ├── combined_statistics.csv                        # Aggregated cleaning statistics
│   │   │   └── cleaned_cons_combined.fasta                    # Final cleaned consensus sequences
│   │   ├── logs/
│   │   │   ├── mge/
│   │   │   │   ├── alignment_files.log                        # List of alignment files
│   │   │   │   ├── mge_stats.log                              # MGE statistics log
│   │   │   │   └── {sample}_r_{r}_s_{s}/                      # MGE vulgar files per sample
│   │   │   ├── fasta_cleaner/
│   │   │   │   └── fasta_cleaner_complete.txt                 # Cleaner completion flag
│   │   │   ├── rename_consensus/
│   │   │   │   └── rename_fasta.log                           # Header renaming logs
│   │   │   ├── fasta_cleaner_complete.txt                     # Main cleaner completion
│   │   │   └── exonerate_int_cleanup_complete.txt             # Intermediate cleanup completion
│   │   ├── out/                                               # MGE output files per sample/parameter
│   │   ├── err/                                               # MGE error logs per sample/parameter
│   │   └── {run_name}_merge-stats.csv                         # Mode-specific statistics
│   └── concat_mode/
│       ├── alignment/
│       │   └── {sample}_r_{r}_s_{s}_align_{sample}.fas        # MGE alignment files
│       ├── consensus/
│       │   ├── {sample}_r_{r}_s_{s}_con_{sample}.fas          # Individual consensus files
│       │   └── {run_name}_cons_combined-concat.fasta          # Combined consensus sequences
│       ├── fasta_cleaner/
│       │   ├── 01_human_filtered/
│       │   │   ├── human_filtered.txt                         # List of filtered files
│       │   │   └── human_filter_metrics.csv                   # Human filtering metrics
│       │   ├── 02_at_filtered/
│       │   │   ├── at_filtered_sequences/                     # Individual filtered files
│       │   │   ├── at_filtered.txt                            # List of filtered files
│       │   │   └── at_filter_summary.csv                      # AT filtering summary
│       │   ├── 03_outlier_filtered/
│       │   │   ├── outlier_filtered.txt                       # List of filtered files
│       │   │   ├── outlier_filter_summary_metrics.csv         # Summary metrics
│       │   │   └── outlier_filter_individual_metrics.csv      # Individual metrics
│       │   ├── 04_reference_filtered/                         # Optional - if reference filtering enabled
│       │   │   ├── reference_filtered.txt                     # List of filtered files
│       │   │   └── reference_filter_metrics.csv               # Reference filtering metrics
│       │   ├── 05_cleaned_consensus/
│       │   │   └── cleaned_cons_metrics-concat.csv            # Consensus generation metrics
│       │   ├── combined_statistics.csv                        # Aggregated cleaning statistics
│       │   └── cleaned_cons_combined.fasta                    # Final cleaned consensus sequences
│       ├── logs/
│       │   ├── mge/
│       │   │   ├── alignment_files.log                        # List of alignment files
│       │   │   ├── mge_stats.log                              # MGE statistics log
│       │   │   └── {sample}_r_{r}_s_{s}/                      # MGE vulgar files per sample
│       │   ├── fasta_cleaner/
│       │   │   └── fasta_cleaner_complete.txt                 # Cleaner completion flag
│       │   ├── rename_consensus/
│       │   │   └── rename_fasta.log                           # Header renaming logs
│       │   ├── fasta_cleaner_complete.txt                     # Main cleaner completion
│       │   └── exonerate_int_cleanup_complete.txt             # Intermediate cleanup completion
│       ├── out/                                               # MGE output files per sample/parameter
│       ├── err/                                               # MGE error logs per sample/parameter
│       └── {run_name}_concat-stats.csv                        # Mode-specific statistics
│   ├── {run_name}_bgee_stats.csv                              # Combined statistics from both modes
│   └── {run_name}_all_cons_combined.fasta                     # All consensus sequences from both modes
│
├── 04_barcode_validation/
│   ├── structural/                                            # Only if run_structural_validation = true
│   │   ├── structural_validation.csv                          # Structural validation results
│   │   ├── {run_name}_full_sequences.fasta                    # Full sequences passing validation
│   │   └── {run_name}_barcode_sequences.fasta                 # Barcode sequences passing validation
│   └── taxonomic/                                             # Only if run_taxonomic_validation = true (& run_structural_validation = true)
│       ├── 01_local_blast_output.csv                          # BLAST results
│       ├── 02_taxonomic_validation.csv                        # Taxonomic validation results
│       └── {run_name}_barcode_sequences.fasta                 # Final validated barcode sequences
│
├── {run_name}_final_validated_barcodes.fasta                  # Only if both validations run
├── {run_name}_final_stats.csv                                 # Only if both validations run
└── logs/                                                      # Top-level logs directory
```


# Contributing #
- Please feel free to submit issues, fork the repository, and create pull requests for any improvements.
- This snakemake pipeline was produced by Dan Parsons @ NHMUK for the BioDiversity Genomics Europe (BGE) consortium.
- Since BGEE uses [MitogeneExtractor](https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.14075) at its core, please cite:
  Brasseur, M.V., Astrin, J.J., Geiger, M.F., Mayer, C., 2023. MitoGeneExtractor: Efficient extraction of mitochondrial genes from next-generation sequencing libraries. Methods in Ecology and Evolution.


# To do #
- Split Snakefile into modular .smk files.
- Increase flexibility of input CSV headers (e.g. ID column in sample.csv and process_id column in sequence_references.csv could be ID/id/Process ID/PROCESS ID/process_id/sample/sample_id/SAMPLE ID/etc).
- Update 01_human_cox1_filter.py to not just filter against human coi, but the whole human genome/mitogenome.
  
