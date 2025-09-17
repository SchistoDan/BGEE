#!/usr/bin/env python3
"""
val_csv_merger.py

This script merges structural validation (structval) and taxonomic validation (taxval) 
CSV data files into a BGEE metrics CSV file. It performs left joins based on 
sequence ID matching to create a comprehensive dataset combining gene expression data 
with structural and taxonomic validation metrics.

EXPECTED DATA STRUCTURE
=======================
Structval CSV columns:
- seq_id: Sequence identifier (matches cleaned FASTA headers)
- barcode_length, barcode_ambiguous_bases_original, barcode_ambiguous_bases
- barcode_longest_stretch, reading_frame, stop_codons, barcode_rank
- best_sequence, selected_barcode_fasta

Taxval CSV columns:
- seq_id: Sequence identifier (matches cleaned FASTA headers)  
- taxval_rank, expected_taxonomy, expected_taxonomy_rank, obs_taxonomy
- match_taxonomy, matched_rank, top_matching_hit, pident
- length (renamed to aligned_length), mismatch, evalue

BGEE CSV:
- fasta_header: FASTA sequence headers (may include '>' prefix)
- Additional gene expression and metadata columns

USAGE
=====
Command line interface:

python val_csv_merger.py --input-structval structval_data.csv \\
                     --input-taxval taxval_data.csv \\
                     --input-bgee-csv bgee_expression.csv \\
                     --output merged_results.csv


OUTPUT
======
Creates a merged CSV file containing:
- All original columns from the Bgee dataset
- Structural validation metrics (where sequence matches exist)
- Taxonomic validation results (where sequence matches exist)
- Merge statistics printed to console

DEPENDENCIES
============
- pandas
- Python 3.6+

AUTHORS
-------
Dan Parsons @NHMUK
License: MIT
Version: 1.0

"""


import pandas as pd
import argparse
import sys
from pathlib import Path


def clean_fasta_header(header):
    """Remove '>' prefix from fasta header if present."""
    if pd.isna(header):
        return header  # Return NaN as-is
    return str(header).lstrip('>')


def merge_csv_files(input_structval, input_taxval, input_bgee_csv, output_file):
    """
    Merge structval and taxval CSV files into bgee CSV file.
    
    Args:
        input_structval (str): Path to structval CSV file
        input_taxval (str): Path to taxval CSV file  
        input_bgee_csv (str): Path to bgee CSV file
        output_file (str): Path for output merged CSV file
    """
    
    print("Loading CSV files...")
    
    # Load CSV files
    try:
        structval_df = pd.read_csv(input_structval)
        taxval_df = pd.read_csv(input_taxval)
        bgee_df = pd.read_csv(input_bgee_csv)
    except FileNotFoundError as e:
        print(f"Error: Could not find file {e.filename}")
        sys.exit(1)
    except pd.errors.EmptyDataError as e:
        print(f"Error: Empty CSV file - {e}")
        sys.exit(1)
    
    print(f"Loaded {len(structval_df)} rows from structval")
    print(f"Loaded {len(taxval_df)} rows from taxval") 
    print(f"Loaded {len(bgee_df)} rows from bgee")
    
    # Clean fasta headers in bgee_df for matching
    bgee_df['clean_fasta_header'] = bgee_df['fasta_header'].apply(clean_fasta_header)
    
    # Define columns to merge from each file
    structval_columns = [
        'seq_id', 'barcode_length', 
        'barcode_ambiguous_bases_original', 'barcode_ambiguous_bases', 
        'barcode_base_count', 'reading_frame', 'stop_codons', 
        'barcode_rank', 'best_sequence', 'selected_barcode_fasta'
    ]
    
    taxval_columns = [
        'seq_id', 'taxval_rank', 'expected_taxonomy', 'expected_taxonomy_rank', 
        'obs_taxonomy', 'match_taxonomy', 'matched_rank', 'top_matching_hit', 
        'pident', 'length', 'mismatch', 'evalue'
    ]
    
    # Check if required columns exist
    missing_structval = [col for col in structval_columns if col not in structval_df.columns]
    missing_taxval = [col for col in taxval_columns if col not in taxval_df.columns]
    
    if missing_structval:
        print(f"Warning: Missing columns in structval file: {missing_structval}")
    if missing_taxval:
        print(f"Warning: Missing columns in taxval file: {missing_taxval}")
    
    # No renaming needed for structval (no length column), rename taxval length to aligned_length
    structval_df_subset = structval_df[structval_columns].copy()
    
    taxval_df_subset = taxval_df[taxval_columns].copy()  
    taxval_df_subset = taxval_df_subset.rename(columns={'length': 'aligned_length'})
    
    # Merge bgee with structval
    print("Merging with structval data...")
    merged_df = pd.merge(
        bgee_df, 
        structval_df_subset, 
        left_on='clean_fasta_header', 
        right_on='seq_id',
        how='left',
        suffixes=('', '_structval')
    )
    
    # Merge with taxval
    print("Merging with taxval data...")
    final_df = pd.merge(
        merged_df,
        taxval_df_subset,
        left_on='clean_fasta_header',
        right_on='seq_id',
        how='left',
        suffixes=('', '_taxval')
    )
    
    # Clean up temporary columns
    final_df = final_df.drop(columns=['clean_fasta_header'], errors='ignore')
    final_df = final_df.drop(columns=['seq_id_structval', 'seq_id_taxval'], errors='ignore')
    
    # Report merge statistics
    structval_matches = final_df['barcode_length'].notna().sum()
    taxval_matches = final_df['taxval_rank'].notna().sum()
    
    print(f"\nMerge Statistics:")
    print(f"Total rows in final output: {len(final_df)}")
    print(f"Rows with structval matches: {structval_matches}")
    print(f"Rows with taxval matches: {taxval_matches}")
    print(f"Rows with both matches: {final_df[final_df['barcode_length'].notna() & final_df['taxval_rank'].notna()].shape[0]}")
    
    # Save merged file
    final_df.to_csv(output_file, index=False)
    print(f"\nMerged data saved to: {output_file}")
    
    return final_df


def main():
    parser = argparse.ArgumentParser(
        description="Merge structval and taxval CSV files into BGEE CSV file"
    )
    parser.add_argument(
        '--input-structval', 
        required=True,
        help='Path to structval CSV file'
    )
    parser.add_argument(
        '--input-taxval', 
        required=True,
        help='Path to taxval CSV file'
    )
    parser.add_argument(
        '--input-bgee-csv', 
        required=True,
        help='Path to bgee CSV file'
    )
    parser.add_argument(
        '--output', 
        required=True,
        help='Path for output merged CSV file'
    )
    
    args = parser.parse_args()
    
    # Verify input files exist
    for file_path in [args.input_structval, args.input_taxval, args.input_bgee_csv]:
        if not Path(file_path).exists():
            print(f"Error: File not found: {file_path}")
            sys.exit(1)
    
    # Run the merge
    merge_csv_files(
        args.input_structval,
        args.input_taxval, 
        args.input_bgee_csv,
        args.output
    )


if __name__ == "__main__":
    main()
