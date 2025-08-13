"""
FASTA Sequence Analysis and Summary Tool
--------------------------------------

This script analyses a log file containing a list of alignment FASTA files, .out files
containing 'raw' summary stats from MGE, and cleaning statistics from a CSV file.
It generates comprehensive summary statistics in CSV format with separate rows for 
alignment data (from FASTA files) and cleaning data (from cleaning CSV).

Usage:
    python mge_stats.py -a/--alignment_log <log_file> -o/--output <output_csv> -od/--out_file_dir <out_file_dir> -c/--cleaning_csv <cleaning_csv>

Arguments:
    -a, --alignment_log : str
        Path to a text file containing a list of FASTA file paths (one per line)
    -o, --output : str
        Name of the output CSV file
    -od, --out_file_dir : str
        Directory containing .out files with additional sequence statistics
    -c, --cleaning_csv : str (optional)
        Path to CSV file containing cleaning statistics
    --ref_seqs : str (optional)
        Path to CSV file containing reference sequence information (taxid, protein_accession, matched_rank)
    --pre-fasta : str (optional)
        Path to multi-FASTA file containing pre-cleaning sequences for header extraction
    --post-fasta : str (optional)
        Path to multi-FASTA file containing post-cleaning sequences for header extraction

Auto-detection:
    The script automatically detects "merge_mode" in file paths and appends "_merge" to 
    the mge_params column accordingly.

Outputs:
    - <output_file>.csv: Main summary file containing all statistics
    - mge_stats.log: Log file with processing information

The script generates the following metrics for each sample:
    - Filename: Identifier extracted from the filename
    - ID: Identifier extracted from the filename
    - mge_params: Parameters used for MGE (e.g., r_1.3_s_50)
    - fasta_header: FASTA header from pre/post-fasta files (if provided)
    - sample_taxid: Taxonomic ID of the sample from reference sequences
    - ref_accession: Protein accession number from reference sequences  
    - ref_rank: Matched taxonomic rank from reference sequences
    - n_reads_in: Number of input sequences (from .out file)
    - n_reads_aligned: Number of aligned sequences (merged: FASTA alignment count or cleaning kept_reads)
    - n_reads_skipped: Number of sequences that were successfully aligned but not in the FASTA file
    - ref_length: Length of alignment (from .out file)
    - cov_min: Minimum coverage depth (merged: alignment or cleaning data)
    - cov_max: Maximum coverage depth (merged: alignment or cleaning data)
    - cov_avg: Average coverage depth (merged: alignment or cleaning data)
    - cov_med: Median coverage depth (merged: alignment or cleaning data)
    - cleaning_removed_human: Number of sequences removed due to human similarity
    - cleaning_removed_at: Number of sequences removed due to AT content
    - cleaning_removed_outlier: Number of sequences removed as statistical outliers
    - cleaning_ambig_bases: Number of ambiguous bases after cleaning
    - cleaning_cov_percent: Coverage percentage after cleaning
"""



import os
import csv
import argparse
import re
import logging
from Bio import SeqIO
import numpy as np
import glob


# Set up logging
logger = logging.getLogger('mge_stats')
logger.setLevel(logging.INFO)
# Use mode='w' to overwrite the log file each time
file_handler = logging.FileHandler('mge_stats.log', mode='w')
file_handler.setLevel(logging.INFO)
console_handler = logging.StreamHandler()
console_handler.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
file_handler.setFormatter(formatter)
console_handler.setFormatter(formatter)
logger.addHandler(file_handler)
logger.addHandler(console_handler)

def extract_sample_info(filename):
    """
    Extract the process ID and parameters from FASTA filenames.
    E.g., from "BGSNL096-23_r_1.3_s_50_align_BGSNL096-23.fas"
    Returns:
        - base_id: "BGSNL096-23"
        - full_id: "BGSNL096-23_r_1.3_s_50"
        - params: "r_1.3_s_50"
    """
    filename_no_ext = os.path.splitext(filename)[0]
    
    # Extract base ID (everything before first underscore)
    base_id_match = re.match(r'^([^_]+)', filename_no_ext)
    base_id = base_id_match.group(1) if base_id_match else None
    
    # Extract parameters (r_X_s_Y pattern)
    params_match = re.search(r'(r_[0-9.]+_s_[0-9.]+)', filename_no_ext)
    params = params_match.group(1) if params_match else ""
    
    # Combine to form the full ID used for matching files
    full_id = f"{base_id}_{params}" if params else base_id
    
    return base_id, full_id, params

def extract_cleaning_info(sample_name):
    """
    Extract the process ID and parameters from cleaning CSV sample names.
    E.g., from "BSNHM002-24_r_1.3_s_50_BSNHM002-24_fcleaner_merge"
    Returns:
        - base_id: "BSNHM002-24"  
        - mge_params: "r_1.3_s_50_fcleaner_merge"
    """
    logger.debug(f"extract_cleaning_info called with: '{sample_name}'")
    
    # Extract base ID (everything before first underscore)
    base_id_match = re.match(r'^([^_]+)', sample_name)
    base_id = base_id_match.group(1) if base_id_match else None
    
    # Extract r_X_s_Y pattern
    params_match = re.search(r'(r_[0-9.]+_s_[0-9.]+)', sample_name)
    basic_params = params_match.group(1) if params_match else ""
    
    # Debug for decimal samples
    if 'r_' in sample_name:
        logger.info(f"Cleaned consensus mge_params extracting for: '{sample_name}'")
        logger.info(f"  base_id: '{base_id}'")
        logger.info(f"  basic_params: '{basic_params}'")
    
    if not basic_params:
        logger.warning(f"No r_X_s_Y pattern found in cleaning sample: {sample_name}")
        return base_id, ""
    
    # Find the suffix (everything after the second occurrence of base_id)
    base_id_pattern = f"_{base_id}_"
    second_occurrence = sample_name.find(base_id_pattern)
    
    if second_occurrence != -1:
        suffix = sample_name[second_occurrence + len(base_id_pattern):]
        mge_params = f"{basic_params}_{suffix}" if suffix else basic_params
        logger.debug(f"Found suffix: '{suffix}', final mge_params: '{mge_params}'")
    else:
        mge_params = basic_params
        logger.warning(f"Could not find base_id pattern '{base_id_pattern}' in {sample_name}")
        
    return base_id, mge_params

def extract_process_id(filename):
    """Extract the process ID from the filename (for backward compatibility)."""
    base_id, _, _ = extract_sample_info(filename)
    return base_id

def extract_id_and_params_from_header(header):
    """
    Extract ID and mge_params from a FASTA header.
    E.g., from ">BSNHM004-24_r_1.5_s_50_BSNHM004-24_merge"
    Returns: ('BSNHM004-24', 'r_1.5_s_50_merge')
    """
    # Remove the '>' if present
    clean_header = header.lstrip('>')
    
    # Extract base ID (everything before first underscore)
    base_id_match = re.match(r'^([^_]+)', clean_header)
    base_id = base_id_match.group(1) if base_id_match else None
    
    if not base_id:
        return None, None
    
    # Extract r_X_s_Y pattern and everything after it
    params_match = re.search(r'(r_[0-9.]+_s_[0-9.]+.*)$', clean_header)
    if not params_match:
        return base_id, None
    
    # The full params string includes everything after r_X_s_Y
    full_params = params_match.group(1)
    
    # Remove the duplicate base_id if it appears in the params
    # Handle both cases: with trailing underscore (merge mode) and without (concat mode)
    base_id_with_underscore = f"_{base_id}_"  # "_BSNHM004-24_"
    base_id_at_end = f"_{base_id}$"           # "_BSNHM004-24" at end of string
    
    if base_id_with_underscore in full_params:
        # Merge mode: r_X_s_Y_base_id_suffix -> r_X_s_Y_suffix
        parts = full_params.split(base_id_with_underscore)
        if len(parts) >= 2:
            r_params = parts[0]  # e.g., "r_1.5_s_50"
            suffix = parts[1]    # e.g., "merge" or "fcleaner_merge"
            mge_params = f"{r_params}_{suffix}" if suffix else r_params
        else:
            mge_params = full_params
    elif re.search(base_id_at_end, full_params):
        # Concat mode: r_X_s_Y_base_id -> r_X_s_Y
        mge_params = re.sub(base_id_at_end, '', full_params)
    else:
        mge_params = full_params
    
    return base_id, mge_params

def parse_fasta_headers(fasta_file_path):
    """
    Parse FASTA file and create a lookup dictionary based on (ID, mge_params) combination.
    Returns a dictionary mapping (ID, mge_params) -> full_header_with_>.
    """
    header_lookup = {}
    
    if not fasta_file_path or not os.path.exists(fasta_file_path):
        return header_lookup
    
    try:
        with open(fasta_file_path, 'r', encoding='utf-8', errors='replace') as handle:
            for record in SeqIO.parse(handle, 'fasta'):
                # Extract ID and mge_params from the header
                header_id, mge_params = extract_id_and_params_from_header(record.id)
                
                if header_id and mge_params:
                    # Create lookup key from (ID, mge_params)
                    lookup_key = (header_id, mge_params)
                    full_header = f">{record.id}"
                    header_lookup[lookup_key] = full_header
                    logger.debug(f"Parsed header: {record.id} -> ID: {header_id}, params: {mge_params}")
                else:
                    logger.warning(f"Could not parse header: {record.id}")
        
        logger.info(f"Parsed {len(header_lookup)} headers from {fasta_file_path}")
        
    except Exception as e:
        logger.error(f"Error parsing FASTA headers from {fasta_file_path}: {e}")
    
    return header_lookup

def parse_out_file(out_file_path):
    """Parse the .out file for additional statistics."""
    try:
        with open(out_file_path, 'r') as file:
            content = file.read()
    except Exception as e:
        logger.error(f"Error reading .out file {out_file_path}: {e}")
        return None

    process_id_data = {
        'n_reads_in': None,
        'ref_length': None,
        'successful_aligned': None
    }

    # Extract number of input sequences (updated pattern)
    reads_in_match = re.search(r'Number of input sequences:\s*(\d+)', content)
    if reads_in_match:
        process_id_data['n_reads_in'] = int(reads_in_match.group(1))
    
    # Simplified pattern - just target the line with the number
    aligned_match = re.search(r'to the amino acid sequence, see vulgar file:\s*(\d+)', content)
    if aligned_match:
        process_id_data['successful_aligned'] = int(aligned_match.group(1))
    
    # Reference length pattern
    ref_length_match = re.search(r'Length of alignment:\s*(\d+)', content)
    if ref_length_match:
        process_id_data['ref_length'] = int(ref_length_match.group(1))
 
    return process_id_data

def parse_cleaning_csv(file_path):
    """Parse the cleaning CSV file for filtering statistics."""
    cleaning_data = {}
    try:
        with open(file_path, 'r', newline='') as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                # Extract sample name which should match with full_id
                sample_name = row.get('sample_name', '')
                if sample_name:
                    # Map the CSV columns to our output structure
                    cleaning_stats = {
                        'removed_human': int(row.get('removed_human', 0)),
                        'removed_at': int(row.get('removed_at_distance', 0)),
                        'removed_outlier': int(row.get('removed_outliers', 0)),
                        'removed_reference': int(row.get('removed_reference', 0)),
                        'kept_reads': int(row.get('cleaned_reads', 0)),
                        'ambig_bases': int(row.get('final_ambig_bases', 0)),
                        'cov_percent': float(row.get('cov_percent', 0)),
                        'cov_avg': float(row.get('cov_avg', 0)),
                        'cov_max': float(row.get('cov_max', 0)),
                        'cov_min': float(row.get('cov_min', 0)),
                        'cov_med': float(row.get('cov_med', 0))
                    }
                    
                    # Store under original sample name only
                    cleaning_data[sample_name] = cleaning_stats
            
        logger.info(f"Parsed cleaning data for {len(cleaning_data)} samples from CSV")
        return cleaning_data
            
    except Exception as e:
        logger.error(f"Error reading cleaning CSV file {file_path}: {e}")
        return {}
        
def parse_ref_seqs_csv(file_path):
    """Parse the reference sequences CSV file for taxid, accession, and rank information."""
    ref_data = {}
    try:
        with open(file_path, 'r', newline='') as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                process_id = row.get('process_id', '').strip()
                if process_id:
                    ref_stats = {
                        'taxid': row.get('taxid', '').strip(),
                        'protein_accession': row.get('protein_accession', '').strip(),
                        'matched_rank': row.get('matched_rank', '').strip()
                    }
                    ref_data[process_id] = ref_stats
        
        logger.info(f"Parsed reference sequence data for {len(ref_data)} samples from CSV")
        return ref_data
            
    except Exception as e:
        logger.error(f"Error reading reference sequences CSV file {file_path}: {e}")
        return {}


def process_fasta_file(file_path):
    """Process a FASTA file and extract sequence statistics."""
    logger.info(f"Processing file: {file_path}")  # Log which file is being processed
    
    if os.path.getsize(file_path) == 0:
        base_id, _, _ = extract_sample_info(os.path.basename(file_path))
        logger.warning(f"Empty file: {file_path}")
        return {
            'ID': base_id,
            'n_reads_aligned': 0,
            'cov_min': 0.00,
            'cov_max': 0.00,
            'cov_avg': 0.00,
            'cov_med': 0.00,
        }

    try:
        with open(file_path, 'r', encoding='utf-8', errors='replace') as handle:
            sequences = list(SeqIO.parse(handle, 'fasta'))
    except Exception as e:
        logger.error(f"Error reading file {file_path}: {e}")
        base_id, _, _ = extract_sample_info(os.path.basename(file_path))
        return {
            'ID': base_id,
            'n_reads_aligned': 0,
            'cov_min': 0.00,
            'cov_max': 0.00,
            'cov_avg': 0.00,
            'cov_med': 0.00,
        }

    if not sequences:
        base_id, _, _ = extract_sample_info(os.path.basename(file_path))
        logger.warning(f"No sequences found in file: {file_path}")
        return {
            'ID': base_id,
            'n_reads_aligned': 0,
            'cov_min': 0.00,
            'cov_max': 0.00,
            'cov_avg': 0.00,
            'cov_med': 0.00,
        }

    unique_sequences = {}
    for seq in sequences:
        if seq.id not in unique_sequences:
            unique_sequences[seq.id] = seq
    
    # Detect sequence length mismatches
    sequence_lengths = [len(seq.seq) for seq in unique_sequences.values()]
    if len(set(sequence_lengths)) > 1:
        logger.error(f"Error in file {file_path}: Found sequences with different lengths: {set(sequence_lengths)}")
        # Log some example sequence IDs with their lengths
        for seq_id, seq in list(unique_sequences.items())[:5]:  # Log up to 5 examples
            logger.error(f"  Sequence {seq_id}: length {len(seq.seq)}")
        
        # Use the most common length
        from collections import Counter
        most_common_length = Counter(sequence_lengths).most_common(1)[0][0]
        logger.info(f"Using most common length: {most_common_length} for file {file_path}")
        
        # Filter sequences to only include those with the most common length
        filtered_sequences = {seq_id: seq for seq_id, seq in unique_sequences.items() 
                             if len(seq.seq) == most_common_length}
        
        if not filtered_sequences:
            logger.error(f"No sequences with common length in {file_path}")
            base_id, _, _ = extract_sample_info(os.path.basename(file_path))
            return {
                'ID': base_id,
                'n_reads_aligned': len(unique_sequences),
                'cov_min': 0.00,
                'cov_max': 0.00,
                'cov_avg': 0.00,
                'cov_med': 0.00,
            }
        
        logger.info(f"Filtered from {len(unique_sequences)} to {len(filtered_sequences)} sequences")
        unique_sequences = filtered_sequences

    sequence_count = len(unique_sequences)
    
    try:
        first_seq = next(iter(unique_sequences.values()))
        coverage = np.zeros(len(first_seq.seq))
        
        for seq_id, seq in unique_sequences.items():
            try:
                seq_array = np.array([1 if base != '-' else 0 for base in seq.seq])
                if coverage.shape != seq_array.shape:
                    logger.error(f"Shape mismatch in {file_path}: Expected {coverage.shape}, got {seq_array.shape} for sequence {seq_id}")
                    # Skip this sequence
                    continue
                coverage += seq_array
            except Exception as e:
                logger.error(f"Error processing sequence {seq_id} in file {file_path}: {e}")
                # Skip this sequence and continue with others
                continue

        # Calculate coverage statistics
        if len(coverage) > 0:
            min_coverage = np.min(coverage)
            max_coverage = np.max(coverage)
            mean_coverage = np.mean(coverage)
            median_coverage = np.median(coverage)
        else:
            min_coverage = max_coverage = mean_coverage = median_coverage = 0
            
    except Exception as e:
        logger.error(f"Error calculating coverage for file {file_path}: {e}")
        base_id, _, _ = extract_sample_info(os.path.basename(file_path))
        return {
            'ID': base_id,
            'n_reads_aligned': sequence_count,
            'cov_min': 0.00,
            'cov_max': 0.00,
            'cov_avg': 0.00,
            'cov_med': 0.00,
        }

    base_id, _, _ = extract_sample_info(os.path.basename(file_path))

    return {
        'ID': base_id,
        'n_reads_aligned': sequence_count,
        'cov_min': round(min_coverage, 2),
        'cov_max': round(max_coverage, 2),
        'cov_avg': round(mean_coverage, 2),
        'cov_med': round(median_coverage, 2),
    }

def summarise_fasta(log_file, output_file, out_file_dir, cleaning_csv=None, ref_seqs_csv=None, pre_fasta=None, post_fasta=None):
    """Summarise FASTA files based on log file with cleaning statistics from CSV and auto-detect merge mode."""
    
    # Auto-detect merge mode from file paths
    merge_mode = False
    paths_to_check = [output_file, log_file, out_file_dir]
    if cleaning_csv:
        paths_to_check.append(cleaning_csv)
    if ref_seqs_csv:
        paths_to_check.append(ref_seqs_csv)
    if pre_fasta:
        paths_to_check.append(pre_fasta)
    if post_fasta:
        paths_to_check.append(post_fasta)
    
    for path in paths_to_check:
        if path and 'merge_mode' in str(path).lower():
            merge_mode = True
            logger.debug(f"Merge mode detected from path: {path}")
            break
    try:
        with open(log_file, 'r') as f:
            file_paths = [line.strip() for line in f if line.strip().endswith(('.fasta', '.fas'))]
    except Exception as e:
        logger.error(f"Error reading log file {log_file}: {e}")
        return

    # Log file information
    logger.info(f"The input log_file contains paths to {len(file_paths)} files for processing")
    
    if not file_paths:
        logger.warning(f"No valid FASTA files found in log file: {log_file}")
        return

    # Extract and log sample names
    sample_info = [extract_sample_info(os.path.basename(path)) for path in file_paths]
    base_ids = [info[0] for info in sample_info if info[0]]
    
    logger.info(f"Sample base IDs found: {', '.join(set(base_ids))}")
    logger.info(f"Total number of samples being processed: {len(file_paths)}")

    # Parse FASTA headers if provided - NEW LOGIC
    fasta_headers_lookup = {}
    if pre_fasta:
        logger.info(f"Parsing pre-fasta headers from: {pre_fasta}")
        pre_headers = parse_fasta_headers(pre_fasta)
        fasta_headers_lookup.update(pre_headers)
    
    if post_fasta:
        logger.info(f"Parsing post-fasta headers from: {post_fasta}")
        post_headers = parse_fasta_headers(post_fasta)
        # Post-fasta headers override pre-fasta if there are duplicates
        fasta_headers_lookup.update(post_headers)
    
    logger.info(f"Total FASTA headers available for matching: {len(fasta_headers_lookup)}")
    
    # Log merge mode status
    if merge_mode:
        logger.info("Merge mode auto-detected from file paths: will append '_merge' to mge_params column")
    else:
        logger.info("Merge mode not detected: no suffix modification")

    # Get alignment stats from .out files
    out_files = [f for f in os.listdir(out_file_dir) if f.endswith('.out')]
    logger.info(f"The out_file_dir contains {len(out_files)} .out files for processing")
    
    out_file_data = {}
    for file in out_files:
        _, full_id, _ = extract_sample_info(file)
        if full_id:
            out_data = parse_out_file(os.path.join(out_file_dir, file))
            if out_data:
                out_file_data[full_id] = out_data
    
    # Get cleaning stats from CSV file
    cleaning_data = {}
    if cleaning_csv and os.path.exists(cleaning_csv):
        cleaning_data = parse_cleaning_csv(cleaning_csv)
    else:
        logger.warning(f"Cleaning CSV file not provided or does not exist")

    # Get reference sequence stats from CSV file
    ref_data = {}
    if ref_seqs_csv and os.path.exists(ref_seqs_csv):
        ref_data = parse_ref_seqs_csv(ref_seqs_csv)
    else:
        logger.warning(f"Reference sequences CSV file not provided or does not exist")

    # Updated fieldnames with fasta_header column added after mge_params
    fieldnames = [
        'Filename', 'ID', 'mge_params', 'fasta_header', 'sample_taxid', 'ref_accession', 'ref_rank',
        'n_reads_in', 'n_reads_aligned', 'n_reads_skipped', 'ref_length', 
        'cov_min', 'cov_max', 'cov_avg', 'cov_med',
        'cleaning_removed_human', 'cleaning_removed_at', 'cleaning_removed_outlier', 
        'cleaning_removed_reference', 'cleaning_ambig_bases', 'cleaning_cov_percent'
    ]

    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        # Process FASTA files - create rows with alignment stats and cleaning columns set to "null"
        for file_path in file_paths:
            base_id, full_id, params = extract_sample_info(os.path.basename(file_path))
            if not base_id:
                continue

            result = process_fasta_file(file_path)
            if result:
                filename = os.path.basename(file_path).replace('.fasta', '').replace('.fas', '').replace('_align_', '_')
                result['Filename'] = filename
                result['ID'] = base_id
                ref_stats = ref_data.get(base_id, {})
                result['sample_taxid'] = ref_stats.get('taxid', 'null')
                result['ref_accession'] = ref_stats.get('protein_accession', 'null')
                result['ref_rank'] = ref_stats.get('matched_rank', 'null')
                result['mge_params'] = params if params else 'null'

                # Apply merge mode suffix if requested
                if merge_mode:
                    # Append '_merge' to mge_params if not already present
                    if params and not params.endswith('_merge'):
                        result['mge_params'] = f"{params}_merge"

                # NEW LOGIC: Look up FASTA header using (ID, mge_params) combination
                lookup_key = (base_id, result['mge_params'])
                result['fasta_header'] = fasta_headers_lookup.get(lookup_key, 'null')
                logger.debug(f"Looking up header for {lookup_key}: {result['fasta_header']}")

                # Add alignment stats with updated keys - use full_id for matching
                out_data = out_file_data.get(full_id, {})
                result['n_reads_in'] = out_data.get('n_reads_in', 'null')
                result['ref_length'] = out_data.get('ref_length', 'null')
                
                # Calculate n_reads_skipped (successful_aligned - n_reads_aligned)
                successful_aligned = out_data.get('successful_aligned', None)
                n_reads_aligned = result.get('n_reads_aligned', 0)
                
                if successful_aligned is not None and n_reads_aligned is not None:
                    result['n_reads_skipped'] = max(0, successful_aligned - n_reads_aligned)
                else:
                    result['n_reads_skipped'] = 'null'

                # Set all cleaning columns to "null" for FASTA-derived rows
                result['cleaning_removed_human'] = 'null'
                result['cleaning_removed_at'] = 'null'
                result['cleaning_removed_outlier'] = 'null'
                result['cleaning_removed_reference'] = 'null'
                result['cleaning_ambig_bases'] = 'null'
                result['cleaning_cov_percent'] = 'null'

                # Write FASTA-derived row to the CSV file
                writer.writerow(result)

        # Process cleaning CSV entries - create separate rows with cleaning stats and alignment columns set to "null"
        logger.info(f"Starting to process {len(cleaning_data)} cleaning data entries")
        
        processed_count = 0
        skipped_count = 0
        
        for sample_name, cleaning_stats in cleaning_data.items():
            logger.debug(f"Processing cleaning entry: {sample_name}")
            
            # Skip the duplicate keys we created during parsing - only process original sample names
            if not ('_fcleaner' in sample_name or '_merge' in sample_name):
                logger.debug(f"Skipping {sample_name} - no _fcleaner or _merge")
                skipped_count += 1
                continue
                
            # Extract base_id and mge_params using extract_cleaning_info function
            base_id, mge_params = extract_cleaning_info(sample_name)
            
            if not base_id:
                logger.warning(f"Could not extract base_id from cleaning sample name: {sample_name}")
                continue

            processed_count += 1

            # Create cleaning-derived result
            cleaning_result = {}
            cleaning_result['Filename'] = sample_name
            cleaning_result['ID'] = base_id
            
            # Get reference data if available
            ref_stats = ref_data.get(base_id, {})
            cleaning_result['sample_taxid'] = ref_stats.get('taxid', 'null')
            cleaning_result['ref_accession'] = ref_stats.get('protein_accession', 'null')
            cleaning_result['ref_rank'] = ref_stats.get('matched_rank', 'null')
            cleaning_result['mge_params'] = mge_params if mge_params else 'null'

            # Apply merge mode suffix if requested
            if merge_mode:
                # Append '_merge' to mge_params if not already present
                if mge_params and not mge_params.endswith('_merge'):
                    cleaning_result['mge_params'] = f"{mge_params}_merge"

            # NEW LOGIC: Look up FASTA header using (ID, mge_params) combination
            lookup_key = (base_id, cleaning_result['mge_params'])
            cleaning_result['fasta_header'] = fasta_headers_lookup.get(lookup_key, 'null')
            logger.debug(f"Looking up header for {lookup_key}: {cleaning_result['fasta_header']}")

            # Set alignment-specific columns to "null" for cleaning-derived rows
            cleaning_result['n_reads_in'] = 'null'
            cleaning_result['n_reads_skipped'] = 'null'
            cleaning_result['ref_length'] = 'null'

            # MERGED COLUMNS: Use cleaning values in the merged column names
            # cleaning_kept_reads -> n_reads_aligned
            cleaning_result['n_reads_aligned'] = cleaning_stats.get('kept_reads', 'null')
            
            # Coverage columns: cleaning_cov_* -> cov_*
            cleaning_result['cov_min'] = round(cleaning_stats.get('cov_min', 0), 2) if 'cov_min' in cleaning_stats else 'null'
            cleaning_result['cov_max'] = round(cleaning_stats.get('cov_max', 0), 2) if 'cov_max' in cleaning_stats else 'null'
            cleaning_result['cov_avg'] = round(cleaning_stats.get('cov_avg', 0), 2) if 'cov_avg' in cleaning_stats else 'null'
            cleaning_result['cov_med'] = round(cleaning_stats.get('cov_med', 0), 2) if 'cov_med' in cleaning_stats else 'null'

            # SEPARATE CLEANING COLUMNS: Keep these as separate columns
            cleaning_result['cleaning_removed_human'] = cleaning_stats.get('removed_human', 'null')
            cleaning_result['cleaning_removed_at'] = cleaning_stats.get('removed_at', 'null')
            cleaning_result['cleaning_removed_outlier'] = cleaning_stats.get('removed_outlier', 'null')
            cleaning_result['cleaning_removed_reference'] = cleaning_stats.get('removed_reference', 'null')
            cleaning_result['cleaning_ambig_bases'] = cleaning_stats.get('ambig_bases', 'null')
            cleaning_result['cleaning_cov_percent'] = round(cleaning_stats.get('cov_percent', 0), 2) if 'cov_percent' in cleaning_stats else 'null'

            # Write cleaning-derived row to the CSV file
            writer.writerow(cleaning_result)

        logger.info(f"Finished processing cleaning data: {processed_count} processed, {skipped_count} skipped")

    output_file_abs = os.path.abspath(output_file)
    logger.info(f"CSV summary created: {output_file_abs}")
    print(f"CSV summary created: '{output_file}'.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Summarise FASTA file statistics including sequence counts, coverage, and cleaning results.')
    parser.add_argument('-a', '--alignment_log', required=True, type=str, help='The log file containing the paths to FASTA files.')
    parser.add_argument('-o', '--output', required=True, type=str, help='The output CSV file name.')
    parser.add_argument('-od', '--out_file_dir', required=True, type=str, help='The directory containing .out files with additional statistics.')
    parser.add_argument('-c', '--cleaning_csv', type=str, help='Path to CSV file containing cleaning statistics')
    parser.add_argument('--ref_seqs', type=str, help='Path to CSV file containing reference sequence information (taxid, accession, rank)')
    parser.add_argument('--pre-fasta', type=str, help='Path to multi-FASTA file containing pre-cleaning sequences for header extraction')
    parser.add_argument('--post-fasta', type=str, help='Path to multi-FASTA file containing post-cleaning sequences for header extraction')

    args = parser.parse_args()

    if not os.path.isfile(args.alignment_log):
        parser.error(f"The log file '{args.alignment_log}' does not exist.")
    if not os.path.isdir(args.out_file_dir):
        parser.error(f"The directory '{args.out_file_dir}' does not exist or is not a directory.")
    if args.cleaning_csv and not os.path.isfile(args.cleaning_csv):
        parser.error(f"The cleaning CSV file '{args.cleaning_csv}' does not exist.")
    if args.ref_seqs and not os.path.isfile(args.ref_seqs):
        parser.error(f"The reference sequences CSV file '{args.ref_seqs}' does not exist.")
    if getattr(args, 'pre_fasta') and not os.path.isfile(getattr(args, 'pre_fasta')):
        parser.error(f"The pre-fasta file '{getattr(args, 'pre_fasta')}' does not exist.")
    if getattr(args, 'post_fasta') and not os.path.isfile(getattr(args, 'post_fasta')):
        parser.error(f"The post-fasta file '{getattr(args, 'post_fasta')}' does not exist.")

    summarise_fasta(args.alignment_log, args.output, args.out_file_dir, args.cleaning_csv, args.ref_seqs, getattr(args, 'pre_fasta', None), getattr(args, 'post_fasta', None))
