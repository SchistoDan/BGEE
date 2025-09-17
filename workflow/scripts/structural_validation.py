#!/usr/bin/env python3
"""
structural_validation.py

This tool processes FASTA files containing DNA sequences and applies stringent quality 
criteria to identify the best barcode sequences for each biological sample. It combines 
structural analysis (gaps, ambiguous bases, sequence length) with functional analysis 
(reading frame determination, stop codon counting, protein translation) to ensure selected 
barcodes are suitable for species identification and phylogenetic studies.

FEATURES
=============
• Multi-file batch processing with progress tracking
• HMM-based barcode region extraction using nhmmer alignment
• Quality assessment across multiple metrics (see below) and selection of the best barcode sequence
• Automated reading frame analysis and genetic code translation
• CSV reporting with all quality metrics
• Robust gap and ambiguous base handling
• Automatic N-trimming and sequence formatting

PROCESS
=============
1. Input Processing:
   - Parse multiple FASTA files
   - Extract process IDs and parameters from sequence headers
   - Validate sequence integrity and format
2. Full Sequence Analysis:
   - Calculate structural metrics (length, gaps, ambiguous bases)
   - Determine longest continuous high-quality barcode regions
3. Barcode Extraction (nhmmer-based):
   - Remove tilde characters (~) representing missing gene regions
   - Replace gap characters (-) with N to maintain structure
   - Align sequences against COI-5P HMM profile using nhmmer
   - Construct barcode sequences in HMM coordinate space
   - Trim leading/trailing Ns while preserving internal structure
4. Translation Analysis:
   - Analyse all three possible reading frames (0, 1, 2)
   - Translate sequences using specified genetic code table
   - Count stop codons in each frame
   - Select optimal reading frame with zero stop codons
5. Quality Assessment:
   - Select best representative for each process ID (see below)
6. Output Generation:
   - Write high-quality full sequences to FASTA
   - Write corresponding barcode regions to separate FASTA
   - Generate comprehensive CSV report with all metrics


SELECTION 
===================
For each unique process ID, sequences are filtered and selected as follows:

Required Filters (all must pass):
• barcode_ambiguous_bases_original == 0 (no original N's)
• stop_codons == 0 (valid protein-coding sequence)
• barcode_base_count > 200 (sufficient length for reliable identification)
• barcode_ambiguous_bases < 30% of barcode_base_count (acceptable quality after processing)

Selection Criterion:
• Among qualifying sequences: select the one with highest barcode_base_count
• This prioritizes actual nucleotide content over total length


RANKING CRITERIA
================
Barcode Quality Ranks (1-6, lower = better):
Rank 1: Perfect sequences
  • No original ambiguous bases (N's)
  • No stop codons in optimal reading frame
  • Correct reading frame identified
  • ≥500 actual nucleotide bases
Rank 2: High quality, slightly shorter
  • No original ambiguous bases
  • No stop codons
  • Correct reading frame
  • 400-499 actual nucleotide bases
Rank 3: Good quality, moderate length
  • No original ambiguous bases
  • No stop codons
  • Correct reading frame
  • 300-399 actual nucleotide bases
Rank 4: Acceptable quality, shorter
  • No original ambiguous bases
  • No stop codons
  • Correct reading frame
  • 200-299 actual nucleotide bases
Rank 5: Minimal acceptable quality
  • No original ambiguous bases
  • No stop codons
  • Correct reading frame
  • 1-199 actual nucleotide bases
Rank 6: Problematic sequences
  • Contains original ambiguous bases (regardless of other qualities)
  • May have translation issues

Full Sequence Ranks (1-3, lower = better):
Rank 1: No ambiguous bases anywhere in sequence
Rank 2: Contains some ambiguous bases
Rank 3: Other quality issues


INPUTS
================
--output-csv/-o: Path where the analysis results CSV file will be saved
--output-fasta/-of: Path where the best full sequences FASTA file will be saved
--output-barcode/-ob: Path where the best barcode sequences FASTA file will be saved
--input/-i: One or more input FASTA files to analyse (space-separated)
--hmm: Path to HMM profile file for nhmmer alignment

Optional:
--code/-c: Genetic code table for translation (default: 1 for standard genetic code)
--log-file LOG_FILE: Specify a custom path for the log file (default: creates timestamped log)
--verbose, -v: Enable detailed debug logging

OUTPUTS
================
1. CSV report: Contains detailed metrics including translation analysis for all sequences
2. Best sequences FASTA: Contains the highest quality full sequences (one per process_id)
3. Best barcode FASTA: Contains the highest quality barcode regions (one per process_id)
4. Log file: Records the analysis process, translation details, warnings, and errors

DEPENDENIES
================
- BioPython: For parsing/manipulating FASTA files and genetic code translation
- HMMER3 (nhmmer): For HMM-based barcode region alignment
- Python 3.7+

CSV COLUMNS
================
- file, seq_id, process_id, parameters: Sequence identification
- length, leading_gaps, trailing_gaps, internal_gaps: Full sequence structure
- ambiguous_bases, longest_stretch: Full sequence quality metrics
- barcode_length, barcode_ambiguous_bases_original: Original barcode metrics
- barcode_ambiguous_bases, barcode_base_count: Processed barcode metrics
- reading_frame, stop_codons: Translation analysis results
- barcode_rank, full_rank: Quality rankings
- best_sequence, selected_full_fasta, selected_barcode_fasta: Selection flags

AUTHORS 
================
Created by Dan Parsons & Ben Prioce @ NHMUK for the BGE consotrium
Version: 2.3
License: MIT
"""

import sys
import csv
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import re
import os
import argparse
from pathlib import Path
import logging
from datetime import datetime
import subprocess
import tempfile
from typing import Optional, Dict, List
from collections import defaultdict


def check_nhmmer_available():
    """
    Check if nhmmer is available in the system PATH.
    
    Returns:
        bool: True if nhmmer is available, False otherwise
    """
    try:
        result = subprocess.run(['nhmmer', '-h'], capture_output=True, text=True)
        return result.returncode == 0
    except FileNotFoundError:
        return False


def parse_hmm_length(hmm_file):
    """
    Parse the HMM profile file to extract the length (LENG field).
    
    Parameters:
        hmm_file (str): Path to HMM profile file
        
    Returns:
        int: HMM length or None if not found
    """
    try:
        with open(hmm_file, 'r') as f:
            for line in f:
                if line.startswith('LENG'):
                    fields = line.strip().split()
                    if len(fields) >= 2:
                        return int(fields[1])
                # Stop reading after the header section
                if line.startswith('//'):
                    break
        return None
    except Exception as e:
        logging.error(f"Error parsing HMM file {hmm_file}: {str(e)}")
        return None


def setup_logging(log_file=None):
    """
    Set up logging configuration.
    
    Parameters:
        log_file (str): Optional path to log file
    """
    log_format = '%(asctime)s - %(levelname)s - %(message)s'
    logging.basicConfig(
        level=logging.INFO,
        format=log_format,
        handlers=[
            logging.StreamHandler(),
            logging.FileHandler(log_file if log_file else f"fasta_analysis_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log")
        ]
    )


def validate_files(input_files, output_csv, output_fasta, hmm_file):
    """
    Validate input and output file paths.
    
    Parameters:
        input_files (list): List of input FASTA file paths
        output_csv (str): Path to output CSV file
        output_fasta (str): Path to output FASTA file
        hmm_file (str): Path to HMM profile file
        
    Returns:
        bool: True if validation passes, False otherwise
    """
    # Check input files exist
    for file in input_files:
        if not os.path.exists(file):
            logging.error(f"Input file does not exist: {file}")
            return False
    
    # Check HMM file exists
    if not os.path.exists(hmm_file):
        logging.error(f"HMM profile file does not exist: {hmm_file}")
        return False
        
    # Check output directories exist, or create them
    for output_file in [output_csv, output_fasta]:
        output_dir = os.path.dirname(output_file)
        if output_dir and not os.path.exists(output_dir):
            try:
                os.makedirs(output_dir)
                logging.info(f"Created output directory: {output_dir}")
            except OSError as e:
                logging.error(f"Cannot create output directory {output_dir}: {e}")
                return False
                
    # Check output files don't already exist, if so overwrite
    for output_file in [output_csv, output_fasta]:
        if os.path.exists(output_file):
            logging.warning(f"Output file already exists and will be overwritten: {output_file}")
            
    return True


def trim_n_characters(sequence):
    """
    Trim leading and trailing N characters from a sequence.
    
    Parameters:
        sequence (str): The sequence to trim
        
    Returns:
        str: Trimmed sequence
    """
    # First trim leading N's
    start = 0
    while start < len(sequence) and sequence[start] in ['N', 'n']:
        start += 1
    
    # If sequence is all N's, return empty string
    if start == len(sequence):
        return ""
    
    # Then trim trailing N's
    end = len(sequence) - 1
    while end >= 0 and sequence[end] in ['N', 'n']:
        end -= 1
    
    # Return trimmed sequence
    return sequence[start:end + 1]


def calculate_longest_stretch_full_seq(sequence):
    """
    Calculate longest continuous stretch without gaps AND N's for full sequences.
    
    Parameters:
        sequence (str): The sequence to analyze
        
    Returns:
        int: Length of longest stretch without gaps or N's
    """
    if not sequence:
        return 0
    
    # Split by both gap characters and N's (case insensitive)
    fragments = re.split(r'[-~Nn]', sequence)
    return max(len(fragment) for fragment in fragments) if fragments else 0


def calculate_barcode_base_count(sequence):
    """
    Calculate the number of actual nucleotide bases in the barcode sequence.
    This is the length of the sequence minus all N characters.
    
    Parameters:
        sequence (str): The barcode sequence to analyze
        
    Returns:
        int: Number of actual nucleotide bases (length - N count)
    """
    if not sequence:
        return 0
    
    # Count all N characters (case insensitive)
    n_count = sequence.upper().count('N')
    
    # Return total length minus N's
    return len(sequence) - n_count


def replace_gaps_with_n(sequence):
    """
    Replace all gap characters (-) with N characters in the sequence.
    
    Parameters:
        sequence (str): The sequence to process
        
    Returns:
        str: Sequence with gaps replaced by N's
    """
    return sequence.replace('-', 'N')


def get_complete_codons(seq, offset):
    """
    Extract complete codons from sequence starting at given offset.
    
    Parameters:
        seq (str): Input sequence
        offset (int): Reading frame offset (0, 1, or 2)
        
    Returns:
        str: String of complete codons (no gaps or Ns)
    """
    complete_codons = ""
    # Loop over the sequence in steps of 3, starting from the given offset.
    for i in range(offset, len(seq) - 2, 3):
        codon = seq[i:i + 3]
        # Only add the codon if it is complete and contains no gap characters or Ns
        if '-' not in codon and 'N' not in codon:
            complete_codons += codon
    return complete_codons
    
def run_nhmmer_on_sequence(sequence, seq_id, hmm_file):
    """
    Run nhmmer on a single sequence and return the best alignment result.
    
    Parameters:
        sequence (str): The sequence to align
        seq_id (str): Sequence identifier
        hmm_file (str): Path to HMM file
        
    Returns:
        dict: Dictionary with alignment coordinates or None if no significant match
    """
    try:
        # Create FASTA file with the sequence
        with tempfile.NamedTemporaryFile(mode='w+', suffix='.fasta', delete=False) as temp_input, \
             tempfile.NamedTemporaryFile(mode='w+', suffix='.tbl', delete=False) as temp_tblout:
            
            temp_input.write(f">{seq_id}\n{sequence}\n")
            temp_input.flush()
            
            # Run nhmmer with separate tabular output file
            nhmmer_cmd = [
                'nhmmer',
                '--tblout', temp_tblout.name,
                '--incE', '1e-3',
                '--cpu', '1',
                str(hmm_file),
                temp_input.name
            ]
            
            logging.debug(f"Running nhmmer on sequence {seq_id}: {' '.join(nhmmer_cmd)}")
            result = subprocess.run(nhmmer_cmd, capture_output=True, text=True)
            
            if result.returncode != 0:
                logging.error(f"nhmmer failed for {seq_id}: {result.stderr}")
                return None

            # Log complete nhmmer output for each sample
            logging.debug(f"=== COMPLETE nhmmer OUTPUT for {seq_id} ===")
            logging.debug(result.stdout)
            
            # Parse the clean tabular output file
            with open(temp_tblout.name, 'r') as f:
                tabular_content = f.read()
            
            logging.debug(f"=== TABULAR OUTPUT for {seq_id} ===")
            logging.debug(tabular_content)
            
            # Parse tabular output to get best alignment
            alignment_result = parse_nhmmer_result(tabular_content, seq_id)
            
            return alignment_result

    except FileNotFoundError:
        logging.error("nhmmer not found in PATH")
        return None
    except Exception as e:
        logging.error(f"Error running nhmmer on sequence {seq_id}: {str(e)}")
        return None


def parse_nhmmer_result(tabular_content, seq_id):
    """
    Parse nhmmer tabular output for the best alignment result.
    
    Parameters:
        tabular_content (str): Content of the .tbl file
        seq_id (str): Original sequence ID
        
    Returns:
        dict: Dictionary with alignment coordinates or None if no significant match
    """
    best_result = None
    best_evalue = float('inf')
    
    try:
        # Parse each line of tabular output
        for line in tabular_content.split('\n'):
            line = line.strip()
            
            # Skip comments and empty lines
            if not line or line.startswith('#'):
                continue
            
            # Split fields
            fields = line.split()
            if len(fields) >= 15:  # Ensure we have all required fields
                try:
                    target_name = fields[0]
                    hmm_from = int(fields[4])
                    hmm_to = int(fields[5])
                    seq_from = int(fields[6])
                    seq_to = int(fields[7])
                    evalue = float(fields[12])
                    score = float(fields[13])
                    bias = float(fields[14])
                    
                    logging.info(f"Found alignment for {target_name}: HMM {hmm_from}-{hmm_to}, "
                                f"Seq {seq_from}-{seq_to}, E-value: {evalue}, Score: {score}, Bias: {bias}")
                    
                    # Only include significant matches (should already be filtered by nhmmer --incE)
                    if evalue <= 1e-3 and target_name == seq_id:
                        # Keep the best (lowest E-value) result
                        if evalue < best_evalue:
                            best_evalue = evalue
                            best_result = {
                                'target_name': target_name,
                                'hmm_from': hmm_from,
                                'hmm_to': hmm_to,
                                'seq_from': seq_from,
                                'seq_to': seq_to,
                                'score': score,
                                'bias': bias,
                                'evalue': evalue
                            }
                            logging.info(f"New best alignment for {seq_id}: E-value={evalue}, Score={score}")
                    else:
                        logging.info(f"Rejected alignment for {seq_id}: E-value {evalue} > threshold 1e-3 or name mismatch")
                
                except (ValueError, IndexError) as e:
                    logging.warning(f"Could not parse tabular line: {line} ({e})")
                    continue
            else:
                logging.debug(f"Insufficient fields in tabular line: {line} (got {len(fields)}, need 15)")
        
    except Exception as e:
        logging.error(f"Error parsing nhmmer tabular output: {str(e)}")
    
    # Log final result summary
    if best_result:
        logging.debug(f"=== BEST ALIGNMENT SUMMARY for {seq_id} ===")
        logging.debug(f"E-value: {best_result['evalue']}")
        logging.debug(f"Score: {best_result['score']}")
        logging.debug(f"Bias: {best_result['bias']}")
        logging.debug(f"HMM coordinates: {best_result['hmm_from']}-{best_result['hmm_to']}")
        logging.debug(f"Sequence coordinates: {best_result['seq_from']}-{best_result['seq_to']}")
        logging.debug(f"HMM alignment length: {best_result['hmm_to'] - best_result['hmm_from'] + 1}")
        logging.debug(f"Sequence alignment length: {best_result['seq_to'] - best_result['seq_from'] + 1}")
    else:
        logging.debug(f"=== NO SIGNIFICANT ALIGNMENT FOUND for {seq_id} ===")
    
    return best_result


def construct_hmm_space_from_alignment(alignment_result, sequence, hmm_length):
    """
    Construct HMM coordinate space sequence from single alignment result.
    
    Parameters:
        alignment_result (dict): Dictionary with alignment coordinates
        sequence (str): The original N-padded sequence
        hmm_length (int): Length of the HMM profile
        
    Returns:
        str: Sequence string in HMM coordinate space
    """
    # Initialize HMM sequence with gaps
    hmm_sequence = ['-'] * hmm_length
    
    logging.debug(f"Constructing HMM space sequence from alignment")
    
    # Extract alignment coordinates
    hmm_start = alignment_result['hmm_from'] - 1  # Convert to 0-based
    hmm_end = alignment_result['hmm_to']
    seq_start = alignment_result['seq_from'] - 1  # Convert to 0-based
    seq_end = alignment_result['seq_to']
    
    logging.debug(f"Placing sequence at HMM positions: {alignment_result['hmm_from']}-{alignment_result['hmm_to']} (1-based)")
    
    # Extract the aligned portion of the sequence
    aligned_sequence_portion = sequence[seq_start:seq_end]
    
    # Place sequence at HMM coordinates
    total_coverage = 0
    seq_pos = 0
    for hmm_pos in range(hmm_start, min(hmm_end, hmm_length)):
        if seq_pos < len(aligned_sequence_portion):
            hmm_sequence[hmm_pos] = aligned_sequence_portion[seq_pos]
            total_coverage += 1
            seq_pos += 1
    
    # Convert to string
    final_sequence = ''.join(hmm_sequence)
    
    # Log coverage statistics
    coverage_percentage = (total_coverage / hmm_length) * 100
    
    logging.debug(f"=== HMM SPACE CONSTRUCTION for {alignment_result['target_name']} ===")
    logging.debug(f"Total HMM length: {hmm_length} positions")
    logging.debug(f"Covered positions: {total_coverage} ({coverage_percentage:.1f}%)")
    logging.debug(f"Missing positions: {hmm_length - total_coverage} ({100-coverage_percentage:.1f}%)")
    logging.debug(f"Final sequence length: {len(final_sequence)}")
    logging.debug(f"Final sequence preview: {final_sequence[:100]}{'...' if len(final_sequence) > 100 else ''}")
    
    return final_sequence


def trim_sequence_ends(sequence):
    """
    Trim leading and trailing Ns while preserving internal Ns.
    
    Parameters:
        sequence (str): Input sequence string
        
    Returns:
        str: Trimmed sequence string
    """
    # Remove leading Ns (gaps already converted to N's at this point)
    start = 0
    while start < len(sequence) and sequence[start] in 'N':
        start += 1
    
    # Remove trailing Ns
    end = len(sequence)
    while end > start and sequence[end-1] in 'N':
        end -= 1
    
    trimmed = sequence[start:end]
    
    logging.debug(f"=== SEQUENCE TRIMMING ===")
    logging.debug(f"Original length: {len(sequence)}")
    logging.debug(f"Trimmed length: {len(trimmed)}")
    logging.debug(f"Removed from start: {start} characters")
    logging.debug(f"Removed from end: {len(sequence) - end} characters")
    
    return trimmed


def align_sequence_with_nhmmer(record, hmm_file, hmm_length):
    """
    Process sequence using N-padding nhmmer approach.
    
    Process sequence by:
    1. Removing tilde characters (preserve gaps)
    2. Replacing gap characters with N characters
    3. Running nhmmer on the single N-padded sequence
    4. Constructing HMM space sequence from alignment results
    5. Replacing HMM coordinate space gaps with N characters
    6. Trimming leading/trailing Ns while preserving internal Ns
    
    Parameters:
        record (SeqRecord): The DNA sequence record to align
        hmm_file (str): Path to HMM profile file
        hmm_length (int): Length of the HMM profile
        
    Returns:
        tuple: (aligned_seq, original_seq_before_nhmmer) where aligned_seq is SeqRecord in HMM coordinate space 
               and original_seq_before_nhmmer is the sequence before nhmmer processing (after step 1)
    """
    try:
        # Step 1: Remove tilde characters (preserve gaps - they're biologically meaningful!)
        seq_str = str(record.seq).replace('~', '')
        logging.debug(f"After tilde removal: {seq_str[:100]}...")
        
        # Store the original sequence after tilde removal but before N-padding (for original N counting)
        original_seq_before_nhmmer = seq_str
        
        # Step 2: Replace all gap characters with N characters
        n_padded_seq = seq_str.replace('-', 'N')
        logging.debug(f"After gap-to-N replacement: {n_padded_seq[:100]}...")
        
        # Step 3: Run nhmmer on the single N-padded sequence
        nhmmer_result = run_nhmmer_on_sequence(n_padded_seq, record.id, hmm_file)
        if not nhmmer_result:
            logging.warning(f"No significant nhmmer alignment found for sequence {record.id}")
            return None, original_seq_before_nhmmer
        
        # Step 4: Construct HMM space sequence using alignment coordinates
        hmm_sequence = construct_hmm_space_from_alignment(nhmmer_result, n_padded_seq, hmm_length)
        
        # Step 5: Replace HMM coordinate space gaps with N characters
        hmm_sequence_with_n = replace_gaps_with_n(hmm_sequence)
        
        # Count gaps replaced for logging
        gap_count = hmm_sequence.count('-')
        logging.debug(f"Replaced {gap_count} HMM coordinate space gaps with N's")
        
        # Step 6: Trim leading/trailing Ns while preserving internal Ns
        trimmed_sequence = trim_sequence_ends(hmm_sequence_with_n)
        
        aligned_record = SeqRecord(
            Seq(trimmed_sequence),
            id=record.id,
            description=record.description
        )
        
        return aligned_record, original_seq_before_nhmmer

    except Exception as e:
        logging.error(f"Error in N-padding nhmmer alignment for {record.id}: {str(e)}")
        return None, str(record.seq).replace('~', '')
        
def determine_reading_frame(aligned_seq, trans_table):
    """
    Determine the reading frame from HMM-aligned sequence.

    Parameters:
        aligned_seq (str): HMM-aligned sequence string
        trans_table (int): Translation table to use
        
    Returns:
        tuple: (reading_frame, stop_codons_count, protein_sequence)
    """
    best_frame = 0
    best_protein = ""
    min_stops = float('inf')

    logging.info(f"=== Analysing reading frames ===")
    
    # Analyze all 3 reading frames
    frame_results = []
    for frame in range(3):
        coding_seq = get_complete_codons(aligned_seq, frame)
        if coding_seq:
            protein = Seq(coding_seq).translate(table=trans_table)
            stops = protein.count('*')
        else:
            protein = Seq("")
            stops = 0
        
        frame_results.append({
            'frame': frame,
            'stops': stops,
            'protein_length': len(protein),
            'coding_length': len(coding_seq)
        })
        
        logging.info(f"Frame {frame}: {stops} stop codons, {len(protein)} amino acids, {len(coding_seq)} coding nucleotides")

        if stops < min_stops:
            min_stops = stops
            best_frame = frame
            best_protein = protein

    # Log results
    logging.info(f"Best reading frame: {best_frame}")
    logging.info(f"Stop codons in best frame: {min_stops}")
    
    # Log translation details
    logging.info(f"=== Translating sequence ===")
    logging.info(f"Protein length: {len(best_protein)} amino acids")
    logging.info(f"Amino acid sequence: {str(best_protein)}")
    
    # Log stop codon analysis
    logging.info(f"=== Analysing stop codons ===")
    logging.info(f"Translation table: {trans_table}")
    logging.info(f"Reading frame: {best_frame}")
    logging.info(f"Searching for stop codons in protein sequence (length: {len(best_protein)})")
    logging.info(f"Total stop codons found: {min_stops}")
    
    # Find stop codon positions (excluding natural terminal stop)
    stop_positions = []
    for i, aa in enumerate(best_protein[:-1]):  # Exclude last amino acid (natural terminal)
        if aa == '*':
            nuc_pos = (i * 3) + best_frame
            stop_positions.append(nuc_pos)
    
    logging.info(f"Stop codon positions (nucleotide coordinates): {stop_positions}")
    
    return best_frame, min_stops, str(best_protein)


def calculate_barcode_rank(barcode_ambiguous_bases_original, stop_codons, reading_frame_valid, barcode_base_count):
    """
    Calculate barcode rank based on updated criteria.
    
    Parameters:
        barcode_ambiguous_bases_original (int): Original N's before nhmmer processing
        stop_codons (int): Number of stop codons in best reading frame
        reading_frame_valid (bool): True if reading frame has zero stop codons
        barcode_base_count (int): Number of actual nucleotide bases in barcode
        
    Returns:
        int: Barcode rank (1-6)
    """
    # Rank 6: Has original N's
    if barcode_ambiguous_bases_original > 0:
        return 6
    
    # For ranks 1-5, must have no original N's, no stop codons, and correct reading frame
    if stop_codons == 0 and reading_frame_valid:
        if barcode_base_count >= 500:
            return 1
        elif 400 <= barcode_base_count <= 499:
            return 2
        elif 300 <= barcode_base_count <= 399:
            return 3
        elif 200 <= barcode_base_count <= 299:
            return 4
        elif 1 <= barcode_base_count <= 199:
            return 5
    
    # If conditions not met, default to worst rank
    return 6


def calculate_full_rank(ambiguous_bases):
    """Calculate full rank based on criteria."""
    if ambiguous_bases == 0:
        return 1
    elif ambiguous_bases >= 1:
        return 2
    return 3


def passes_quality_criteria(result):
    """
    Check if a sequence result passes the updated quality criteria.
    
    Parameters:
        result (dict): Sequence result dictionary
        
    Returns:
        bool: True if sequence meets all quality criteria
    """
    # Basic structural validation
    if result['barcode_ambiguous_bases_original'] != 0 or result['stop_codons'] != 0:
        return False
    
    # Additional quality criteria
    if result['barcode_base_count'] <= 200:
        return False
    
    # Check if barcode_ambiguous_bases < 30% of barcode_base_count
    if result['barcode_base_count'] > 0:
        ambiguous_percentage = result['barcode_ambiguous_bases'] / result['barcode_base_count']
        if ambiguous_percentage >= 0.30:
            return False
    
    return True


def select_sequences_for_process(results):
    """
    Select the best sequence for a process based on updated quality criteria.
    
    Filter: 
    1. barcode_ambiguous_bases_original == 0 AND stop_codons == 0 (REQUIRED)
    2. barcode_base_count > 200
    3. barcode_ambiguous_bases < 30% of barcode_base_count
    Select: Highest barcode_base_count (actual nucleotide content)
    
    Parameters:
        results (list): List of sequence results for a process
        
    Returns:
        tuple: (best_sequence, best_sequence) or (None, None) if no qualifying sequences
    """
    # Filter sequences that meet all quality criteria
    qualifying_sequences = [
        result for result in results
        if passes_quality_criteria(result)
    ]
    
    if not qualifying_sequences:
        return None, None
    
    # Select sequence with highest barcode_base_count
    best_sequence = max(qualifying_sequences, key=lambda x: x['barcode_base_count'])
    
    return best_sequence, best_sequence


def format_sequence(seq_record, trim_gaps=True, convert_internal_gaps=True):
    """
    Format a sequence according to specified criteria.
    
    Parameters:
        seq_record (SeqRecord): The sequence record to format
        trim_gaps (bool): Whether to trim leading and trailing gaps
        convert_internal_gaps (bool): Whether to convert internal gaps to N
        
    Returns:
        SeqRecord: Formatted sequence record
    """
    sequence = str(seq_record.seq)
    
    if trim_gaps:
        # Trim leading and trailing gaps
        sequence = sequence.strip('-').strip('~')
    
    if convert_internal_gaps:
        # First remove ~ characters and stitch sequence
        sequence = sequence.replace('~', '')
        # Then replace remaining gaps (-) with N
        sequence = re.sub(r'-', 'N', sequence)
    
    # Trim leading and trailing N's
    sequence = trim_n_characters(sequence)
    
    # Create new sequence record with formatted sequence
    new_record = SeqRecord(
        Seq(sequence),
        id=seq_record.id,
        description=seq_record.description
    )
    
    return new_record


def format_barcode_sequence(aligned_barcode_record):
    """
    Format the nhmmer-aligned barcode sequence.
    
    Parameters:
        aligned_barcode_record (SeqRecord): The nhmmer-aligned barcode sequence record
        
    Returns:
        SeqRecord: Formatted barcode sequence record
    """
    if not aligned_barcode_record:
        return None
        
    # The sequence is already processed with gaps replaced by N's and trimmed
    sequence = str(aligned_barcode_record.seq)
    
    # Format gaps according to rules for barcode output
    formatted_sequence = format_barcode_gaps(sequence)
    
    # Trim leading and trailing N's
    formatted_sequence = trim_n_characters(formatted_sequence)
    
    new_record = SeqRecord(
        Seq(formatted_sequence),
        id=aligned_barcode_record.id,
        description=aligned_barcode_record.description
    )
    
    return new_record


def format_barcode_gaps(sequence):
    """
    Format gaps in barcode sequence according to specified rules:
    - Remove internal ~ characters and stitch sequence back together
    - For gaps marked with -, if ≤ 6 bases fill with N
    - For gaps > 6 bases marked with -, keep the longest fragment
    
    Note: In the new workflow, gaps are already replaced with N's during HMM processing,
    so this function primarily handles any remaining formatting needs.
    
    Parameters:
        sequence (str): The sequence to process
        
    Returns:
        str: Formatted sequence
    """
    # First handle ~ characters by removing them and stitching sequence
    sequence = sequence.replace('~', '')
    
    # Since gaps are already replaced with N's in the new workflow,
    # we primarily need to handle any remaining gap characters
    # This preserves the original logic for backwards compatibility
    fragments = []
    current_fragment = []
    gap_count = 0
    
    for char in sequence:
        if char == '-':
            gap_count += 1
            if gap_count > 6:
                if current_fragment:
                    fragments.append(''.join(current_fragment))
                current_fragment = []
                gap_count = 0
        else:
            if gap_count > 0 and gap_count <= 6:
                # Fill small gaps with N
                current_fragment.extend(['N'] * gap_count)
            gap_count = 0
            current_fragment.append(char)
    
    # Add last fragment if exists
    if current_fragment:
        fragments.append(''.join(current_fragment))
    
    # Find and get longest fragment if we have fragments
    if not fragments:
        return sequence  # Return original if no fragments (shouldn't happen in new workflow)
    
    longest_fragment = max(fragments, key=len)
    return longest_fragment


def format_sequence_id(process_id, parameters):
    """
    Format sequence ID according to specified format.
    
    Parameters:
        process_id (str): Process ID
        parameters (str): Parameters
        
    Returns:
        str: Formatted sequence ID
    """
    return f"{process_id}_{parameters}" if parameters else process_id
    
def analyse_fasta(file_path, hmm_file, hmm_length, trans_table):
    """
    Analyse a FASTA file and extract various metrics from each sequence using nhmmer.

    Parameters:
        file_path (str): The path to the FASTA file to be analysed.
        hmm_file (str): Path to HMM profile file
        hmm_length (int): Length of the HMM profile
        trans_table (int): Translation table for genetic code

    Returns:
        dict: A dictionary containing analysis results for each sequence in the FASTA file.
    """
    try:
        # Initialise dictionaries
        results = {}
        sequences = {}
        
        # Count total sequences for progress reporting
        total_sequences = sum(1 for _ in SeqIO.parse(file_path, "fasta"))
        processed = 0

        # Parse the FASTA file
        for record in SeqIO.parse(file_path, "fasta"):
            try:
                processed += 1
                if processed % 100 == 0:  # Log progress every 100 sequences
                    logging.info(f"Processing {file_path}: {processed}/{total_sequences} sequences")

                # Get the sequence ID (header) and the sequence itself
                seq_id = record.id
                seq = str(record.seq)
                
                # Log start of processing for this sequence
                logging.info(f"")
                logging.info(f"{'='*60}")
                logging.info(f"PROCESSING SEQUENCE: {seq_id}")
                logging.info(f"{'='*60}")
                
                # Store the complete record
                sequences[seq_id] = record

                # Clean and split seq_id into process_id and parameters
                seq_id_clean = seq_id.replace("Consensus_", "") if seq_id.startswith("Consensus_") else seq_id

                # Initialize process_id and parameters
                process_id = seq_id_clean
                parameters = ""

                # Extract process_id and parameters
                if '_r_' in seq_id_clean:
                    param_start_idx = seq_id_clean.index('_r_')
                    process_id = seq_id_clean[:param_start_idx]
                    parameters = seq_id_clean[param_start_idx + 1:]

                logging.debug(f"Detected process_id: {process_id}, parameters: {parameters}")
                
                # Skip empty sequences
                if not seq or seq.strip() == "":
                    logging.warning(f"Skipping empty sequence: {seq_id}")
                    # Still add to results with null/zero values for barcode metrics
                    results[seq_id] = {
                        'file': file_path,
                        'length': 0,
                        'leading_gaps': 0,
                        'trailing_gaps': 0,
                        'internal_gaps': 0,
                        'ambiguous_bases': 0,
                        'longest_stretch': 0,
                        'barcode_length': 0,
                        'barcode_ambiguous_bases_original': 0,
                        'barcode_ambiguous_bases': 0,
                        'barcode_base_count': 0,
                        'reading_frame': -1,  # Invalid reading frame
                        'stop_codons': 0,
                        'barcode_rank': 6,  # Worst rank for empty sequences
                        'full_rank': 3,     # Worst rank for empty sequences
                        'process_id': process_id,  # Use parsed process_id
                        'parameters': parameters,  # Use parsed parameters
                        'sequence_record': record,
                        'aligned_barcode_record': None
                    }
                    continue

                # Calculate full sequence metrics
                length = len(seq)
                leading_gaps = len(seq) - len(seq.lstrip('-'))
                trailing_gaps = len(seq) - len(seq.rstrip('-'))
                internal_gaps = seq.count('-') + seq.count('~') - leading_gaps - trailing_gaps
                ambiguous_bases = seq.count('N')

                # Calculate longest stretch without gaps AND N's (full seq)
                longest_stretch = calculate_longest_stretch_full_seq(seq)

                # Use nhmmer to extract and align barcode region
                aligned_barcode, original_seq_before_nhmmer = align_sequence_with_nhmmer(record, hmm_file, hmm_length)
                
                # Initialize translation variables
                reading_frame = -1
                stop_codons = 0
                reading_frame_valid = False
                
                if aligned_barcode:
                    # Calculate barcode metrics from aligned sequence
                    barcode_seq = str(aligned_barcode.seq)
                    barcode_length = len(barcode_seq)
                    
                    # Count original N's in the sequence before nhmmer processing (after tilde removal)
                    barcode_ambiguous_bases_original = original_seq_before_nhmmer.count('N')
                    
                    # Count N's in final barcode sequence
                    barcode_ambiguous_bases = barcode_seq.count('N')
                    
                    # Calculate actual nucleotide base count (length - N count)
                    barcode_base_count = calculate_barcode_base_count(barcode_seq)
                    
                    # Count gaps in barcode for logging
                    barcode_gaps = barcode_seq.count('-')  # Should be 0 in new workflow
                    
                    # Determine reading frame and analyze translation
                    if barcode_seq:  # Only if we have sequence to translate
                        reading_frame, stop_codons, protein_seq = determine_reading_frame(barcode_seq, trans_table)
                        reading_frame_valid = (stop_codons == 0)
                else:
                    # No alignment found - set barcode metrics to None/0
                    barcode_seq = ""
                    barcode_length = 0
                    barcode_ambiguous_bases_original = 0
                    barcode_ambiguous_bases = 0
                    barcode_base_count = 0
                    barcode_gaps = 0

                # Determine ranks using new criteria
                barcode_rank = calculate_barcode_rank(
                    barcode_ambiguous_bases_original, 
                    stop_codons, 
                    reading_frame_valid, 
                    barcode_base_count
                )
                full_rank = calculate_full_rank(ambiguous_bases)

                # Log final sequence analysis summary
                logging.info(f"=== SEQUENCE ANALYSIS SUMMARY for {seq_id} ===")
                logging.info(f"Process ID: {process_id}")
                logging.info(f"Full sequence length: {length}")
                logging.info(f"Full sequence ambiguous bases: {ambiguous_bases}")
                logging.info(f"Full sequence longest stretch (no gaps/N's): {longest_stretch}")
                logging.info(f"Full sequence rank: {full_rank}")
                logging.info(f"Barcode length: {barcode_length}")
                logging.info(f"Barcode ambiguous bases (original): {barcode_ambiguous_bases_original}")
                logging.info(f"Barcode ambiguous bases (final): {barcode_ambiguous_bases}")
                if aligned_barcode:
                    logging.info(f"Barcode gaps after HMM extraction: {barcode_gaps}")
                logging.info(f"Barcode base count (nucleotides): {barcode_base_count}")
                logging.info(f"Reading frame: {reading_frame}")
                logging.info(f"Stop codons: {stop_codons}")
                logging.info(f"Reading frame valid: {reading_frame_valid}")
                logging.info(f"Barcode rank: {barcode_rank}")
                logging.info(f"Alignment successful: {'Yes' if aligned_barcode else 'No'}")
                logging.info("=" * 50)

                # Store results
                results[seq_id] = {
                    'file': file_path,
                    'length': length,
                    'leading_gaps': leading_gaps,
                    'trailing_gaps': trailing_gaps,
                    'internal_gaps': internal_gaps,
                    'ambiguous_bases': ambiguous_bases,
                    'longest_stretch': longest_stretch,
                    'barcode_length': barcode_length,
                    'barcode_ambiguous_bases_original': barcode_ambiguous_bases_original,
                    'barcode_ambiguous_bases': barcode_ambiguous_bases,
                    'barcode_base_count': barcode_base_count,
                    'reading_frame': reading_frame,
                    'stop_codons': stop_codons,
                    'barcode_rank': barcode_rank,
                    'full_rank': full_rank,
                    'process_id': process_id,
                    'parameters': parameters,
                    'sequence_record': sequences[seq_id],
                    'aligned_barcode_record': aligned_barcode
                }

            except Exception as e:
                logging.error(f"Error processing sequence {seq_id} in {file_path}: {str(e)}")
                continue

        return results

    except Exception as e:
        logging.error(f"Error analysing file {file_path}: {str(e)}")
        return {}
        
def write_best_sequences(best_sequences, output_fasta, output_barcode_fasta):
    """
    Write the best sequences that meet specific quality criteria to new FASTA files.
    Sequences are written without line breaks.
    
    Parameters:
        best_sequences (dict): Dictionary containing the best sequences for each process_id
        output_fasta (str): Path to output FASTA file for full sequences
        output_barcode_fasta (str): Path to output FASTA file for barcode regions
        
    Returns:
        dict: Selection records tracking which sequences were written to output files
    """
    try:
        selected_full_records = []
        selected_barcode_records = []
        selection_records = {}
        
        for process_id, best_result in best_sequences.items():
            # Create unique identifier for the sequence
            unique_key = f"{best_result['file']}_{best_result['seq_id']}"
            
            # Process full sequence
            seq_record = best_result['sequence_record']
            new_id = format_sequence_id(best_result['process_id'], best_result['parameters'])
            
            full_seq_record = format_sequence(seq_record)
            full_seq_record.id = new_id
            full_seq_record.description = ""
            selected_full_records.append(full_seq_record)
            
            selection_records[process_id] = {
                'full_selected_seq': unique_key,
                'barcode_selected_seq': None
            }
            
            # Process barcode sequence, if we have one
            if best_result['aligned_barcode_record']:
                aligned_barcode_record = best_result['aligned_barcode_record']
                
                barcode_seq_record = format_barcode_sequence(aligned_barcode_record)
                if barcode_seq_record:
                    barcode_seq_record.id = new_id
                    barcode_seq_record.description = ""
                    selected_barcode_records.append(barcode_seq_record)
                    selection_records[process_id]['barcode_selected_seq'] = unique_key
            
            logging.info(f"Process {process_id}:")
            logging.info(f"  Selected sequence {best_result['seq_id']} from {best_result['file']} "
                        f"(full_rank: {best_result['full_rank']}, barcode_rank: {best_result['barcode_rank']}, "
                        f"base_count: {best_result['barcode_base_count']})")

        # Write selected sequences to output files without line breaks
        def write_fasta_no_wrap(records, filename):
            with open(filename, 'w') as f:
                for record in records:
                    f.write(f">{record.id}\n{str(record.seq)}\n")
                    
        if selected_full_records:
            write_fasta_no_wrap(selected_full_records, output_fasta)
            logging.info(f"Wrote {len(selected_full_records)} sequences to {output_fasta}")
        else:
            logging.warning("No sequences met the criteria for full sequence output")
            
        if selected_barcode_records:
            write_fasta_no_wrap(selected_barcode_records, output_barcode_fasta)
            logging.info(f"Wrote {len(selected_barcode_records)} sequences to {output_barcode_fasta}")
        else:
            logging.warning("No sequences met the criteria for barcode sequence output")
            
        return selection_records
            
    except Exception as e:
        logging.error(f"Error writing sequences: {str(e)}")
        raise


def group_sequences_by_process(all_results):
    """
    Group sequences by process_id for selection.
    
    Parameters:
        all_results (list): List of all sequence results
        
    Returns:
        dict: Dictionary with process_id as key and list of results as value
    """
    process_groups = defaultdict(list)
    
    for result in all_results:
        process_id = result['process_id']
        process_groups[process_id].append(result)
    
    return dict(process_groups)


def select_best_sequences_for_all_processes(process_groups):
    """
    Select the best sequence for each process using the centralized selection logic.
    
    Parameters:
        process_groups (dict): Dictionary with process_id as key and list of results as value
        
    Returns:
        dict: Dictionary with process_id as key and best result as value
    """
    best_sequences = {}
    processes_with_no_qualifying_sequences = []
    
    for process_id, results in process_groups.items():
        best_full, best_barcode = select_sequences_for_process(results)
        
        if best_full:
            # Store the best sequence (they should be the same since we return the same sequence twice)
            best_sequences[process_id] = best_full
            logging.info(f"Process {process_id}: Selected sequence with base_count {best_full['barcode_base_count']}")
        else:
            # Log warning for processes with no qualifying sequences
            processes_with_no_qualifying_sequences.append(process_id)
            logging.warning(f"Process {process_id}: No sequences passed quality criteria "
                          f"(requires: no original N's, no stop codons, >200 base count, <30% ambiguous bases)")
    
    if processes_with_no_qualifying_sequences:
        logging.warning(f"Total processes with no qualifying sequences: {len(processes_with_no_qualifying_sequences)}")
        logging.info(f"Processes with no qualifying sequences: {', '.join(processes_with_no_qualifying_sequences)}")
    
    return best_sequences


def annotate_results_with_selection(all_results, best_sequences, selection_records):
    """
    Annotate all results with selection flags.
    
    Parameters:
        all_results (list): List of all sequence results
        best_sequences (dict): Dictionary of best sequences by process_id
        selection_records (dict): Dictionary of selection records from write_best_sequences
    """
    for result in all_results:
        process_id = result['process_id']
        unique_key = f"{result['file']}_{result['seq_id']}"
        
        # Mark as best sequence if it matches the selected sequence for this process
        result['best_sequence'] = 'yes' if (
            process_id in best_sequences and 
            result == best_sequences[process_id]
        ) else 'no'
        
        # Mark selection flags based on what was actually written to output files
        if process_id in selection_records:
            result['selected_full_fasta'] = 'yes' if selection_records[process_id]['full_selected_seq'] == unique_key else 'no'
            result['selected_barcode_fasta'] = 'yes' if selection_records[process_id]['barcode_selected_seq'] == unique_key else 'no'
        else:
            result['selected_full_fasta'] = 'no'
            result['selected_barcode_fasta'] = 'no'


def clean_results_for_csv_output(all_results):
    """
    Remove sequence record objects from results before CSV output.
    
    Parameters:
        all_results (list): List of all sequence results
    """
    for result in all_results:
        # Remove sequence records after FASTA files are written
        result.pop('sequence_record', None)
        result.pop('aligned_barcode_record', None)
        
def main():
    # Arg parser
    parser = argparse.ArgumentParser(description='Analyse FASTA files and select best COI-5P sequences using nhmmer-based barcode extraction.')
    
    # Required arguments with flags
    parser.add_argument('--output-csv', '-o', required=True, help='Path to output CSV file')
    parser.add_argument('--output-fasta', '-of', required=True, help='Path to output FASTA file for best sequences')
    parser.add_argument('--output-barcode', '-ob', required=True, help='Path to output FASTA file for barcode regions')
    parser.add_argument('--input', '-i', required=True, nargs='+', help='Input FASTA files to analyse')
    parser.add_argument('--hmm', required=True, help='Path to HMM profile file for nhmmer alignment')
    
    # Optional arguments
    parser.add_argument('--code', '-c', type=int, default=1, help='Genetic code table for translation (default: 1 for standard code)')
    parser.add_argument('--log-file', help='Path to log file (optional)')
    parser.add_argument('--verbose', '-v', action='store_true', help='Enable verbose output')
    
    args = parser.parse_args()
    
    # Setup logging FIRST
    setup_logging(args.log_file)
    
    # Set verbosity level
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    # Check nhmmer availability
    if not check_nhmmer_available():
        logging.error("nhmmer (HMMER3) is not available in the system PATH. Please install HMMER3 and ensure nhmmer is accessible.")
        sys.exit(1)
    
    # Parse HMM length
    hmm_length = parse_hmm_length(args.hmm)
    if not hmm_length:
        logging.error(f"Could not parse HMM length from profile file: {args.hmm}")
        sys.exit(1)
    
    # Log configuration info
    logging.info(f"Using HMM profile: {args.hmm} (length: {hmm_length})")
    logging.info(f"Using updated quality-based selection criteria:")
    logging.info(f"  - No original N's AND no stop codons (required)")
    logging.info(f"  - Barcode base count > 200")
    logging.info(f"  - Ambiguous bases < 30% of base count")
    logging.info(f"  - Select highest base count among qualifying sequences")
    logging.info(f"Using genetic code table: {args.code}")
        
    # Validate input and output files
    if not validate_files(args.input, args.output_csv, args.output_fasta, args.hmm):
        sys.exit(1)
    
    try:
        # Initialise results list
        all_results = []
        
        # Analyse each FASTA file
        for file in args.input:
            logging.info(f"Processing file: {file}")
            results = analyse_fasta(file, args.hmm, hmm_length, args.code)
            for seq_id, result in results.items():
                result['seq_id'] = seq_id
                all_results.append(result)

        # Group sequences by process_id
        logging.info("Grouping sequences by process_id for selection...")
        process_groups = group_sequences_by_process(all_results)
        logging.info(f"Found {len(process_groups)} unique process_ids")

        # Select best sequences for each process using centralized logic
        logging.info("Selecting best sequences for each process...")
        best_sequences = select_best_sequences_for_all_processes(process_groups)
        
        # Write best sequences to FASTA files
        logging.info("Writing selected sequences to output files...")
        selection_records = write_best_sequences(best_sequences, args.output_fasta, args.output_barcode)

        # Annotate all results with selection flags
        logging.info("Annotating results with selection flags...")
        annotate_results_with_selection(all_results, best_sequences, selection_records)

        # Clean results for CSV output
        clean_results_for_csv_output(all_results)

        # Write CSV with results
        logging.info("Writing CSV results...")
        fieldnames = [
            'file', 'seq_id', 'process_id', 'parameters', 'length',
            'leading_gaps', 'trailing_gaps', 'internal_gaps',
            'ambiguous_bases', 'longest_stretch', 'barcode_length',
            'barcode_ambiguous_bases_original', 'barcode_ambiguous_bases', 'barcode_base_count',
            'reading_frame', 'stop_codons', 'barcode_rank', 'full_rank', 'best_sequence',
            'selected_full_fasta', 'selected_barcode_fasta'
        ]

        with open(args.output_csv, 'w', newline='') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(all_results)

        logging.info(f"Analysis results saved to {args.output_csv}")
        
        # Log selection summary
        total_processes = len(process_groups)
        selected_processes = len(best_sequences)
        total_sequences = len(all_results)
        
        # Count successful barcode extractions
        successful_barcode_extractions = sum(1 for process_id, records in selection_records.items() 
                                           if records['barcode_selected_seq'] is not None)
        
        logging.info("="*60)
        logging.info("FINAL SELECTION SUMMARY")
        logging.info("="*60)
        logging.info(f"Total sequences processed: {total_sequences}")
        logging.info(f"Total unique process_ids: {total_processes}")
        logging.info(f"Processes with qualifying sequences: {selected_processes}")
        logging.info(f"Processes with no qualifying sequences: {total_processes - selected_processes}")
        logging.info(f"Full sequences written to output: {len(selection_records)}")
        logging.info(f"Barcode sequences written to output: {successful_barcode_extractions}")
        logging.info("="*60)
        
        # Log quality criteria summary
        if total_processes > selected_processes:
            logging.warning(f"Note: {total_processes - selected_processes} process(es) had no sequences meeting the quality criteria:")
            logging.warning("  - No original N's (barcode_ambiguous_bases_original == 0)")
            logging.warning("  - No stop codons (stop_codons == 0)")
            logging.warning("  - Barcode base count > 200")
            logging.warning("  - Ambiguous bases < 30% of base count")
     
    except Exception as e:
        logging.error(f"An error occurred during execution: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main()
    
