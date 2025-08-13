#!/usr/bin/env python3
"""
AT Content FASTA Sequence Filter

A tool for filtering FASTA sequence alignments based on AT content.
Extracts sequences that have AT content within specified thresholds compared to consensus.

Features:
---------
- AT content calculation using vectorised operations
- Consensus sequence generation with configurable threshold
- Three filtering modes: absolute, higher, lower
- Parallel processing for multiple files with memory-aware batching
- Comprehensive metrics and logging
- Memory-efficient processing with automatic chunking

Filtering Modes:
---------------
- absolute: Remove sequences if AT content differs from consensus by more than threshold in either direction
- higher: Remove only sequences with AT content above consensus + threshold  
- lower: Remove only sequences with AT content below consensus - threshold

Output Structure:
----------------
The tool creates the following output structure in the specified output directory:

output_dir/
├── at_filtered_sequences/          # Individual filtered FASTA files
│   ├── sample1_at_filtered.fasta
│   ├── sample2_at_filtered.fasta
│   └── ...
├── metrics/                       # Individual sample metrics
│   ├── sample1_at_metrics.csv
│   ├── sample2_at_metrics.csv
│   └── ...
├── logs/                          # Processing logs
│   ├── sample1_log.txt
│   ├── sample2_log.txt
│   ├── at_filter_main.log         # Main processing log
│   └── ...
├── at_filtered_paths.log          # Paths to all filtered files (one per line) [customisable name]
└── at_filter_summary.csv          # Summary statistics

Key Output Files:
----------------
- {filtered_files_list}: Contains full paths to all successfully filtered alignment files,
  one per line. Use this file as input to downstream tools.
- at_filter_summary.csv: Sample-level summary statistics

Note: The filtered files list name can be customised using --filtered-files-list argument.

Authors: B. Price & D. Parsons @ NHMUK (AT Content Filter Version)
Version: 1.2.0
Licence: MIT
"""

# =============================================================================
# Imports and Setup
# =============================================================================
import os
import sys
import csv
import argparse
from datetime import datetime
from typing import List, Dict, Tuple, Optional, Any, NamedTuple
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO
from Bio.Align import AlignInfo
import numpy as np
from collections import Counter, defaultdict
import traceback
import multiprocessing as mp
from functools import partial
import concurrent.futures
import gc
import time
import tarfile
import shutil

# =============================================================================
# Core Constants and Configuration
# =============================================================================
# File size threshold for chunking optimisation (100 MB)
CHUNK_THRESHOLD = 100 * 1024 * 1024  # 100 MB in bytes

# File size threshold for memory-aware processing (350 MB)
LARGE_FILE_THRESHOLD = 350 * 1024 * 1024  # 350 MB in bytes

# Process large files sequentially (1 at a time) to prevent memory exhaustion
SEQUENTIAL_LARGE_FILES = True

# =============================================================================
# File Size Classification
# =============================================================================
class FileInfo(NamedTuple):
    """Container for file information"""
    path: str
    size_bytes: int
    size_str: str
    basename: str

def classify_files_by_size(file_paths: List[str]) -> Tuple[List[FileInfo], List[FileInfo], List[str]]:
    """
    Classify files into large and normal based on size threshold.
    
    Returns:
        Tuple[List[FileInfo], List[FileInfo], List[str]]: 
        (large_files, normal_files, empty_files)
    """
    large_files = []
    normal_files = []
    empty_files = []
    
    print(f"\nClassifying {len(file_paths)} files by size...")
    
    for i, file_path in enumerate(file_paths):
        if i % 100 == 0:  # Progress update every 100 files
            print(f"Classifying file {i+1}/{len(file_paths)}")
        
        try:
            # Check if file is empty first
            is_empty, empty_reason = is_file_empty(file_path)
            if is_empty:
                empty_files.append(os.path.basename(file_path))
                continue
            
            # Get file size
            size_bytes = os.path.getsize(file_path)
            size_str = format_file_size(size_bytes)
            basename = os.path.basename(file_path)
            
            file_info = FileInfo(
                path=file_path,
                size_bytes=size_bytes,
                size_str=size_str,
                basename=basename
            )
            
            # Classify by size
            if size_bytes >= LARGE_FILE_THRESHOLD:
                large_files.append(file_info)
            else:
                normal_files.append(file_info)
                
        except OSError as e:
            print(f"Warning: Could not access file {file_path}: {e}")
            empty_files.append(os.path.basename(file_path))
    
    # Sort large files by size (largest first for better load balancing)
    large_files.sort(key=lambda x: x.size_bytes, reverse=True)
    
    print(f"\nFile classification complete:")
    print(f"- Large files (≥{format_file_size(LARGE_FILE_THRESHOLD)}): {len(large_files)}")
    print(f"- Normal files (<{format_file_size(LARGE_FILE_THRESHOLD)}): {len(normal_files)}")
    print(f"- Empty/invalid files: {len(empty_files)}")
    
    if large_files:
        print(f"\nLarge files to process:")
        for file_info in large_files[:10]:  # Show first 10
            print(f"  {file_info.basename} ({file_info.size_str})")
        if len(large_files) > 10:
            print(f"  ... and {len(large_files) - 10} more large files")
    
    return large_files, normal_files, empty_files

# =============================================================================
# Utility Functions
# =============================================================================
def memory_cleanup():
    """Simplified memory cleanup function."""
    gc.collect()
    try:
        import numpy as np
        if hasattr(np, 'clear_cache'):
            np.clear_cache()
    except:
        pass

def sequences_to_numpy_array(sequences: List[str]) -> np.ndarray:
    """Convert sequences to numpy array for vectorised operations."""
    if not sequences:
        return np.array([])
    
    max_len = max(len(seq) for seq in sequences)
    padded_sequences = [seq.ljust(max_len, '-') for seq in sequences]
    array = np.array([list(seq.upper()) for seq in padded_sequences])
    
    # Clean up intermediate data
    del padded_sequences
    gc.collect()
    
    return array

def get_base_filename(filepath: str) -> str:
    """Extract base filename without extension and without '_align' or '_human_filtered' substrings."""
    basename = os.path.basename(filepath)
    for ext in ['.fasta', '.fas', '.fa']:
        if basename.lower().endswith(ext):
            basename = basename[:-len(ext)]
            break
    
    if '_align' in basename:
        basename = basename.replace('_align', '')
    
    if '_human_filtered' in basename:
        basename = basename.replace('_human_filtered', '')
    
    return basename

def get_default_threads() -> int:
    """Get default number of threads from SLURM environment or system CPU count."""
    slurm_cpus = os.environ.get('SLURM_CPUS_PER_TASK')
    if slurm_cpus:
        try:
            return int(slurm_cpus)
        except ValueError:
            pass
    return mp.cpu_count()

def create_output_subdirectories(base_output_dir: str) -> Dict[str, str]:
    """Create subdirectories for different output file types."""
    subdirs = {
        'filtered_seqs': os.path.join(base_output_dir, 'at_filtered_sequences'),
        'metrics': os.path.join(base_output_dir, 'metrics'),
        'logs': os.path.join(base_output_dir, 'logs')
    }
    
    for subdir in subdirs.values():
        os.makedirs(subdir, exist_ok=True)
    
    return subdirs

def format_file_size(size_bytes):
    """Convert bytes to human-readable format."""
    if size_bytes == 0:
        return "0 B"
    
    size_names = ["B", "KB", "MB", "GB"]
    import math
    i = int(math.floor(math.log(size_bytes, 1024)))
    p = math.pow(1024, i)
    s = round(size_bytes / p, 2)
    return f"{s} {size_names[i]}"

def log_message(message: str, log_file=None, stdout=False, error=False) -> None:
    """Log message to file and optionally to stdout."""
    if log_file:
        log_file.write(message + '\n')
        log_file.flush()
    
    if stdout or error:
        if error:
            print(f"\nERROR: {message}", flush=True)
        else:
            print(message, flush=True)

def read_input_files_from_log(log_file_path: str) -> List[str]:
    """Read FASTA file paths from input log file."""
    file_paths = []
    
    try:
        with open(log_file_path, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if line:  # Skip empty lines
                    if not os.path.exists(line):
                        print(f"Warning: File not found (line {line_num}): {line}")
                        continue
                    file_paths.append(line)
        
        if not file_paths:
            raise ValueError("No valid file paths found in log file")
            
    except Exception as e:
        raise ValueError(f"Error reading input files log: {str(e)}")
    
    return file_paths

def is_file_empty(file_path: str) -> Tuple[bool, str]:
    """
    Check if a FASTA file is empty or contains no sequences.
    
    Returns:
        Tuple[bool, str]: (is_empty, reason)
    """
    try:
        # Check if file is zero bytes
        if os.path.getsize(file_path) == 0:
            return True, "zero_bytes"
        
        # Check if file contains any sequences
        with open(file_path, 'r') as f:
            # Read first few lines to check for sequence content
            content = f.read(1000)  # Read first 1KB
            if not content.strip():
                return True, "no_content"
        
        # Try to parse sequences
        sequences = list(SeqIO.parse(file_path, "fasta"))
        if len(sequences) == 0:
            return True, "no_sequences"
        
        # Check if all sequences are empty
        non_empty_sequences = [seq for seq in sequences if len(str(seq.seq).replace('-', '')) > 0]
        if len(non_empty_sequences) == 0:
            return True, "all_sequences_empty"
        
        return False, "not_empty"
        
    except Exception as e:
        # If we can't read the file, consider it problematic
        return True, f"read_error: {str(e)}"

# =============================================================================
# AT Content Analysis Functions
# =============================================================================
def calculate_at_content(sequence: str) -> float:
    """Calculate AT content for a sequence."""
    if not sequence:
        return 0.0
    sequence_clean = sequence.replace('-', '').upper()
    if not sequence_clean:
        return 0.0
    at_count = sequence_clean.count('A') + sequence_clean.count('T')
    return at_count / len(sequence_clean)

def vectorised_at_content(sequences_array: np.ndarray) -> np.ndarray:
    """Calculate AT content using vectorised operations."""
    if sequences_array.size == 0:
        return np.array([])
    
    non_gap_mask = sequences_array != '-'
    at_mask = (sequences_array == 'A') | (sequences_array == 'T')
    
    at_counts = np.sum(at_mask & non_gap_mask, axis=1)
    total_bases = np.sum(non_gap_mask, axis=1)
    
    at_content = np.divide(at_counts, total_bases, 
                          out=np.zeros_like(at_counts, dtype=float), 
                          where=total_bases > 0)
    
    # Clean up intermediate arrays
    del non_gap_mask, at_mask, at_counts, total_bases
    gc.collect()
    
    return at_content

def vectorised_compare_at_content(consensus_array: np.ndarray, sequences_array: np.ndarray) -> Tuple[np.ndarray, float, np.ndarray]:
    """Compare AT content using vectorised operations."""
    if sequences_array.size == 0 or consensus_array.size == 0:
        return np.array([]), 0.0, np.array([])
    
    cons_valid = consensus_array != '-'
    seq_valid = sequences_array != '-'
    overlap_mask = cons_valid & seq_valid
    
    # Calculate consensus AT content
    cons_at_mask = (consensus_array == 'A') | (consensus_array == 'T')
    cons_at_bases = np.sum(cons_at_mask & cons_valid)
    cons_total_bases = np.sum(cons_valid)
    consensus_at = cons_at_bases / cons_total_bases if cons_total_bases > 0 else 0.0
    
    # Calculate sequence AT content
    seq_at_mask = (sequences_array == 'A') | (sequences_array == 'T')
    seq_at_counts = np.sum(seq_at_mask & overlap_mask, axis=1)
    seq_total_counts = np.sum(overlap_mask, axis=1)
    
    query_at = np.divide(seq_at_counts, seq_total_counts,
                        out=np.zeros_like(seq_at_counts, dtype=float),
                        where=seq_total_counts > 0)
    
    at_differences = np.abs(query_at - consensus_at)
    
    # Clean up intermediate arrays
    del cons_valid, seq_valid, overlap_mask, cons_at_mask, seq_at_mask
    del seq_at_counts, seq_total_counts
    gc.collect()
    
    return query_at, consensus_at, at_differences

def calculate_position_frequencies(sequences: List[str]) -> List[Dict[str, float]]:
    """Calculate residue frequencies at each position using vectorised operations."""
    if not sequences:
        raise ValueError("No sequences provided")
    
    sequences_array = sequences_to_numpy_array(sequences)
    align_length = sequences_array.shape[1]
    position_freqs = []
    
    for i in range(align_length):
        column = sequences_array[:, i]
        residues = column[column != '-']
        
        if len(residues) > 0:
            unique, counts = np.unique(residues, return_counts=True)
            total = np.sum(counts)
            frequencies = {res: count/total for res, count in zip(unique, counts)}
        else:
            frequencies = {}
            
        position_freqs.append(frequencies)
    
    return position_freqs

def generate_consensus_sequence(alignment: MultipleSeqAlignment, 
                              threshold: float = 0.5) -> Tuple[str, List[Dict[str, float]]]:
    """Generate consensus sequence using optimised operations."""
    if not alignment:
        raise ValueError("Empty alignment provided")
    
    sequences = [str(record.seq) for record in alignment]
    frequencies = calculate_position_frequencies(sequences)
    
    consensus = []
    for pos_freqs in frequencies:
        if pos_freqs:
            most_common = max(pos_freqs.items(), key=lambda x: x[1])
            if most_common[1] >= threshold:
                consensus.append(most_common[0])
            else:
                consensus.append('-')
        else:
            consensus.append('-')
    
    return ''.join(consensus), frequencies

# =============================================================================
# AT Content Filtering Functions
# =============================================================================
def filter_at_content(records: List[SeqRecord],
                     consensus_seq: str,
                     threshold: float = 0.1,
                     mode: str = 'absolute') -> Tuple[List[SeqRecord], List[SeqRecord], Dict[str, float]]:
    """Filter sequences based on AT content with automatic chunking."""
    if not records or not consensus_seq:
        return records, [], {}
    
    # Calculate reasonable chunk size for memory management
    avg_length = len(consensus_seq)
    target_memory_mb = 500
    chunk_size = max(10, min(200, int(target_memory_mb / (avg_length * 8 * 3 / (1024 * 1024)))))
    
    # For very small files, use smaller chunks to avoid unnecessary overhead
    if len(records) < 50:
        chunk_size = len(records)
    
    # Process in chunks
    kept = []
    removed = []
    all_query_at = []
    all_at_differences = []
    
    for chunk_start in range(0, len(records), chunk_size):
        chunk_end = min(chunk_start + chunk_size, len(records))
        chunk_records = records[chunk_start:chunk_end]
        
        # Process this chunk with vectorised operations
        sequences = [str(record.seq).upper() for record in chunk_records]
        max_len = len(consensus_seq)
        padded_sequences = [seq.ljust(max_len, '-')[:max_len] for seq in sequences]
        sequences_array = sequences_to_numpy_array(padded_sequences)
        consensus_array = np.array(list(consensus_seq.upper()))
        
        # Calculate AT content differences
        query_at, cons_at, at_differences = vectorised_compare_at_content(
            consensus_array, sequences_array
        )
        
        # Store for statistics
        all_query_at.extend(query_at.tolist())
        all_at_differences.extend(at_differences.tolist())
        
        signed_differences = query_at - cons_at
        
        # Determine removal criteria based on mode
        if mode == 'absolute':
            remove_mask = at_differences > threshold
        elif mode == 'higher':
            remove_mask = signed_differences > threshold
        elif mode == 'lower':
            remove_mask = signed_differences < -threshold
        else:
            remove_mask = np.zeros(len(chunk_records), dtype=bool)
        
        # Filter this chunk
        for i, record in enumerate(chunk_records):
            if i < len(remove_mask) and remove_mask[i]:
                removed.append(record)
            else:
                kept.append(record)
        
        # Delete arrays
        del sequences_array, consensus_array, query_at, at_differences
        
        # Cleanup every few chunks
        if chunk_start % (chunk_size * 5) == 0:
            memory_cleanup()
    
    # Calculate summary statistics
    stats = {}
    if all_query_at:
        stats['consensus_at_content'] = cons_at
        stats['mean_query_at_content'] = np.mean(all_query_at)
        stats['std_query_at_content'] = np.std(all_query_at)
        stats['mean_at_difference'] = np.mean(all_at_differences)
        stats['max_at_difference'] = np.max(all_at_differences)
    
    # Final cleanup
    memory_cleanup()
    
    return kept, removed, stats

# =============================================================================
# Sequence Metrics Functions
# =============================================================================
def calculate_sequence_metrics(record: SeqRecord, 
                             consensus_seq: str,
                             consensus_at: float) -> Dict[str, float]:
    """Calculate AT content metrics for a sequence."""
    sequence = str(record.seq).upper()
    sequence_no_gaps = sequence.replace('-', '')
    
    # Basic metrics
    metrics = {
        'length_no_gaps': len(sequence_no_gaps),
        'length_with_gaps': len(sequence),
        'at_content': calculate_at_content(sequence_no_gaps),
        'consensus_at_content': consensus_at
    }
    
    # AT difference from consensus
    metrics['at_difference'] = abs(metrics['at_content'] - consensus_at)
    metrics['at_difference_signed'] = metrics['at_content'] - consensus_at
    
    # Coverage metrics
    if len(consensus_seq) > 0:
        non_gap_count = sum(1 for char in sequence if char != '-')
        coverage_percent = (non_gap_count / len(consensus_seq)) * 100
        metrics['coverage_percent'] = coverage_percent
    else:
        metrics['coverage_percent'] = 0.0
    
    return metrics

# =============================================================================
# Core Processing Functions
# =============================================================================
class ATFilterResult(NamedTuple):
    """Container for AT content filtering results"""
    kept_records: List[SeqRecord]
    removed_records: List[SeqRecord]
    consensus_at: float
    filter_stats: Dict[str, float]
    sequence_metrics: Dict[str, Dict[str, float]]

def process_at_filtering(alignment: MultipleSeqAlignment,
                        consensus_threshold: float = 0.5,
                        at_threshold: float = 0.1,
                        at_mode: str = 'absolute') -> ATFilterResult:
    """Apply AT content filtering to an alignment."""
    
    # Generate consensus sequence internally for AT content comparison
    consensus_seq, frequencies = generate_consensus_sequence(alignment, consensus_threshold)
    consensus_at = calculate_at_content(consensus_seq)
    
    # Apply AT content filtering
    kept_records, removed_records, filter_stats = filter_at_content(
        list(alignment), consensus_seq, threshold=at_threshold, mode=at_mode
    )
    
    # Calculate metrics for all sequences (using internal consensus)
    sequence_metrics = {}
    for record in alignment:
        sequence_metrics[record.id] = calculate_sequence_metrics(
            record, consensus_seq, consensus_at
        )
    
    return ATFilterResult(
        kept_records=kept_records,
        removed_records=removed_records,
        consensus_at=consensus_at,
        filter_stats=filter_stats,
        sequence_metrics=sequence_metrics
    )

# =============================================================================
# File I/O Functions
# =============================================================================
def write_sequences(records: List[SeqRecord], output_path: str) -> None:
    """Write sequences to FASTA file."""
    try:
        # Build entire content in memory first, then write once
        content = []
        for record in records:
            content.append(f">{record.id}\n{str(record.seq)}\n")
        
        with open(output_path, 'w') as handle:
            handle.write(''.join(content))
    except Exception as e:
        raise Exception(f"Error writing file {output_path}: {str(e)}")

def write_metrics_report(filter_result: ATFilterResult, output_path: str) -> None:
    """Write comprehensive AT content metrics report."""
    with open(output_path, 'w', newline='') as csvfile:
        # Define output fields
        fields = [
            'sequence_id',
            'length_no_gaps',
            'length_with_gaps',
            'at_content',
            'consensus_at_content',
            'at_difference',
            'at_difference_signed',
            'coverage_percent',
            'filter_status'
        ]
        
        writer = csv.writer(csvfile)
        writer.writerow(fields)
        
        # Create removal mapping
        removed_ids = {record.id for record in filter_result.removed_records}
        
        # Write all sequence metrics
        for seq_id, metrics in filter_result.sequence_metrics.items():
            row = [seq_id]
            
            # Add metrics
            for field in fields[1:-1]:  # Skip sequence_id and filter_status
                value = metrics.get(field, '')
                row.append(f"{value:.4f}" if isinstance(value, float) else str(value))
            
            # Add filter status
            status = 'removed' if seq_id in removed_ids else 'kept'
            row.append(status)
            
            writer.writerow(row)

def write_consensus_sequence(consensus_seq: str, base_name: str, output_path: str) -> None:
    """Write consensus sequence to FASTA file."""
    consensus_record = SeqRecord(
        Seq(consensus_seq),
        id=f"{base_name}_consensus",
        description="Generated by AT Content Filter"
    )
    write_sequences([consensus_record], output_path)

# =============================================================================
# Main Processing Functions
# =============================================================================
def process_fasta_file(input_file: str, 
                      output_dir: str,
                      consensus_threshold: float = 0.5,
                      at_threshold: float = 0.1,
                      at_mode: str = 'absolute') -> Dict[str, Any]:
    """Process a single FASTA file for AT content filtering."""
    
    # Get base name for clear identification
    base_name = get_base_filename(input_file)
    original_filename = os.path.basename(input_file)
    
    # Get file size for logging
    try:
        file_size = os.path.getsize(input_file)
        file_size_str = format_file_size(file_size)
    except OSError:
        file_size_str = "unknown"
    
    # Check if file is empty before any processing
    is_empty, empty_reason = is_file_empty(input_file)
    if is_empty:
        return {
            'status': 'skipped',
            'reason': empty_reason,
            'filename': original_filename,
            'base_name': base_name,
            'file_size': file_size_str,
            'sequence_count': 0
        }
    
    # Create output subdirectories
    output_subdirs = create_output_subdirectories(output_dir)
    
    # Set up output paths
    output_paths = {
        'filtered': os.path.join(output_subdirs['filtered_seqs'], f"{base_name}_at_filtered.fasta"),
        'metrics': os.path.join(output_subdirs['metrics'], f"{base_name}_at_metrics.csv"),
        'log': os.path.join(output_subdirs['logs'], f"{base_name}_log.txt")
    }
    
    # Initialise statistics
    stats = {
        'status': 'success',
        'input_sequences': 0,
        'kept_sequences': 0,
        'removed_sequences': 0,
        'processing_time': None,
        'filename': original_filename,
        'base_name': base_name,
        'file_size': file_size_str,
        'sequence_count': 0,
        'consensus_at_content': 0.0,
        'mean_query_at_content': 0.0,
        'filter_mode': at_mode,
        'at_threshold': at_threshold
    }
    
    start_time = datetime.now()
    
    with open(output_paths['log'], 'w') as log_file:
        try:
            # Write log header
            separator = "=" * 80
            log_message(separator, log_file)
            log_message(f"AT Content Filter - Processing Sample: {base_name}", log_file)
            log_message(separator, log_file)
            log_message(f"Started at: {start_time}", log_file)
            log_message(f"Input file: {input_file}", log_file)
            log_message(f"File size: {file_size_str}\n", log_file)
            
            # Configuration section
            log_message("Configuration:", log_file)
            log_message(f"- Consensus threshold: {consensus_threshold}", log_file)
            log_message(f"- AT threshold: {at_threshold}", log_file)
            log_message(f"- AT mode: {at_mode}\n", log_file)
            
            # Read input alignment
            alignment = AlignIO.read(input_file, "fasta")
            stats['input_sequences'] = len(alignment)
            stats['sequence_count'] = len(alignment)
            log_message(f"Input sequences: {stats['input_sequences']}", log_file)
            
            # Process AT filtering
            log_message(separator, log_file)
            log_message("AT Content Analysis:", log_file)
            log_message(separator, log_file)
            
            filter_result = process_at_filtering(
                alignment,
                consensus_threshold=consensus_threshold,
                at_threshold=at_threshold,
                at_mode=at_mode
            )
            
            # Update statistics
            stats['kept_sequences'] = len(filter_result.kept_records)
            stats['removed_sequences'] = len(filter_result.removed_records)
            stats['consensus_at_content'] = filter_result.consensus_at
            
            if filter_result.filter_stats:
                stats['mean_query_at_content'] = filter_result.filter_stats.get('mean_query_at_content', 0.0)
            
            # Log results
            log_message(f"Consensus AT content: {filter_result.consensus_at:.4f}", log_file)
            if filter_result.filter_stats:
                log_message(f"Mean query AT content: {filter_result.filter_stats.get('mean_query_at_content', 0.0):.4f}", log_file)
                log_message(f"Mean AT difference: {filter_result.filter_stats.get('mean_at_difference', 0.0):.4f}", log_file)
                log_message(f"Max AT difference: {filter_result.filter_stats.get('max_at_difference', 0.0):.4f}", log_file)
            
            log_message(f"\nFiltering mode: {at_mode}", log_file)
            log_message(f"AT threshold: {at_threshold}", log_file)
            log_message(f"Sequences kept: {stats['kept_sequences']}", log_file)
            log_message(f"Sequences removed: {stats['removed_sequences']}", log_file)
            
            # Write output files
            log_message(separator, log_file)
            log_message("Writing Output Files:", log_file)
            log_message(separator, log_file)
            
            if filter_result.kept_records:
                write_sequences(filter_result.kept_records, output_paths['filtered'])
                log_message(f"Filtered sequences written to: {output_paths['filtered']}", log_file)
            
            write_metrics_report(filter_result, output_paths['metrics'])
            log_message(f"Metrics report written to: {output_paths['metrics']}", log_file)
            
            # Final timing
            stats['processing_time'] = (datetime.now() - start_time).total_seconds()
            log_message(f"\nProcessing completed in {stats['processing_time']:.2f} seconds", log_file)
            log_message(separator + "\n", log_file)
            
            # Cleanup
            del alignment, filter_result
            memory_cleanup()
            
        except Exception as e:
            stats['status'] = 'error'
            error_msg = f"Error processing {original_filename}: {str(e)}"
            log_message(error_msg, log_file, stdout=False)
            log_message(f"Traceback:\n{traceback.format_exc()}", log_file)
            stats['error_message'] = error_msg
    
    return stats

def process_large_files_sequentially(large_files: List[FileInfo], 
                                    args: argparse.Namespace) -> Tuple[List[str], List[str]]:
    """Process large files one at a time to prevent memory exhaustion."""
    processed_files = []
    error_files = []
    
    if not large_files:
        return processed_files, error_files
    
    print(f"\nProcessing {len(large_files)} large files sequentially (1 at a time)")
    print("Sequential processing - using single worker per large file for memory management")
    
    # Show file details for large files
    total_size = sum(f.size_bytes for f in large_files)
    print(f"Total size of large files: {format_file_size(total_size)}")
    
    # Create partial function with common arguments
    process_func = partial(
        process_fasta_file,
        output_dir=args.output_dir,
        consensus_threshold=args.consensus_threshold,
        at_threshold=args.at_threshold,
        at_mode=args.at_mode
    )
    
    # Process each large file individually
    for i, file_info in enumerate(large_files, 1):
        filename = file_info.basename
        
        print(f"\n--- Processing large file {i}/{len(large_files)} ---")
        
        # Log start with file size info
        timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        start_msg = f"[{timestamp}] STARTING {filename} | Size: {file_info.size_str}"
        print(start_msg, flush=True)
        
        try:
            # Process single file
            result = process_func(file_info.path)
            
            timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
            status = result.get('status', 'unknown')
            
            if status == 'error':
                error_files.append(filename)
                error_msg = result.get('error_message', 'Unknown error')
                log_msg = f"[{timestamp}] ERROR processing {filename}: {error_msg}"
                print(f"\nERROR: {log_msg}", flush=True)
                
            elif status == 'skipped':
                reason = result.get('reason', 'unknown')
                log_msg = f"[{timestamp}] SKIPPED {filename} ({reason})"
                print(log_msg, flush=True)
                
            else:  # success
                processed_files.append(filename)
                kept = result.get('kept_sequences', 0)
                removed = result.get('removed_sequences', 0)
                at_content = result.get('consensus_at_content', 0.0)
                processing_time = result.get('processing_time', 0)
                log_msg = f"[{timestamp}] COMPLETED {filename} | Size: {file_info.size_str} | Kept: {kept} | Removed: {removed} | Consensus AT: {at_content:.3f} | Time: {processing_time:.1f}s"
                print(log_msg, flush=True)
                
        except Exception as e:
            error_files.append(filename)
            timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
            log_msg = f"[{timestamp}] ERROR processing {filename}: {str(e)}"
            print(f"\nERROR: {log_msg}", flush=True)
        
        # Memory cleanup after each large file
        if i < len(large_files):  # Don't cleanup after the last file
            print("Performing memory cleanup...")
            memory_cleanup()
            time.sleep(0.5)  # Brief pause to allow memory cleanup
    
    print(f"\nLarge files batch completed: {len(processed_files)} processed, {len(error_files)} errors")
    return processed_files, error_files
def process_normal_files_parallel(normal_files: List[FileInfo], 
                                  num_threads: int,
                                  args: argparse.Namespace) -> Tuple[List[str], List[str]]:
    """Process normal files with full parallelism."""
    processed_files = []
    error_files = []
    
    batch_size = len(normal_files)
    if batch_size == 0:
        return processed_files, error_files
    
    print(f"\nProcessing {batch_size} normal files with {num_threads} workers")
    
    # Create partial function with common arguments
    process_func = partial(
        process_fasta_file,
        output_dir=args.output_dir,
        consensus_threshold=args.consensus_threshold,
        at_threshold=args.at_threshold,
        at_mode=args.at_mode
    )
    
    # Process files with ProcessPoolExecutor
    with concurrent.futures.ProcessPoolExecutor(max_workers=num_threads) as executor:
        futures = {}
        
        # Submit jobs
        for file_info in normal_files:
            filename = file_info.basename
            
            # Log start with file size info
            timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
            start_msg = f"[{timestamp}] STARTING {filename} | Size: {file_info.size_str}"
            print(start_msg, flush=True)
            
            future = executor.submit(process_func, file_info.path)
            futures[future] = (filename, file_info.size_str)
        
        # Process results as they complete
        try:
            for future in concurrent.futures.as_completed(futures):
                filename, file_size_str = futures[future]
                timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
                
                try:
                    result = future.result()
                    status = result.get('status', 'unknown')
                    
                    if status == 'error':
                        error_files.append(filename)
                        error_msg = result.get('error_message', 'Unknown error')
                        log_msg = f"[{timestamp}] ERROR processing {filename}: {error_msg}"
                        print(f"\nERROR: {log_msg}", flush=True)
                        
                    elif status == 'skipped':
                        reason = result.get('reason', 'unknown')
                        log_msg = f"[{timestamp}] SKIPPED {filename} ({reason})"
                        print(log_msg, flush=True)
                        
                    else:  # success
                        processed_files.append(filename)
                        kept = result.get('kept_sequences', 0)
                        removed = result.get('removed_sequences', 0)
                        at_content = result.get('consensus_at_content', 0.0)
                        log_msg = f"[{timestamp}] COMPLETED {filename} | Size: {file_size_str} | Kept: {kept} | Removed: {removed} | Consensus AT: {at_content:.3f}"
                        print(log_msg, flush=True)
                    
                except Exception as e:
                    error_files.append(filename)
                    log_msg = f"[{timestamp}] ERROR processing {filename}: {str(e)}"
                    print(f"\nERROR: {log_msg}", flush=True)
                    
        except KeyboardInterrupt:
            print(f"\nProcessing interrupted by user. Shutting down...")
            executor.shutdown(wait=True)
            raise
    
    print(f"Normal files batch completed: {len(processed_files)} processed, {len(error_files)} errors")
    return processed_files, error_files

def parallel_process_files(args: argparse.Namespace) -> None:
    """Process multiple FASTA files in parallel with memory-aware batching."""
    
    # Create output subdirectories
    output_subdirs = create_output_subdirectories(args.output_dir)
    
    # Read input files from log file
    try:
        all_fasta_files = read_input_files_from_log(args.input_files)
    except Exception as e:
        print(f"Error reading input files: {str(e)}")
        sys.exit(1)
    
    # Classify files by size
    large_files, normal_files, empty_files = classify_files_by_size(all_fasta_files)
    
    print("=" * 80)
    
    if not large_files and not normal_files:
        print("No valid files to process!")
        return
    
    print(f"AT Content Filter Processing - {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 80)
    
    # Log configuration
    print("Configuration:")
    print(f"- Consensus threshold: {args.consensus_threshold}")
    print(f"- AT threshold: {args.at_threshold}")
    print(f"- AT mode: {args.at_mode}")
    print(f"- Standard threads: {args.threads}")
    print(f"- Large file threshold: {format_file_size(LARGE_FILE_THRESHOLD)}")
    print(f"- Large file processing: Sequential (1 at a time)")
    print("=" * 80)
    
    # Log empty files if any
    if empty_files:
        print(f"Pre-filtered empty files ({len(empty_files)}):")
        for filename in empty_files[:10]:  # Show first 10
            print(f"- {filename}")
        if len(empty_files) > 10:
            print(f"... and {len(empty_files) - 10} more")
        print("=" * 80)
    
    # Process files in batches
    start_time = datetime.now()
    all_processed_files = []
    all_error_files = []
    
    try:
        # Process large files sequentially if any exist
        if large_files:
            print(f"\nPHASE 1: Processing {len(large_files)} large files sequentially (≥{format_file_size(LARGE_FILE_THRESHOLD)})")
            print("=" * 80)
            
            large_processed, large_errors = process_large_files_sequentially(
                large_files, 
                args
            )
            all_processed_files.extend(large_processed)
            all_error_files.extend(large_errors)
            
            # Memory cleanup between phases
            print("\nPerforming memory cleanup between processing phases...")
            memory_cleanup()
            time.sleep(1)  # Brief pause to allow memory cleanup
        
        # Process normal files with full parallelism
        if normal_files:
            print(f"\nPHASE 2: Processing {len(normal_files)} normal files in parallel (<{format_file_size(LARGE_FILE_THRESHOLD)})")
            print("=" * 80)
            
            normal_processed, normal_errors = process_normal_files_parallel(
                normal_files, 
                args.threads, 
                args
            )
            all_processed_files.extend(normal_processed)
            all_error_files.extend(normal_errors)
    
    except KeyboardInterrupt:
        print("Processing interrupted by user")
        all_processed_files, all_error_files = [], []
    
    # Final summary
    elapsed_time = (datetime.now() - start_time).total_seconds()
    hours = int(elapsed_time // 3600)
    minutes = int((elapsed_time % 3600) // 60)
    seconds = int(elapsed_time % 60)
    
    print("=" * 80)
    print(f"Processing complete in {hours:02d}:{minutes:02d}:{seconds:02d}")
    print(f"Total input files: {len(all_fasta_files)}")
    print(f"Empty files skipped: {len(empty_files)}")
    print(f"Large files processed: {len([f for f in all_processed_files if f in [lf.basename for lf in large_files]])}/{len(large_files)}")
    print(f"Normal files processed: {len([f for f in all_processed_files if f in [nf.basename for nf in normal_files]])}/{len(normal_files)}")
    print(f"Total files processed successfully: {len(all_processed_files)}/{len(large_files) + len(normal_files)}")
    
    if all_error_files:
        print(f"Errors occurred in {len(all_error_files)} files:")
        for filename in all_error_files[:20]:  # Show first 20 errors
            print(f"- {filename}")
        if len(all_error_files) > 20:
            print(f"... and {len(all_error_files) - 20} more errors")
    
    print("=" * 80)
    
    # Generate output files
    print("\nGenerating output files...")
    
    # Write filtered paths log
    if all_processed_files:
        filtered_paths = write_filtered_paths_log(args.output_dir, all_processed_files, args.filtered_files_list)
    
    # Generate summary statistics
    print("Generating summary statistics...")
    generate_summary_statistics(args.output_dir)
    
    # Compress metrics directory
    print("\nCompressing metrics files...")
    compression_success = compress_metrics_directory(args.output_dir)
    
    if not compression_success:
        print("Warning: Metrics directory compression failed - individual files remain available")

def write_filtered_paths_log(output_dir: str, processed_files: List[str], filename: str = 'at_filtered_paths.log') -> None:
    """Write filtered paths log to exact location Snakemake expects."""
    filtered_seqs_dir = os.path.join(output_dir, 'at_filtered_sequences')
    
    # Handle absolute vs relative paths correctly
    if os.path.isabs(filename):
        paths_log = filename
    else:
        paths_log = os.path.join(output_dir, filename)
    
    successful_paths = []
    
    # Check which files were successfully processed and have output
    for file_basename in processed_files:
        # Extract base name properly
        base_name = get_base_filename(file_basename)
        filtered_file = os.path.join(filtered_seqs_dir, f"{base_name}_at_filtered.fasta")
        
        if os.path.exists(filtered_file) and os.path.getsize(filtered_file) > 0:
            successful_paths.append(filtered_file)
    
    # Write paths log file
    try:
        with open(paths_log, 'w') as f:
            for path in successful_paths:
                f.write(path + '\n')
        
        print(f"Filtered alignment paths written to: {paths_log}")
        print(f"Contains {len(successful_paths)} successfully filtered alignment files")
        
    except Exception as e:
        print(f"Error writing filtered paths log: {str(e)}")
    
    return successful_paths

def compress_metrics_directory(output_dir: str) -> bool:
    """
    Compress the metrics directory into a tar.gz file and remove the original directory.
    
    Returns:
        bool: True if compression was successful, False otherwise
    """
    metrics_dir = os.path.join(output_dir, 'metrics')
    tar_path = os.path.join(output_dir, 'concatenated_individual_at_metrics.tar.gz')
    
    if not os.path.exists(metrics_dir):
        print(f"Metrics directory not found: {metrics_dir}")
        return False
    
    # Check if directory has any files
    try:
        files_in_dir = os.listdir(metrics_dir)
        csv_files = [f for f in files_in_dir if f.endswith('.csv')]
        
        if not csv_files:
            print("No CSV files found in metrics directory to compress")
            return False
        
        print(f"Compressing metrics directory...")
        print(f"Found {len(csv_files)} metrics files to compress")
        
        # Create tar.gz file
        with tarfile.open(tar_path, 'w:gz') as tar:
            # Add the entire metrics directory, but use 'metrics' as the arcname
            tar.add(metrics_dir, arcname='metrics')
        
        # Check if tar file was created successfully
        if os.path.exists(tar_path) and os.path.getsize(tar_path) > 0:
            # Get file sizes for reporting
            original_size = sum(
                os.path.getsize(os.path.join(metrics_dir, f)) 
                for f in os.listdir(metrics_dir) 
                if os.path.isfile(os.path.join(metrics_dir, f))
            )
            compressed_size = os.path.getsize(tar_path)
            compression_ratio = (1 - compressed_size / original_size) * 100 if original_size > 0 else 0
            
            print(f"Compression successful!")
            print(f"Original size: {format_file_size(original_size)}")
            print(f"Compressed size: {format_file_size(compressed_size)}")
            print(f"Compression ratio: {compression_ratio:.1f}%")
            print(f"Compressed file: {tar_path}")
            
            # Remove original directory
            shutil.rmtree(metrics_dir)
            print(f"Original metrics directory removed")
            
            return True
        else:
            print("Error: Compressed file was not created successfully")
            return False
            
    except Exception as e:
        print(f"Error compressing metrics directory: {str(e)}")
        return False

def generate_summary_statistics(output_dir: str) -> None:
    """Generate summary statistics from processing logs."""
    logs_dir = os.path.join(output_dir, 'logs')
    output_file = os.path.join(output_dir, 'at_filter_summary.csv')
    
    if not os.path.exists(logs_dir):
        print(f"Logs directory not found: {logs_dir}")
        return
    
    # Get list of log files (exclude main log)
    log_files = [f for f in os.listdir(logs_dir) 
                if f.endswith('_log.txt') and f != 'at_filter_main.log']
    
    if not log_files:
        print("No individual log files found")
        return
    
    # Extract statistics from each log file
    all_stats = []
    
    for log_file in log_files:
        log_path = os.path.join(logs_dir, log_file)
        stats = extract_log_statistics(log_path)
        if stats:
            all_stats.append(stats)
    
    # Write summary CSV
    if all_stats:
        fieldnames = [
            'sample_name', 'input_sequences', 'kept_sequences', 'removed_sequences',
            'consensus_at_content', 'mean_query_at_content', 'filter_mode', 
            'at_threshold', 'processing_time_seconds'
        ]
        
        with open(output_file, 'w', newline='') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            
            for stats in all_stats:
                # Ensure all fields are present
                row = {field: stats.get(field, '') for field in fieldnames}
                writer.writerow(row)
        
        print(f"Summary statistics written to: {output_file}")
        
        # Print basic summary to console
        total_input = sum(stats.get('input_sequences', 0) for stats in all_stats)
        total_kept = sum(stats.get('kept_sequences', 0) for stats in all_stats)
        total_removed = sum(stats.get('removed_sequences', 0) for stats in all_stats)
        
        print(f"\nOverall Summary:")
        print(f"- Samples processed: {len(all_stats)}")
        print(f"- Total input sequences: {total_input}")
        print(f"- Total kept sequences: {total_kept}")
        print(f"- Total removed sequences: {total_removed}")
        if total_input > 0:
            print(f"- Retention rate: {(total_kept/total_input)*100:.1f}%")
    else:
        print("No statistics to write")

def extract_log_statistics(log_file_path: str) -> Dict[str, Any]:
    """Extract statistics from a log file."""
    import re
    
    stats = {}
    
    try:
        with open(log_file_path, 'r') as f:
            log_content = f.read()
        
        # Extract sample name
        sample_match = re.search(r"AT Content Filter - Processing Sample: (.+)$", log_content, re.MULTILINE)
        if sample_match:
            stats['sample_name'] = sample_match.group(1)
        
        # Extract numerical statistics
        patterns = {
            'input_sequences': r"Input sequences: (\d+)",
            'kept_sequences': r"Sequences kept: (\d+)",
            'removed_sequences': r"Sequences removed: (\d+)",
            'consensus_at_content': r"Consensus AT content: ([\d.]+)",
            'mean_query_at_content': r"Mean query AT content: ([\d.]+)",
            'processing_time_seconds': r"Processing completed in ([\d.]+) seconds"
        }
        
        for key, pattern in patterns.items():
            match = re.search(pattern, log_content)
            if match:
                value = match.group(1)
                if key in ['input_sequences', 'kept_sequences', 'removed_sequences']:
                    stats[key] = int(value)
                else:
                    stats[key] = float(value)
        
        # Extract configuration
        mode_match = re.search(r"Filtering mode: (\w+)", log_content)
        if mode_match:
            stats['filter_mode'] = mode_match.group(1)
        
        threshold_match = re.search(r"AT threshold: ([\d.]+)", log_content)
        if threshold_match:
            stats['at_threshold'] = float(threshold_match.group(1))
    
    except Exception as e:
        print(f"Error extracting statistics from {log_file_path}: {str(e)}")
        return {}
    
    return stats

# =============================================================================
# Command Line Interface
# =============================================================================
def parse_arguments() -> argparse.Namespace:
    """Set up command line argument parsing."""
    parser = argparse.ArgumentParser(
        description='AT Content FASTA Sequence Filter - Specialized tool for filtering '
                   'sequences based on AT content analysis compared to consensus sequence. '
                   'Features memory-aware processing for large files.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # Get default number of threads
    default_threads = get_default_threads()
    
    # Required arguments
    parser.add_argument(
        '-i', '--input_files',
        required=True,
        help='Log file containing paths to input FASTA alignment files (one per line)'
    )
    
    parser.add_argument(
        '-o', '--output_dir',
        required=True,
        help='Directory where output files will be saved'
    )
    
    # File naming options
    parser.add_argument(
        '--filtered-files-list',
        type=str,
        default='at_filtered_paths.log',
        help='Name for the output file containing paths to filtered alignment files'
    )
    
    # Processing parameters
    parser.add_argument(
        '-t', '--threads',
        type=int,
        default=default_threads,
        help='Number of threads to use for parallel processing. '
             'Defaults to SLURM_CPUS_PER_TASK if available, '
             'otherwise uses all available CPU cores. '
             f'Large files (≥{format_file_size(LARGE_FILE_THRESHOLD)}) are processed '
             f'sequentially (1 at a time) to prevent memory exhaustion.'
    )
    
    # AT content filtering parameters
    parser.add_argument(
        '-d', '--at_threshold',
        type=float,
        default=0.1,
        help='Maximum allowed AT content difference from consensus (0.0-1.0)'
    )
    
    parser.add_argument(
        '-m','--at_mode',
        type=str,
        choices=['absolute', 'higher', 'lower'],
        default='absolute',
        help='AT content filtering mode: "absolute" removes sequences if AT content differs '
             'from consensus by more than threshold in either direction, "higher" removes only '
             'sequences with AT content above consensus + threshold, "lower" removes only '
             'sequences with AT content below consensus - threshold'
    )
    
    # Consensus generation parameters
    parser.add_argument(
        '-c', '--consensus_threshold',
        type=float,
        default=0.5,
        help='Threshold for consensus sequence generation (0.0-1.0)'
    )
    
    args = parser.parse_args()
    
    # Validate input files log
    if not os.path.isfile(args.input_files):
        parser.error(f"Input files log does not exist: {args.input_files}")
    
    # Validate thresholds
    if not 0 <= args.at_threshold <= 1:
        parser.error("AT threshold must be between 0.0 and 1.0")
    
    if not 0 < args.consensus_threshold <= 1:
        parser.error("Consensus threshold must be between 0.0 and 1.0")
    
    if args.threads < 1:
        parser.error("Number of threads must be at least 1")
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    return args

# =============================================================================
# Main Entry Point
# =============================================================================
if __name__ == "__main__":
    try:
        # Parse arguments
        args = parse_arguments()
        
        # Set the multiprocessing start method early
        try:
            mp.set_start_method('spawn', force=True)
        except RuntimeError:
            pass
        
        # Display configuration
        print("AT Content FASTA Sequence Filter")
        print("=" * 50)
        print(f"Input files log: {args.input_files}")
        print(f"Output directory: {args.output_dir}")
        print(f"Filtered files list: {args.filtered_files_list}")
        print(f"AT threshold: {args.at_threshold}")
        print(f"AT mode: {args.at_mode}")
        print(f"Consensus threshold: {args.consensus_threshold}")
        print(f"Standard threads: {args.threads}")
        print(f"Large file threshold: {format_file_size(LARGE_FILE_THRESHOLD)}")
        print(f"Large file processing: Sequential (1 at a time)")
        print("=" * 50)
        
        # Process files
        parallel_process_files(args)
        
        print("\nAT content filtering completed successfully!")
        print("\nOutput Files Generated:")
        print("=" * 50)
        print(f"• Filtered sequences: {args.output_dir}/at_filtered_sequences/")
        print(f"• Individual metrics: {args.output_dir}/concatenated_individual_at_metrics.tar.gz")
        print(f"• Processing logs: {args.output_dir}/logs/")
        print(f"• Filtered paths log: {args.output_dir}/{args.filtered_files_list}")
        print(f"• Summary statistics: {args.output_dir}/at_filter_summary.csv")
        print("=" * 50)
        
    except KeyboardInterrupt:
        print("\nProcessing interrupted by user.")
        sys.exit(1)
    except Exception as e:
        print(f"Error: {str(e)}")
        sys.exit(1)
