#!/usr/bin/env python3
"""
Statistical Outlier Filter
Removes sequences that are statistical outliers compared to consensus.

1. Creates a reference consensus sequence using the most common nucleotide at each position
2. Uses a threshold (default 0.5) to determine if a position has a clear consensus
3. The script calculates two types of deviation scores for each sequence:
    - Unweighted Deviation: Simple count of positions where the sequence differs from consensus. (number of differences) / (total valid positions). Treats all positions equally
    - Weighted Deviation: Considers how conserved each position is in the alignment. (conservation-weighted differences) / (total conservation weight). Differences at highly conserved positions are weighted more heavily
4. Calculates percentile thresholds for both deviation scores across all sequences. Sequences exceeding either threshold are flagged as outliers. Removes sequences that are statistical outliers in either scoring method.
"""


import os
import sys
import csv
import argparse
import numpy as np
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import List, Tuple, Dict
import gc
from datetime import datetime
from collections import Counter

def log_message(message: str, log_file=None, stdout=False):
    """Log message to file and optionally stdout"""
    timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    formatted_msg = f"[{timestamp}] {message}"
    
    if log_file:
        log_file.write(formatted_msg + '\n')
        log_file.flush()
    
    if stdout:
        print(formatted_msg, flush=True)

def calculate_position_frequencies(sequences: List[str]) -> List[Dict[str, float]]:
    """Calculate residue frequencies at each position"""
    if not sequences:
        return []
    
    max_len = max(len(seq) for seq in sequences)
    position_freqs = []
    
    for i in range(max_len):
        residues = []
        for seq in sequences:
            if i < len(seq) and seq[i] != '-':
                residues.append(seq[i])
        
        if residues:
            counter = Counter(residues)
            total = len(residues)
            frequencies = {res: count/total for res, count in counter.items()}
        else:
            frequencies = {}
            
        position_freqs.append(frequencies)
    
    return position_freqs

def generate_consensus_sequence(sequences: List[str], threshold: float = 0.5) -> Tuple[str, List[Dict[str, float]]]:
    """Generate consensus sequence and position frequencies"""
    if not sequences:
        return "", []
    
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

def calculate_unweighted_deviation(sequence: str, reference: str) -> float:
    """Calculate unweighted deviation score"""
    if len(sequence) != len(reference):
        min_len = min(len(sequence), len(reference))
        sequence = sequence[:min_len]
        reference = reference[:min_len]
    
    differences = 0
    valid_positions = 0
    
    for seq_res, ref_res in zip(sequence, reference):
        if seq_res != '-' and ref_res != '-':
            valid_positions += 1
            if seq_res != ref_res:
                differences += 1
    
    return differences / valid_positions if valid_positions > 0 else 0.0

def calculate_weighted_deviation(sequence: str, reference: str, frequencies: List[Dict[str, float]]) -> float:
    """Calculate weighted deviation score based on conservation"""
    if len(sequence) != len(reference):
        min_len = min(len(sequence), len(reference))
        sequence = sequence[:min_len]
        reference = reference[:min_len]
    
    if len(frequencies) < len(sequence):
        frequencies.extend([{}] * (len(sequence) - len(frequencies)))
    
    total_score = 0.0
    total_weight = 0.0
    
    for i, (seq_res, ref_res) in enumerate(zip(sequence, reference)):
        if seq_res != '-' and ref_res != '-' and i < len(frequencies):
            conservation_weight = frequencies[i].get(ref_res, 0)
            total_weight += conservation_weight
            
            if seq_res != ref_res:
                total_score += conservation_weight
    
    return total_score / total_weight if total_weight > 0 else 0.0

def process_single_file(file_path: str, outlier_percentile: float, consensus_threshold: float, output_dir: str) -> Dict:
    """Process a single FASTA file for statistical outlier filtering"""
    try:
        base_name = os.path.splitext(os.path.basename(file_path))[0]
        if '_human_filtered' in base_name:
            base_name = base_name.replace('_human_filtered', '')
        elif '_at_filtered' in base_name:
            base_name = base_name.replace('_at_filtered', '')
        
        # Check if file is empty
        if os.path.getsize(file_path) == 0:
            return {
                'file_path': file_path,
                'base_name': base_name,
                'status': 'skipped',
                'reason': 'empty_file',
                'input_count': 0,
                'kept_count': 0,
                'removed_count': 0,
                'removed_sequences': []
            }
        
        # Read sequences
        try:
            records = list(SeqIO.parse(file_path, "fasta"))
        except Exception as e:
            return {
                'file_path': file_path,
                'base_name': base_name,
                'status': 'error',
                'reason': f'parse_error: {str(e)}',
                'input_count': 0,
                'kept_count': 0,
                'removed_count': 0,
                'removed_sequences': []
            }
        
        if not records:
            return {
                'file_path': file_path,
                'base_name': base_name,
                'status': 'skipped',
                'reason': 'no_sequences',
                'input_count': 0,
                'kept_count': 0,
                'removed_count': 0,
                'removed_sequences': []
            }
        
        input_count = len(records)
        
        # Generate consensus sequence
        sequences = [str(record.seq).upper() for record in records]
        consensus_seq, frequencies = generate_consensus_sequence(sequences, consensus_threshold)
        
        if not consensus_seq:
            return {
                'file_path': file_path,
                'base_name': base_name,
                'status': 'error',
                'reason': 'consensus_generation_failed',
                'input_count': input_count,
                'kept_count': 0,
                'removed_count': 0,
                'removed_sequences': []
            }
        
        # Calculate deviation scores for all sequences
        deviation_scores = []
        for record in records:
            sequence = str(record.seq).upper()
            
            unweighted_dev = calculate_unweighted_deviation(sequence, consensus_seq)
            weighted_dev = calculate_weighted_deviation(sequence, consensus_seq, frequencies)
            
            deviation_scores.append({
                'record': record,
                'unweighted_deviation': unweighted_dev,
                'weighted_deviation': weighted_dev
            })
        
        # Calculate outlier thresholds
        unweighted_scores = [s['unweighted_deviation'] for s in deviation_scores if s['unweighted_deviation'] > 0]
        weighted_scores = [s['weighted_deviation'] for s in deviation_scores if s['weighted_deviation'] > 0]
        
        unweighted_threshold = np.percentile(unweighted_scores, outlier_percentile) if unweighted_scores else float('inf')
        weighted_threshold = np.percentile(weighted_scores, outlier_percentile) if weighted_scores else float('inf')
        
        # Filter sequences
        kept_records = []
        removed_sequences = []
        
        for score_info in deviation_scores:
            record = score_info['record']
            unweighted_dev = score_info['unweighted_deviation']
            weighted_dev = score_info['weighted_deviation']
            
            # Check if sequence is an outlier
            is_outlier = (unweighted_dev > unweighted_threshold or weighted_dev > weighted_threshold)
            
            if is_outlier:
                removal_reason = 'statistical_outlier'
                if unweighted_dev > unweighted_threshold and weighted_dev > weighted_threshold:
                    removal_reason = 'statistical_outlier_both'
                elif unweighted_dev > unweighted_threshold:
                    removal_reason = 'statistical_outlier_unweighted'
                elif weighted_dev > weighted_threshold:
                    removal_reason = 'statistical_outlier_weighted'
                
                removed_sequences.append({
                    'sequence_id': record.id,
                    'removal_reason': removal_reason,
                    'unweighted_deviation': unweighted_dev,
                    'weighted_deviation': weighted_dev,
                    'unweighted_threshold': unweighted_threshold,
                    'weighted_threshold': weighted_threshold
                })
            else:
                kept_records.append(record)
        
        # Write filtered sequences
        if kept_records:
            output_file = os.path.join(output_dir, f"{base_name}_outlier_filtered.fasta")
            with open(output_file, 'w') as handle:
                SeqIO.write(kept_records, handle, "fasta")
        
        return {
            'file_path': file_path,
            'base_name': base_name,
            'status': 'success',
            'reason': 'processed',
            'input_count': input_count,
            'kept_count': len(kept_records),
            'removed_count': len(removed_sequences),
            'removed_sequences': removed_sequences,
            'output_file': os.path.join(output_dir, f"{base_name}_outlier_filtered.fasta") if kept_records else None,
            'thresholds': {
                'unweighted_threshold': unweighted_threshold,
                'weighted_threshold': weighted_threshold,
                'outlier_percentile': outlier_percentile
            }
        }
        
    except Exception as e:
        return {
            'file_path': file_path,
            'base_name': base_name,
            'status': 'error',
            'reason': f'processing_error: {str(e)}',
            'input_count': 0,
            'kept_count': 0,
            'removed_count': 0,
            'removed_sequences': []
        }

def main():
    parser = argparse.ArgumentParser(description='Filter sequences based on statistical outlier detection')
    parser.add_argument('--input-files-list', required=True, help='File containing list of FASTA files to process')
    parser.add_argument('--output-dir', required=True, help='Output directory for filtered files')
    parser.add_argument('--filtered-files-list', required=True, help='Output file listing successfully filtered files')
    parser.add_argument('--summary-csv', required=True, help='Output CSV file with file-level summary metrics')
    parser.add_argument('--metrics-csv', required=True, help='Output CSV file with detailed sequence removal metrics')
    parser.add_argument('--outlier-percentile', type=float, default=90.0, help='Percentile threshold for outlier detection')
    parser.add_argument('--consensus-threshold', type=float, default=0.5, help='Consensus generation threshold')
    parser.add_argument('--threads', type=int, default=1, help='Number of threads for parallel processing')
    
    args = parser.parse_args()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Read input files
    with open(args.input_files_list, 'r') as f:
        input_files = [line.strip() for line in f if line.strip()]
    
    print(f"Processing {len(input_files)} files with outlier percentile {args.outlier_percentile}")
    
    # Process files in parallel
    results = []
    with ProcessPoolExecutor(max_workers=args.threads) as executor:
        future_to_file = {
            executor.submit(process_single_file, file_path, args.outlier_percentile, args.consensus_threshold, args.output_dir): file_path
            for file_path in input_files
        }
        
        for future in as_completed(future_to_file):
            file_path = future_to_file[future]
            try:
                result = future.result()
                results.append(result)
                
                if result['status'] == 'success':
                    thresholds = result.get('thresholds', {})
                    print(f"✓ {result['base_name']}: {result['kept_count']}/{result['input_count']} sequences kept")
                    print(f"  Thresholds: unweighted={thresholds.get('unweighted_threshold', 'N/A'):.4f}, weighted={thresholds.get('weighted_threshold', 'N/A'):.4f}")
                elif result['status'] == 'skipped':
                    print(f"⚠ {result['base_name']}: skipped ({result['reason']})")
                else:
                    print(f"✗ {result['base_name']}: error ({result['reason']})")
                    
            except Exception as e:
                print(f"✗ {file_path}: processing failed - {str(e)}")
                results.append({
                    'file_path': file_path,
                    'base_name': os.path.basename(file_path),
                    'status': 'error',
                    'reason': f'executor_error: {str(e)}',
                    'input_count': 0,
                    'kept_count': 0,
                    'removed_count': 0,
                    'removed_sequences': []
                })
    
    # Write filtered files list
    successful_files = []
    with open(args.filtered_files_list, 'w') as f:
        for result in results:
            if result['status'] == 'success' and result.get('output_file') and os.path.exists(result['output_file']):
                f.write(result['output_file'] + '\n')
                successful_files.append(result['output_file'])
    
    # Write summary CSV (file-level statistics)
    with open(args.summary_csv, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow([
            'file_path', 'base_name', 'unweighted_threshold', 'weighted_threshold',
            'input_count', 'kept_count', 'removed_count'
        ])
        
        for result in results:
            thresholds = result.get('thresholds', {})
            writer.writerow([
                result['file_path'],
                result['base_name'],
                thresholds.get('unweighted_threshold', ''),
                thresholds.get('weighted_threshold', ''),
                result['input_count'],
                result['kept_count'],
                result['removed_count']
            ])
    
    # Write detailed metrics CSV (individual sequence removal details)
    with open(args.metrics_csv, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow([
            'base_name', 'sequence_id', 'removal_reason',
            'unweighted_deviation', 'weighted_deviation', 'unweighted_threshold', 'weighted_threshold'
        ])
        
        for result in results:
            # Write individual removed sequences only
            for seq_info in result['removed_sequences']:
                writer.writerow([
                    result['base_name'],
                    seq_info['sequence_id'],
                    seq_info['removal_reason'],
                    f"{seq_info['unweighted_deviation']:.6f}",
                    f"{seq_info['weighted_deviation']:.6f}",
                    f"{seq_info['unweighted_threshold']:.6f}",
                    f"{seq_info['weighted_threshold']:.6f}"
                ])
    
    # Summary statistics
    total_input = sum(r['input_count'] for r in results)
    total_kept = sum(r['kept_count'] for r in results)
    total_removed = sum(r['removed_count'] for r in results)
    successful_count = len([r for r in results if r['status'] == 'success'])
    
    print(f"\nStatistical Outlier Filtering Summary:")
    print(f"Files processed successfully: {successful_count}/{len(input_files)}")
    print(f"Total input sequences: {total_input}")
    print(f"Sequences kept: {total_kept}")
    print(f"Sequences removed (statistical outliers): {total_removed}")
    print(f"Filtered files written: {len(successful_files)}")

if __name__ == "__main__":
    main()
