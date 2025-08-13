#!/usr/bin/env python3
"""
Consensus Generator
Generates consensus sequences from filtered alignments and creates final outputs
Now includes coverage statistics calculation and no line wrapping for FASTA output
"""

import os
import sys
import csv
import argparse
import statistics
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.SeqIO import FastaIO
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import List, Tuple, Dict
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

def calculate_coverage_statistics(sequences: List[str]) -> Dict[str, float]:
    """Calculate coverage statistics for alignment sequences"""
    if not sequences:
        return {
            'cleaning_cov_percent': 0.0,
            'cleaning_cov_avg': 0.0,
            'cleaning_cov_med': 0.0,
            'cleaning_cov_max': 0,
            'cleaning_cov_min': 0
        }
    
    # Get maximum sequence length
    max_len = max(len(seq) for seq in sequences)
    
    # Calculate coverage depth at each position
    coverage_depths = []
    positions_with_coverage = 0
    
    for i in range(max_len):
        coverage_at_position = 0
        for seq in sequences:
            if i < len(seq) and seq[i] != '-':
                coverage_at_position += 1
        
        coverage_depths.append(coverage_at_position)
        if coverage_at_position > 0:
            positions_with_coverage += 1
    
    # Calculate statistics
    if coverage_depths:
        cleaning_cov_min = min(coverage_depths)
        cleaning_cov_max = max(coverage_depths)
        cleaning_cov_avg = sum(coverage_depths) / len(coverage_depths)
        cleaning_cov_med = statistics.median(coverage_depths)
        cleaning_cov_percent = (positions_with_coverage / len(coverage_depths)) * 100.0
    else:
        cleaning_cov_min = 0
        cleaning_cov_max = 0
        cleaning_cov_avg = 0.0
        cleaning_cov_med = 0.0
        cleaning_cov_percent = 0.0
    
    return {
        'cleaning_cov_percent': cleaning_cov_percent,
        'cleaning_cov_avg': cleaning_cov_avg,
        'cleaning_cov_med': cleaning_cov_med,
        'cleaning_cov_max': cleaning_cov_max,
        'cleaning_cov_min': cleaning_cov_min
    }

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

def calculate_at_content(sequence: str) -> float:
    """Calculate AT content for a sequence"""
    if not sequence:
        return 0.0
    sequence_clean = sequence.replace('-', '').upper()
    if not sequence_clean:
        return 0.0
    at_count = sequence_clean.count('A') + sequence_clean.count('T')
    return at_count / len(sequence_clean)

def count_ambiguous_bases(consensus_sequence: str) -> int:
    """Count the number of ambiguous bases (not G, T, A, or C) in a consensus sequence"""
    sequence = consensus_sequence.upper().replace('-', '')
    non_standard_bases = sum(1 for base in sequence if base not in 'GTAC')
    return non_standard_bases

def get_base_filename(filepath: str) -> str:
    """Extract base filename without extension and filter suffixes"""
    basename = os.path.basename(filepath)
    for ext in ['.fasta', '.fas', '.fa']:
        if basename.lower().endswith(ext):
            basename = basename[:-len(ext)]
            break
    
    # Remove filter suffixes
    for suffix in ['_reference_filtered', '_outlier_filtered', '_at_filtered', '_human_filtered', '_align']:
        basename = basename.replace(suffix, '')
    
    return basename

def process_single_file(file_path: str, consensus_threshold: float, preprocessing_mode: str, output_dir: str) -> Dict:
    """Process a single filtered FASTA file to generate consensus and metrics"""
    try:
        base_name = get_base_filename(file_path)
        
        # Check if file is empty
        if os.path.getsize(file_path) == 0:
            return {
                'base_name': base_name,
                'status': 'skipped',
                'reason': 'empty_file',
                'input_sequences': 0,
                'consensus_sequence': '',
                'consensus_at_content': 0.0,
                'ambiguous_bases': 0,
                'consensus_length': 0,
                'cleaning_cov_percent': 0.0,
                'cleaning_cov_avg': 0.0,
                'cleaning_cov_med': 0.0,
                'cleaning_cov_max': 0,
                'cleaning_cov_min': 0
            }
        
        # Read sequences
        try:
            records = list(SeqIO.parse(file_path, "fasta"))
        except Exception as e:
            return {
                'base_name': base_name,
                'status': 'error',
                'reason': f'parse_error: {str(e)}',
                'input_sequences': 0,
                'consensus_sequence': '',
                'consensus_at_content': 0.0,
                'ambiguous_bases': 0,
                'consensus_length': 0,
                'cleaning_cov_percent': 0.0,
                'cleaning_cov_avg': 0.0,
                'cleaning_cov_med': 0.0,
                'cleaning_cov_max': 0,
                'cleaning_cov_min': 0
            }
        
        if not records:
            return {
                'base_name': base_name,
                'status': 'skipped',
                'reason': 'no_sequences',
                'input_sequences': 0,
                'consensus_sequence': '',
                'consensus_at_content': 0.0,
                'ambiguous_bases': 0,
                'consensus_length': 0,
                'cleaning_cov_percent': 0.0,
                'cleaning_cov_avg': 0.0,
                'cleaning_cov_med': 0.0,
                'cleaning_cov_max': 0,
                'cleaning_cov_min': 0
            }
        
        input_count = len(records)
        
        # Get sequences as strings
        sequences = [str(record.seq).upper() for record in records]
        
        # Calculate coverage statistics
        coverage_stats = calculate_coverage_statistics(sequences)
        
        # Generate consensus sequence
        consensus_seq, frequencies = generate_consensus_sequence(sequences, consensus_threshold)
        
        if not consensus_seq:
            return {
                'base_name': base_name,
                'status': 'error',
                'reason': 'consensus_generation_failed',
                'input_sequences': input_count,
                'consensus_sequence': '',
                'consensus_at_content': 0.0,
                'ambiguous_bases': 0,
                'consensus_length': 0,
                **coverage_stats
            }
        
        # Calculate consensus metrics
        consensus_at_content = calculate_at_content(consensus_seq)
        ambiguous_bases = count_ambiguous_bases(consensus_seq)
        consensus_length = len(consensus_seq.replace('-', ''))
        
        # Create consensus record with appropriate naming
        consensus_id = f"{base_name}_fcleaner"
        if preprocessing_mode == 'merge':
            consensus_id += "_merge"
        
        consensus_record = SeqRecord(
            Seq(consensus_seq),
            id=consensus_id,
            description=""
        )
        
        # Write individual consensus file directly to output_dir (without wrapping)
        consensus_file = os.path.join(output_dir, f"{base_name}_fcleaner.fasta")
        with open(consensus_file, 'w') as handle:
            writer = FastaIO.FastaWriter(handle, wrap=None)
            writer.write_file([consensus_record])
        
        # Write cleaned sequences (without wrapping)
        cleaned_output_dir = os.path.join(output_dir, "filter_pass_seqs")
        os.makedirs(cleaned_output_dir, exist_ok=True)
        cleaned_reads_file = os.path.join(cleaned_output_dir, f"{base_name}_cleaned.fasta")
        with open(cleaned_reads_file, 'w') as handle:
            writer = FastaIO.FastaWriter(handle, wrap=None)
            writer.write_file(records)
        
        return {
            'base_name': base_name,
            'status': 'success',
            'reason': 'processed',
            'input_sequences': input_count,
            'consensus_sequence': consensus_seq,
            'consensus_at_content': consensus_at_content,
            'ambiguous_bases': ambiguous_bases,
            'consensus_length': consensus_length,
            'consensus_record': consensus_record,
            'consensus_file': consensus_file,
            'cleaned_reads_file': cleaned_reads_file,
            **coverage_stats
        }
        
    except Exception as e:
        return {
            'base_name': get_base_filename(file_path),
            'status': 'error',
            'reason': f'processing_error: {str(e)}',
            'input_sequences': 0,
            'consensus_sequence': '',
            'consensus_at_content': 0.0,
            'ambiguous_bases': 0,
            'consensus_length': 0,
            'cleaning_cov_percent': 0.0,
            'cleaning_cov_avg': 0.0,
            'cleaning_cov_med': 0.0,
            'cleaning_cov_max': 0,
            'cleaning_cov_min': 0
        }

def main():
    parser = argparse.ArgumentParser(description='Generate consensus sequences from filtered alignments')
    parser.add_argument('--input-files-list', required=True, help='File containing list of filtered FASTA files')
    parser.add_argument('--output-dir', required=True, help='Base output directory')
    parser.add_argument('--consensus-fasta', required=True, help='Output concatenated consensus FASTA file')
    parser.add_argument('--consensus-metrics', required=True, help='Filename for consensus metrics CSV (will be placed in output directory)')
    parser.add_argument('--consensus-threshold', type=float, default=0.5, help='Consensus generation threshold')
    parser.add_argument('--preprocessing-mode', choices=['merge', 'concat'], default='concat', 
                       help='Preprocessing mode for naming')
    parser.add_argument('--threads', type=int, default=1, help='Number of threads for parallel processing')
    
    args = parser.parse_args()
    
    # Create output directories early
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Ensure consensus FASTA output directory exists
    os.makedirs(os.path.dirname(args.consensus_fasta), exist_ok=True)
    
    # Set consensus metrics path to be in main output directory
    consensus_metrics_path = os.path.join(args.output_dir, os.path.basename(args.consensus_metrics))
    
    # Read input files
    with open(args.input_files_list, 'r') as f:
        input_files = [line.strip() for line in f if line.strip()]
    
    print(f"Generating consensus sequences from {len(input_files)} filtered files")
    print(f"Consensus threshold: {args.consensus_threshold}")
    print(f"Preprocessing mode: {args.preprocessing_mode}")
    print(f"Consensus metrics will be written to: {consensus_metrics_path}")
    
    # Process files in parallel
    results = []
    with ProcessPoolExecutor(max_workers=args.threads) as executor:
        future_to_file = {
            executor.submit(process_single_file, file_path, args.consensus_threshold, 
                          args.preprocessing_mode, args.output_dir): file_path
            for file_path in input_files
        }
        
        for future in as_completed(future_to_file):
            file_path = future_to_file[future]
            try:
                result = future.result()
                results.append(result)
                
                if result['status'] == 'success':
                    print(f"✓ {result['base_name']}: consensus generated from {result['input_sequences']} sequences")
                    print(f"  AT content: {result['consensus_at_content']:.3f}, Ambiguous bases: {result['ambiguous_bases']}")
                    print(f"  Coverage - avg: {result['cleaning_cov_avg']:.1f}, med: {result['cleaning_cov_med']:.1f}, range: {result['cleaning_cov_min']}-{result['cleaning_cov_max']}")
                elif result['status'] == 'skipped':
                    print(f"⚠ {result['base_name']}: skipped ({result['reason']})")
                else:
                    print(f"✗ {result['base_name']}: error ({result['reason']})")
                    
            except Exception as e:
                print(f"✗ {file_path}: processing failed - {str(e)}")
                results.append({
                    'base_name': os.path.basename(file_path),
                    'status': 'error',
                    'reason': f'executor_error: {str(e)}',
                    'input_sequences': 0,
                    'consensus_sequence': '',
                    'consensus_at_content': 0.0,
                    'ambiguous_bases': 0,
                    'consensus_length': 0,
                    'cleaning_cov_percent': 0.0,
                    'cleaning_cov_avg': 0.0,
                    'cleaning_cov_med': 0.0,
                    'cleaning_cov_max': 0,
                    'cleaning_cov_min': 0
                })
    
    # Collect successful consensus records
    successful_results = [r for r in results if r['status'] == 'success']
    
    if successful_results:
        # Write concatenated consensus file (without wrapping)
        consensus_records = [r['consensus_record'] for r in successful_results]
        with open(args.consensus_fasta, 'w') as handle:
            writer = FastaIO.FastaWriter(handle, wrap=None)
            writer.write_file(consensus_records)
        
        print(f"Concatenated consensus sequences written to: {args.consensus_fasta}")
    else:
        # Create empty consensus file
        open(args.consensus_fasta, 'w').close()
        print("No successful consensus sequences generated - empty file created")
    
    # Write consensus metrics CSV to consensus_seqs directory
    with open(consensus_metrics_path, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow([
            'sample_name', 'status', 'cleaned_reads', 
            'consensus_at_content', 'ambiguous_bases', 'consensus_length',
            'cov_percent', 'cov_avg', 'cov_med', 
            'cov_max', 'cov_min',
            'consensus_file', 'cleaned_reads_file'
        ])
        
        for result in results:
            # Create sample_name with _fcleaner suffix
            sample_name = f"{result['base_name']}_fcleaner"
            if args.preprocessing_mode == 'merge':
                sample_name += "_merge"
                
            writer.writerow([
                sample_name,
                result['status'],
                result['input_sequences'],
                f"{result['consensus_at_content']:.6f}",
                result['ambiguous_bases'],
                result['consensus_length'],
                f"{result['cleaning_cov_percent']:.2f}",
                f"{result['cleaning_cov_avg']:.2f}",
                f"{result['cleaning_cov_med']:.2f}",
                result['cleaning_cov_max'],
                result['cleaning_cov_min'],
                result.get('consensus_file', ''),
                result.get('cleaned_reads_file', '')
            ])
    
    print(f"Consensus metrics written to: {consensus_metrics_path}")
    
    # Summary statistics
    total_input_sequences = sum(r['input_sequences'] for r in results)
    successful_count = len(successful_results)
    total_ambiguous = sum(r['ambiguous_bases'] for r in successful_results)
    
    if successful_results:
        avg_at_content = sum(r['consensus_at_content'] for r in successful_results) / len(successful_results)
        avg_ambiguous = total_ambiguous / len(successful_results)
        avg_coverage = sum(r['cleaning_cov_avg'] for r in successful_results) / len(successful_results)
    else:
        avg_at_content = 0.0
        avg_ambiguous = 0.0
        avg_coverage = 0.0
    
    print(f"\nConsensus Generation Summary:")
    print(f"Files processed successfully: {successful_count}/{len(input_files)}")
    print(f"Total input sequences: {total_input_sequences}")
    print(f"Consensus sequences generated: {successful_count}")
    print(f"Average AT content: {avg_at_content:.3f}")
    print(f"Average ambiguous bases per consensus: {avg_ambiguous:.1f}")
    print(f"Average coverage depth: {avg_coverage:.1f}")
    print(f"Total ambiguous bases: {total_ambiguous}")

if __name__ == "__main__":
    main()
