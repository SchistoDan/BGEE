#!/usr/bin/env python3
"""
Human COX1 Similarity Filter
Removes sequences with high similarity to human mitochondrial COX1
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

# Human COX1 reference sequence
HUMAN_COX1 = """ATGTTCGCCGACCGTTGACTATTCTCTACAAACCACAAAGACATTGGAACACTATACCTATTATTCGGCG
CATGAGCTGGAGTCCTAGGCACAGCTCTAAGCCTCCTTATTCGAGCCGAGCTGGGCCAGCCAGGCAACCT
TCTAGGTAACGACCACATCTACAACGTTATCGTCACAGCCCATGCATTTGTAATAATCTTCTTCATAGTA
ATACCCATCATAATCGGAGGCTTTGGCAACTGACTAGTTCCCCTAATAATCGGTGCCCCCGATATGGCGT
TTCCCCGCATAAACAACATAAGCTTCTGACTCTTACCTCCCTCTCTCCTACTCCTGCTCGCATCTGCTAT
AGTGGAGGCCGGAGCAGGAACAGGTTGAACAGTCTACCCTCCCTTAGCAGGGAACTACTCCCACCCTGGA
GCCTCCGTAGACCTAACCATCTTCTCCTTACACCTAGCAGGTGTCTCCTCTATCTTAGGGGCCATCAATT
TCATCACAACAATTATCAATATAAAACCCCCTGCCATAACCCAATACCAAACGCCCCTCTTCGTCTGATC
CGTCCTAATCACAGCAGTCCTACTTCTCCTATCTCTCCCAGTCCTAGCTGCTGGCATCACTATACTACTA
ACAGACCGCAACCTCAACACCACCTTCTTCGACCCCGCCGGAGGAGGAGACCCCATTCTATACCAACACC
TATTCTGATTTTTCGGTCACCCTGAAGTTTATATTCTTATCCTACCAGGCTTCGGAATAATCTCCCATAT
TGTAACTTACTACTCCGGAAAAAAAGAACCATTTGGATACATAGGTATGGTCTGAGCTATGATATCAATT
GGCTTCCTAGGGTTTATCGTGTGAGCACACCATATATTTACAGTAGGAATAGACGTAGACACACGAGCAT
ATTTCACCTCCGCTACCATAATCATCGCTATCCCCACCGGCGTCAAAGTATTTAGCTGACTCGCCACACT
CCACGGAAGCAATATGAAATGATCTGCTGCAGTGCTCTGAGCCCTAGGATTCATCTTTCTTTTCACCGTA
GGTGGCCTGACTGGCATTGTATTAGCAAACTCATCACTAGACATCGTACTACACGACACGTACTACGTTG
TAGCCCACTTCCACTATGTCCTATCAATAGGAGCTGTATTTGCCATCATAGGAGGCTTCATTCACTGATT
TCCCCTATTCTCAGGCTACACCCTAGACCAAACCTACGCCAAAATCCATTTCACTATCATATTCATCGGC
GTAAATCTAACTTTCTTCCCACAACACTTTCTCGGCCTATCCGGAATGCCCCGACGTTACTCGGACTACC
CCGATGCATACACCACATGAAACATCCTATCATCTGTAGGCTCATTCATTTCTCTAACAGCAGTAATATT
AATAATTTTCATGATTTGAGAAGCCTTCGCTTCGAAGCGAAAAGTCCTAATAGTAGAAGAACCCTCCATA
AACCTGGAGTGACTATATGGATGCCCCCCACCCTACCACACATTCGAAGAACCCGTATACATAAAATCTA
GA""".replace("\n", "")

HUMAN_COX1_ARRAY = np.array(list(HUMAN_COX1.upper()))

def log_message(message: str, log_file=None, stdout=False):
    """Log message to file and optionally stdout"""
    timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    formatted_msg = f"[{timestamp}] {message}"
    
    if log_file:
        log_file.write(formatted_msg + '\n')
        log_file.flush()
    
    if stdout:
        print(formatted_msg, flush=True)

def vectorised_sequence_similarity(seq_no_gaps: str, reference_array: np.ndarray, min_overlap: int = 20) -> float:
    """Calculate sequence similarity using vectorised operations"""
    if not seq_no_gaps or len(seq_no_gaps) == 0:
        return 0.0
    
    seq_array = np.array(list(seq_no_gaps.upper()))
    ref_len = min(len(seq_array), len(reference_array))
    
    if ref_len < min_overlap:
        return 0.0
    
    seq_part = seq_array[:ref_len]
    ref_part = reference_array[:ref_len]
    
    matches = np.sum(seq_part == ref_part)
    similarity = matches / ref_len
    
    return similarity

def process_single_file(file_path: str, human_threshold: float, output_dir: str) -> Dict:
    """Process a single FASTA file for human COX1 filtering"""
    try:
        base_name = os.path.splitext(os.path.basename(file_path))[0]
        if '_align' in base_name:
            base_name = base_name.replace('_align', '')
        
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
        kept_records = []
        removed_sequences = []
        
        # Filter sequences
        for record in records:
            sequence_no_gaps = str(record.seq).replace('-', '').upper()
            
            if len(sequence_no_gaps) == 0:
                # Empty sequence - remove
                removed_sequences.append({
                    'sequence_id': record.id,
                    'removal_reason': 'empty_sequence',
                    'human_similarity': 0.0
                })
                continue
            
            # Calculate human similarity
            similarity = vectorised_sequence_similarity(sequence_no_gaps, HUMAN_COX1_ARRAY)
            
            if similarity >= human_threshold:
                removed_sequences.append({
                    'sequence_id': record.id,
                    'removal_reason': 'human_similar',
                    'human_similarity': similarity
                })
            else:
                kept_records.append(record)
        
        # Write filtered sequences
        if kept_records:
            output_file = os.path.join(output_dir, f"{base_name}_human_filtered.fasta")
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
            'output_file': os.path.join(output_dir, f"{base_name}_human_filtered.fasta") if kept_records else None
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
    parser = argparse.ArgumentParser(description='Filter sequences based on human COX1 similarity')
    parser.add_argument('--input-log', required=True, help='File containing list of FASTA alignment files')
    parser.add_argument('--output-dir', required=True, help='Output directory for filtered files')
    parser.add_argument('--filtered-files-list', required=True, help='Output file listing successfully filtered files')
    parser.add_argument('--metrics-csv', required=True, help='Output CSV file with filtering metrics')
    parser.add_argument('--human-threshold', type=float, default=0.95, help='Human COX1 similarity threshold')
    parser.add_argument('--threads', type=int, default=1, help='Number of threads for parallel processing')
    
    args = parser.parse_args()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Read input files
    with open(args.input_log, 'r') as f:
        input_files = [line.strip() for line in f if line.strip()]
    
    print(f"Processing {len(input_files)} files with human COX1 similarity threshold {args.human_threshold}")
    
    # Process files in parallel
    results = []
    with ProcessPoolExecutor(max_workers=args.threads) as executor:
        future_to_file = {
            executor.submit(process_single_file, file_path, args.human_threshold, args.output_dir): file_path
            for file_path in input_files
        }
        
        for future in as_completed(future_to_file):
            file_path = future_to_file[future]
            try:
                result = future.result()
                results.append(result)
                
                if result['status'] == 'success':
                    print(f"✓ {result['base_name']}: {result['kept_count']}/{result['input_count']} sequences kept")
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
    
    # Write detailed metrics CSV
    with open(args.metrics_csv, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['file_path', 'base_name', 'sequence_id', 'removal_reason', 'human_similarity', 'step_name', 'input_count', 'kept_count', 'removed_count'])
        
        for result in results:
            # Write file-level summary
            writer.writerow([
                result['file_path'],
                result['base_name'],
                'FILE_SUMMARY',
                result['reason'],
                '',
                'human_cox1_filter',
                result['input_count'],
                result['kept_count'],
                result['removed_count']
            ])
            
            # Write individual removed sequences
            for seq_info in result['removed_sequences']:
                writer.writerow([
                    result['file_path'],
                    result['base_name'],
                    seq_info['sequence_id'],
                    seq_info['removal_reason'],
                    f"{seq_info['human_similarity']:.4f}",
                    'human_cox1_filter',
                    '',
                    '',
                    ''
                ])
    
    # Summary statistics
    total_input = sum(r['input_count'] for r in results)
    total_kept = sum(r['kept_count'] for r in results)
    total_removed = sum(r['removed_count'] for r in results)
    successful_count = len([r for r in results if r['status'] == 'success'])
    
    print(f"\nHuman COX1 Filtering Summary:")
    print(f"Files processed successfully: {successful_count}/{len(input_files)}")
    print(f"Total input sequences: {total_input}")
    print(f"Sequences kept: {total_kept}")
    print(f"Sequences removed (human similar): {total_removed}")
    print(f"Filtered files written: {len(successful_files)}")

if __name__ == "__main__":
    main()
