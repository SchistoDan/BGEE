#!/usr/bin/env python3
"""
===>BWA-Based Reference Sequence Filter<===

A tool for filtering DNA sequences based on their similarity to reference 
sequences using BWA (Burrows-Wheeler Aligner) mapping. This script can either remove 
contaminating sequences that map to known references or retain only sequences that match 
specific reference patterns.

1. For each input FASTA file, extract the sample name to find corresponding reference
2. Create ungapped versions of input sequences (BWA requires ungapped sequences)
3. Build BWA index for the reference sequence (with thread-safe locking)
4. Align ungapped sequences to reference using BWA MEM algorithm
5. Parse alignment results to determine mapped vs unmapped sequences
6. Filter sequences based on the specified mode:
   - 'remove_similar': Keep unmapped sequences (remove contamination/references)
   - 'keep_similar': Keep mapped sequences (retain sequences matching references)
7. Write filtered sequences and detailed metrics


USAGE:
    python bwa_reference_filter.py \\
        --input-files-list input_files.txt \\
        --output-dir filtered_sequences/ \\
        --filtered-files-list filtered_output.txt \\
        --metrics-csv filtering_metrics.csv \\
        --reference-dir reference_sequences/ \\
        --filter-mode remove_similar \\
        --threads 8

Options:
- filter_mode: 
  * 'remove_similar': Removes sequences that map to reference (contamination removal)
  * 'keep_similar': Keeps sequences that map to reference (reference-based retention)
- threads: Number of parallel processes for batch processing
- BWA parameters: Uses default BWA MEM settings with -M flag for multi-mapping detection

File naming:
Input files should follow pattern: {sample_name}_r_{number}_s_{number}_{additional_info}.fasta
- Sample name (e.g., "BSNHM089-24") used to find reference file
- Full basename preserved in output filenames with "_reference_filtered" suffix

Dependencies:
- BWA: Burrows-Wheeler Aligner (must be in PATH)
- samtools: For alignment statistics (must be in PATH)
- BioPython: For FASTA file parsing and writing
- Python 3.6+: For type hints and modern features

AUTHORS: Dan Parsons & Ben Price @ NHMUK.
VERSION: 2.0.0
LICENSE: MIT
"""

import os
import sys
import csv
import argparse
import subprocess
import tempfile
import fcntl
import time
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import List, Tuple, Dict, Optional
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
        
def get_sample_name_for_reference(filepath: str) -> str:
    """Extract just the sample name for finding reference files"""
    basename = os.path.basename(filepath)
    
    # Remove file extensions
    for ext in ['.fasta', '.fas', '.fa']:
        if basename.lower().endswith(ext):
            basename = basename[:-len(ext)]
            break
    
    # Extract everything before the first _r_ occurrence
    if '_r_' in basename:
        return basename.split('_r_')[0]
    
    # Fallback: remove known filter suffixes if no _r_ pattern found
    for suffix in ['_outlier_filtered', '_at_filtered', '_human_filtered', '_align']:
        basename = basename.replace(suffix, '')
    
    return basename

def get_output_basename(filepath: str) -> str:
    """Extract full basename preserving r/s pattern for output filename"""
    basename = os.path.basename(filepath)
    
    # Remove file extensions
    for ext in ['.fasta', '.fas', '.fa']:
        if basename.lower().endswith(ext):
            basename = basename[:-len(ext)]
            break
    
    # Remove only the filter suffix, preserve r/s pattern
    for suffix in ['_outlier_filtered', '_at_filtered', '_human_filtered']:
        if basename.endswith(suffix):
            basename = basename[:-len(suffix)]
            break
    
    return basename

def run_command(cmd: List[str], check: bool = True, capture_output: bool = False) -> subprocess.CompletedProcess:
    """Run a shell command with error handling"""
    try:
        if capture_output:
            result = subprocess.run(cmd, check=check, capture_output=True, text=True)
        else:
            result = subprocess.run(cmd, check=check, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        return result
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"Command failed: {' '.join(cmd)}\nError: {e}")
    except FileNotFoundError:
        raise RuntimeError(f"Command not found: {cmd[0]}. Please ensure BWA and samtools are installed.")

def create_bwa_index_safely(reference_file: str) -> bool:
    """Safely create BWA index with file locking to prevent conflicts"""
    # Check if index already exists
    index_files = [f"{reference_file}{ext}" for ext in ['.amb', '.ann', '.bwt', '.pac', '.sa']]
    if all(os.path.exists(f) for f in index_files):
        return True  # Index already exists
    
    # Use lock file to prevent simultaneous indexing
    lock_file = f"{reference_file}.bwa_index.lock"
    
    try:
        # Try to acquire lock
        with open(lock_file, 'w') as lock_fd:
            # Try to get exclusive lock (non-blocking)
            try:
                fcntl.flock(lock_fd.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
                
                # We got the lock - check again if index exists (another process might have created it)
                if all(os.path.exists(f) for f in index_files):
                    return True
                
                # Create the index
                run_command(['bwa', 'index', reference_file])
                
                # Verify index was created successfully
                if all(os.path.exists(f) for f in index_files):
                    return True
                else:
                    raise RuntimeError(f"BWA index files were not created properly for {reference_file}")
                    
            except BlockingIOError:
                # Another process is creating the index - wait for it
                print(f"Waiting for BWA index creation for {os.path.basename(reference_file)}...")
                
                # Wait with timeout
                max_wait = 300  # 5 minutes
                wait_interval = 5
                waited = 0
                
                while waited < max_wait:
                    time.sleep(wait_interval)
                    waited += wait_interval
                    
                    # Check if index is ready
                    if all(os.path.exists(f) for f in index_files):
                        return True
                        
                # Timeout
                raise RuntimeError(f"Timeout waiting for BWA index creation for {reference_file}")
                
    except Exception as e:
        raise RuntimeError(f"Failed to create BWA index for {reference_file}: {str(e)}")
    
    finally:
        # Clean up lock file if we created it
        try:
            if os.path.exists(lock_file):
                os.remove(lock_file)
        except:
            pass
    
    return False

def create_ungapped_fasta(input_file: str, output_file: str) -> int:
    """Remove gaps from FASTA sequences and write to output file"""
    sequence_count = 0
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for record in SeqIO.parse(infile, 'fasta'):
            ungapped_seq = str(record.seq).replace('-', '')
            if ungapped_seq:  # Only write sequences with content
                outfile.write(f">{record.id}\n{ungapped_seq}\n")
                sequence_count += 1
    return sequence_count

def parse_sam_file(sam_file: str) -> Dict[str, bool]:
    """Parse SAM file and extract unmapped status for each sequence"""
    mapping_info = {}
    
    with open(sam_file, 'r') as f:
        for line in f:
            if line.startswith('@'):
                continue  # Skip header lines
            
            fields = line.strip().split('\t')
            if len(fields) < 11:
                continue
                
            seq_id = fields[0]
            flag = int(fields[1])
            
            # Determine if sequence is unmapped
            is_unmapped = bool(flag & 4)  # Check unmapped flag
            mapping_info[seq_id] = is_unmapped
    
    return mapping_info

def get_flagstat_info(sam_file: str) -> Dict[str, int]:
    """Get mapping statistics using samtools flagstat"""
    try:
        result = run_command(['samtools', 'flagstat', sam_file], capture_output=True)
        lines = result.stdout.strip().split('\n')
        
        stats = {}
        for line in lines:
            if 'in total' in line:
                stats['total_reads'] = int(line.split()[0])
            elif 'mapped (' in line and 'primary' not in line:
                stats['mapped_reads'] = int(line.split()[0])
        
        stats['unmapped_reads'] = stats.get('total_reads', 0) - stats.get('mapped_reads', 0)
        return stats
    except Exception:
        return {'total_reads': 0, 'mapped_reads': 0, 'unmapped_reads': 0}

def process_single_file(file_path: str, reference_dir: str, filter_mode: str,
                       output_dir: str, bwa_params: Dict) -> Dict:
    """Process a single FASTA file for BWA-based reference filtering"""
    try:
        # Get sample name for finding reference (e.g., "BSNHM089-24")
        sample_name = get_sample_name_for_reference(file_path)
        
        # Get full basename for output filename (e.g., "BSNHM089-24_r_1_s_50_BSNHM089-24")
        output_basename = get_output_basename(file_path)
        
        # Check if file is empty
        if os.path.getsize(file_path) == 0:
            return {
                'file_path': file_path,
                'base_name': sample_name,
                'sample_name': sample_name,
                'output_basename': output_basename,
                'status': 'skipped',
                'reason': 'empty_file',
                'input_count': 0,
                'kept_count': 0,
                'removed_count': 0,
                'removed_sequences': []
            }
        
        # Find reference file using sample name
        reference_file = os.path.join(reference_dir, f"{sample_name}_reference.fasta")
        if not os.path.exists(reference_file):
            return {
                'file_path': file_path,
                'base_name': sample_name,
                'sample_name': sample_name,
                'output_basename': output_basename,
                'status': 'skipped',
                'reason': 'no_reference_file',
                'input_count': 0,
                'kept_count': 0,
                'removed_count': 0,
                'removed_sequences': []
            }
        
        # Read input sequences to get total count
        try:
            records = list(SeqIO.parse(file_path, "fasta"))
        except Exception as e:
            return {
                'file_path': file_path,
                'base_name': sample_name,
                'sample_name': sample_name,
                'output_basename': output_basename,
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
                'base_name': sample_name,
                'sample_name': sample_name,
                'output_basename': output_basename,
                'status': 'skipped',
                'reason': 'no_sequences',
                'input_count': 0,
                'kept_count': 0,
                'removed_count': 0,
                'removed_sequences': []
            }
        
        input_count = len(records)
        
        # Create temporary files
        with tempfile.TemporaryDirectory() as temp_dir:
            ungapped_file = os.path.join(temp_dir, 'sequences_ungapped.fasta')
            sam_file = os.path.join(temp_dir, 'alignments.sam')
            
            # Step 1: Create ungapped sequences for BWA
            ungapped_count = create_ungapped_fasta(file_path, ungapped_file)
            if ungapped_count == 0:
                return {
                    'file_path': file_path,
                    'base_name': sample_name,
                    'sample_name': sample_name,
                    'output_basename': output_basename,
                    'status': 'skipped',
                    'reason': 'no_ungapped_sequences',
                    'input_count': input_count,
                    'kept_count': 0,
                    'removed_count': 0,
                    'removed_sequences': []
                }
            
            # Step 2: Create BWA index safely (with locking to prevent conflicts)
            try:
                if not create_bwa_index_safely(reference_file):
                    return {
                        'file_path': file_path,
                        'base_name': sample_name,
                        'sample_name': sample_name,
                        'output_basename': output_basename,
                        'status': 'error',
                        'reason': 'bwa_index_creation_failed',
                        'input_count': input_count,
                        'kept_count': 0,
                        'removed_count': 0,
                        'removed_sequences': []
                    }
            except Exception as e:
                return {
                    'file_path': file_path,
                    'base_name': sample_name,
                    'sample_name': sample_name,
                    'output_basename': output_basename,
                    'status': 'error',
                    'reason': f'bwa_index_error: {str(e)}',
                    'input_count': input_count,
                    'kept_count': 0,
                    'removed_count': 0,
                    'removed_sequences': []
                }
            
            # Step 3: Run BWA mem alignment
            try:
                bwa_cmd = ['bwa', 'mem']
                
                # Add BWA parameters
                for param, value in bwa_params.items():
                    if value is not None:
                        if isinstance(value, bool) and value:
                            bwa_cmd.append(f'-{param}')
                        else:
                            bwa_cmd.extend([f'-{param}', str(value)])
                
                bwa_cmd.extend([reference_file, ungapped_file])
                
                with open(sam_file, 'w') as sam_out:
                    result = subprocess.run(bwa_cmd, stdout=sam_out, stderr=subprocess.DEVNULL, check=True)
                    
            except Exception as e:
                return {
                    'file_path': file_path,
                    'base_name': sample_name,
                    'sample_name': sample_name,
                    'output_basename': output_basename,
                    'status': 'error',
                    'reason': f'bwa_mem_error: {str(e)}',
                    'input_count': input_count,
                    'kept_count': 0,
                    'removed_count': 0,
                    'removed_sequences': []
                }
            
            # Step 4: Parse SAM file to get mapping information
            mapping_info = parse_sam_file(sam_file)
            flagstat_info = get_flagstat_info(sam_file)
            
            # Step 5: Determine which sequences to keep based on filter mode
            kept_records = []
            removed_sequences = []
            
            for record in records:
                is_unmapped = mapping_info.get(record.id, True)  # Default to unmapped if not found
                
                # Determine if sequence should be removed based on filter mode
                if filter_mode == 'remove_similar':
                    # Remove mapped sequences (keep unmapped)
                    is_removed = not is_unmapped
                    removal_reason = 'reference_mapped' if is_removed else None
                else:  # keep_similar
                    # Remove unmapped sequences (keep mapped)  
                    is_removed = is_unmapped
                    removal_reason = 'reference_unmapped' if is_removed else None
                
                if is_removed:
                    removed_sequences.append({
                        'sequence_id': record.id,
                        'removal_reason': removal_reason
                    })
                else:
                    kept_records.append(record)
            
            # Clean up BWA index files only if they didn't exist before
            # (Don't remove shared index files that other processes might be using)
            # The index files will be reused for other samples with the same reference
        
        # Step 6: Write filtered sequences using full basename
        output_file = None
        if kept_records:
            output_file = os.path.join(output_dir, f"{output_basename}_reference_filtered.fasta")
            with open(output_file, 'w') as handle:
                SeqIO.write(kept_records, handle, "fasta")
        
        return {
            'file_path': file_path,
            'base_name': sample_name,
            'sample_name': sample_name,
            'output_basename': output_basename,
            'status': 'success',
            'reason': 'processed',
            'input_count': input_count,
            'kept_count': len(kept_records),
            'removed_count': len(removed_sequences),
            'removed_sequences': removed_sequences,
            'output_file': output_file,
            'reference_file': reference_file,
            'filter_mode': filter_mode,
            'flagstat_info': flagstat_info,
            'bwa_params': bwa_params
        }
        
    except Exception as e:
        sample_name = get_sample_name_for_reference(file_path)
        output_basename = get_output_basename(file_path)
        return {
            'file_path': file_path,
            'base_name': sample_name,
            'sample_name': sample_name,
            'output_basename': output_basename,
            'status': 'error',
            'reason': f'processing_error: {str(e)}',
            'input_count': 0,
            'kept_count': 0,
            'removed_count': 0,
            'removed_sequences': []
        }

def main():
    parser = argparse.ArgumentParser(description='Filter sequences using BWA alignment to reference')
    parser.add_argument('--input-files-list', required=True, help='File containing list of FASTA files to process')
    parser.add_argument('--output-dir', required=True, help='Output directory for filtered files')
    parser.add_argument('--filtered-files-list', required=True, help='Output file listing successfully filtered files')
    parser.add_argument('--metrics-csv', required=True, help='Output CSV file with filtering metrics')
    parser.add_argument('--reference-dir', required=True, help='Directory containing reference sequences')
    parser.add_argument('--filter-mode', required=True, choices=['keep_similar', 'remove_similar'],
                       help='Filter mode: keep_similar (keep mapped sequences) or remove_similar (keep unmapped sequences)')
    parser.add_argument('--threads', type=int, default=1, help='Number of threads for parallel processing')
    
    # BWA parameters - using defaults except for -M which helps with multi-mapping contamination detection
    
    args = parser.parse_args()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Validate reference directory
    if not os.path.exists(args.reference_dir):
        print(f"Error: Reference directory does not exist: {args.reference_dir}")
        sys.exit(1)
    
    # Check BWA and samtools availability
    try:
        run_command(['bwa'], check=False)
        run_command(['samtools', '--version'], check=False)
    except RuntimeError as e:
        print(f"Error: {e}")
        sys.exit(1)
    
    # Read input files
    with open(args.input_files_list, 'r') as f:
        input_files = [line.strip() for line in f if line.strip()]
    
    # Prepare BWA parameters - using defaults plus -M for multi-mapping reads
    bwa_params = {
        'M': True  # Mark shorter split hits as secondary (helps catch all contamination mappings)
    }
    
    # Print configuration
    mode_description = {
        'keep_similar': 'keeping sequences that map to reference (removing unmapped/divergent)',
        'remove_similar': 'removing sequences that map to reference (keeping unmapped/divergent)'
    }
    
    print(f"Processing {len(input_files)} files with BWA-based reference filtering")
    print(f"Filter mode: {args.filter_mode} ({mode_description[args.filter_mode]})")
    print(f"BWA: using default parameters with -M flag for multi-mapping detection")
    
    # Process files in parallel
    results = []
    with ProcessPoolExecutor(max_workers=args.threads) as executor:
        future_to_file = {
            executor.submit(process_single_file, file_path, args.reference_dir, 
                          args.filter_mode, args.output_dir, bwa_params): file_path
            for file_path in input_files
        }
        
        for future in as_completed(future_to_file):
            file_path = future_to_file[future]
            try:
                result = future.result()
                results.append(result)
                
                if result['status'] == 'success':
                    flagstat = result.get('flagstat_info', {})
                    print(f"✓ {result['base_name']}: {result['kept_count']}/{result['input_count']} sequences kept")
                    print(f"  Reference: {os.path.basename(result.get('reference_file', 'N/A'))}")
                    print(f"  Mode: {result.get('filter_mode', 'N/A')}")
                    print(f"  Mapping stats: {flagstat.get('mapped_reads', 0)}/{flagstat.get('total_reads', 0)} mapped")
                elif result['status'] == 'skipped':
                    print(f"⚠ {result['base_name']}: skipped ({result['reason']})")
                else:
                    print(f"✗ {result['base_name']}: error ({result['reason']})")
                    
            except Exception as e:
                print(f"✗ {file_path}: processing failed - {str(e)}")
                sample_name = get_sample_name_for_reference(file_path)
                output_basename = get_output_basename(file_path)
                results.append({
                    'file_path': file_path,
                    'base_name': sample_name,
                    'sample_name': sample_name,
                    'output_basename': output_basename,
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
        writer.writerow([
            'file_path', 'base_name', 'sequence_id', 'removal_reason', 
            'reference_file', 'filter_mode', 'step_name', 'input_count', 'kept_count', 'removed_count'
        ])
        
        for result in results:
            # Write file-level summary
            writer.writerow([
                result['file_path'],
                result['base_name'],
                'FILE_SUMMARY',
                result['reason'],
                result.get('reference_file', ''),
                result.get('filter_mode', ''),
                'reference_filter',
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
                    result.get('reference_file', ''),
                    result.get('filter_mode', ''),
                    'reference_filter',
                    '',
                    '',
                    ''
                ])
    
    # Summary statistics
    total_input = sum(r['input_count'] for r in results)
    total_kept = sum(r['kept_count'] for r in results)
    total_removed = sum(r['removed_count'] for r in results)
    successful_count = len([r for r in results if r['status'] == 'success'])
    skipped_no_ref = len([r for r in results if r['status'] == 'skipped' and r['reason'] == 'no_reference_file'])
    
    print(f"\nBWA Reference Filtering Summary:")
    print(f"Filter mode: {args.filter_mode}")
    print(f"Files processed successfully: {successful_count}/{len(input_files)}")
    print(f"Files skipped (no reference): {skipped_no_ref}")
    print(f"Total input sequences: {total_input}")
    print(f"Sequences kept: {total_kept}")
    print(f"Sequences removed: {total_removed}")
    print(f"Filtered files written: {len(successful_files)}")

if __name__ == "__main__":
    main()
