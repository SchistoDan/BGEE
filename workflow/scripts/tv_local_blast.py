#!/usr/bin/env python3
"""
Local BLASTn (parallelised)

This script provides parallel processing capabilities for running BLASTn searches on FASTA files.
It can handle single sequences, multi-FASTA files, or entire directories of FASTA files, and
automatically creates BLAST databases from FASTA files when needed.

Database Handling:
The script accepts three types of database input via the -d/--database argument:
1. Directory containing BLAST databases (auto-detects if only one database present)
2. Specific path to an existing BLAST database
3. FASTA file to create a database from (creates nucleotide database using makeblastdb)

When a FASTA file is provided for database creation:
- Checks if a database with the same name already exists in the same directory
- If database exists: uses the existing database (saves time on repeat runs)
- If database doesn't exist: creates a new nucleotide database using makeblastdb
- Database name is derived from the FASTA filename (without extension)

Input Processing:
- Single FASTA files: Processes as single sequence or splits multi-FASTA automatically
- Directories: Processes all .fasta and .fa files found in the directory
- Multi-FASTA files: Automatically splits into individual sequences for parallel processing
- Output files: Named based on sequence headers, organized in subdirectories for multi-FASTA input

Output Format:
- Tab-separated values (TSV) format with standard BLAST fields, ordered by pident (highest first)
- Default output format: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle
- One output file per input sequence
- Files named using sanitized sequence headers
- Optional summary CSV file with top hits per query in a flattened format (includes sequences with no hits)

Requirements:
- Python 3.6+
- BLAST+ suite (blastn and makeblastdb commands must be in PATH)

Usage Examples:
    # Use existing database directory (auto-detect)
    python tv_local_blast.py -i sequences.fasta -d /path/to/databases/ -o results/
    
    # Use specific existing database
    python tv_local_blast.py -i sequences.fasta -d /path/to/databases/nt -o results/
    
    # Create database from FASTA file
    python tv_local_blast.py -i queries.fasta -d reference_genome.fasta -o results/
    
    # Process directory with custom settings and CSV summary
    python tv_local_blast.py -i /fasta_dir/ -d database.fasta -o results/ -p 16 --output-csv summary.csv
    
    # Custom BLAST parameters
    python tv_local_blast.py -i input.fasta -d db.fasta -o out/ --blast-opts "-evalue 1e-10 -max_target_seqs 5"

Author: D. Parsons @NHMUK
License: N/A
Version: 2.2
"""
import os
import sys
import subprocess
import argparse
import tempfile
import shutil
import csv
from pathlib import Path
from multiprocessing import cpu_count
from concurrent.futures import ProcessPoolExecutor, as_completed
import logging
from typing import List, Tuple, Optional, Dict, Set
import re

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

class BLASTRunner:
    def __init__(self, db_path: str, output_dir: str, blast_options: str, num_processes: int = None, output_csv: str = None):
        """
        Initialize BLAST runner with configuration parameters.
        
        Args:
            db_path: Path to BLAST database directory, specific database, or FASTA file to create database from
            output_dir: Directory to save BLAST results
            blast_options: Additional BLAST options string
            num_processes: Number of parallel processes (default: CPU count)
            output_csv: Filename for summary CSV file (optional)
        """
        self.db_input = Path(db_path)
        self.output_dir = Path(output_dir)
        self.blast_options = blast_options
        self.num_processes = num_processes or cpu_count()
        self.output_csv = output_csv
        
        # Track all processed sequences for CSV summary
        self.processed_sequences: Dict[str, str] = {}  # seq_id -> original_header
        
        # Create output directory
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Find/create and validate database
        self.db_path, self.db_name = self._handle_database()
        
        # Validate setup
        self._validate_setup()
    
    def _handle_database(self) -> Tuple[Path, str]:
        """
        Handle database input - either find existing database or create from FASTA.
        
        Returns:
            Tuple of (full_db_path, db_name)
        """
        # Check if input is a FASTA file
        if self.db_input.is_file() and self.db_input.suffix.lower() in ['.fasta', '.fa', '.fas']:
            logger.info(f"FASTA file provided: {self.db_input}")
            return self._handle_fasta_input()
        else:
            # Original behavior - directory or database path
            return self._find_database()
    
    def _handle_fasta_input(self) -> Tuple[Path, str]:
        """
        Handle FASTA file input - check for existing database or create new one.
        
        Returns:
            Tuple of (full_db_path, db_name)
        """
        fasta_path = self.db_input
        db_name = fasta_path.stem  # filename without extension
        db_path = fasta_path.parent / db_name
        
        # Check if database already exists
        if self._is_valid_database(db_path):
            logger.info(f"Database already exists for {fasta_path.name}, using existing database: {db_name}")
            return db_path, db_name
        else:
            logger.info(f"Creating BLAST database from {fasta_path.name}")
            self._create_blast_database(fasta_path, db_path)
            return db_path, db_name
    
    def _create_blast_database(self, fasta_path: Path, db_path: Path) -> None:
        """
        Create BLAST database from FASTA file.
        
        Args:
            fasta_path: Path to input FASTA file
            db_path: Path for output database (without extension)
        """
        try:
            cmd = [
                'makeblastdb',
                '-in', str(fasta_path),
                '-dbtype', 'nucl',
                '-out', str(db_path),
                '-title', f"Database created from {fasta_path.name}"
            ]
            
            logger.info("Running makeblastdb...")
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            
            if result.stdout:
                logger.debug(f"makeblastdb output: {result.stdout}")
            
            logger.info(f"✓ Database created successfully: {db_path}")
            
        except subprocess.CalledProcessError as e:
            error_msg = f"Failed to create BLAST database: {e.stderr}"
            logger.error(error_msg)
            raise RuntimeError(error_msg)
        except FileNotFoundError:
            raise RuntimeError("makeblastdb command not found. Please make sure BLAST+ is installed and in your PATH")
    
    def _find_database(self) -> Tuple[Path, str]:
        """
        Find BLAST database from the given path.
        Can handle both directory (auto-detect) and specific database path.
        
        Returns:
            Tuple of (full_db_path, db_name)
        """
        if self.db_input.is_dir():
            # Directory provided - find databases automatically
            available_dbs = self._find_available_databases_in_dir(self.db_input)
            
            if not available_dbs:
                raise RuntimeError(f"No BLAST databases found in directory: {self.db_input}")
            elif len(available_dbs) == 1:
                # Only one database found - use it automatically
                db_name = available_dbs[0]
                logger.info(f"Auto-detected database: {db_name}")
                return self.db_input / db_name, db_name
            else:
                # Multiple databases found - let user choose
                db_list = '\n'.join([f"  - {db}" for db in available_dbs])
                raise RuntimeError(f"Multiple BLAST databases found in {self.db_input}:\n{db_list}\n"
                                 f"Please specify the exact database path instead of the directory.")
        else:
            # Specific file/database name provided
            if self.db_input.exists():
                # Full path to database file provided
                return self.db_input, self.db_input.name
            else:
                # Check if it's a database name in current directory or if parent dir exists
                parent_dir = self.db_input.parent
                db_name = self.db_input.name
                
                if parent_dir.exists():
                    # Check if database exists in the parent directory
                    potential_db = parent_dir / db_name
                    if self._is_valid_database(potential_db):
                        return potential_db, db_name
                
                raise FileNotFoundError(f"Database not found: {self.db_input}")
    
    def _find_available_databases_in_dir(self, db_dir: Path) -> List[str]:
        """Find available BLAST databases in a directory."""
        if not db_dir.exists():
            return []
        
        db_files = []
        # Look for nucleotide database files
        for pattern in ['*.nhr', '*.00.nhr', '*.nal']:
            db_files.extend(db_dir.glob(pattern))
        
        # Extract database names (remove extensions)
        db_names = set()
        for db_file in db_files:
            name = db_file.name
            # Remove common BLAST database extensions
            name = re.sub(r'\.(nhr|nal)$', '', name)
            name = re.sub(r'\.00\.nhr$', '', name)
            db_names.add(name)
        
        return sorted(list(db_names))
    
    def _is_valid_database(self, db_path: Path) -> bool:
        """Check if a database path points to a valid BLAST database."""
        db_extensions = ['.nhr', '.00.nhr', '.nal']
        return any((db_path.parent / f"{db_path.name}{ext}").exists() for ext in db_extensions)
    
    def _validate_setup(self):
        """Validate BLAST installation and database availability."""
        # Check if blastn is available
        try:
            subprocess.run(['blastn', '-version'], 
                         capture_output=True, check=True)
        except (subprocess.CalledProcessError, FileNotFoundError):
            raise RuntimeError("Error: blastn command not found. Please make sure BLAST+ is installed and in your PATH")
        
        # Database validation is now handled in _handle_database()
        if not self._is_valid_database(self.db_path):
            available_dbs = self._find_available_databases_in_dir(self.db_path.parent)
            raise RuntimeError(f"Error: BLAST database '{self.db_name}' not found\n"
                             f"Available databases in {self.db_path.parent}: {available_dbs}")
    
    @staticmethod
    def sanitize_header(header: str) -> str:
        """
        Sanitize FASTA header for use as filename.
        
        Args:
            header: FASTA header string
            
        Returns:
            Sanitized string safe for use as filename
        """
        # Remove '>' character and replace problematic characters
        sanitized = header.lstrip('>')
        sanitized = re.sub(r'[/|:*"<>? ]', '_', sanitized)
        # Limit length to avoid filename issues
        return sanitized[:100]
    
    def _extract_sequences_from_fasta(self, fasta_file: Path) -> Dict[str, str]:
        """
        Extract sequence IDs and headers from a FASTA file.
        
        Args:
            fasta_file: Path to FASTA file
            
        Returns:
            Dictionary mapping sanitized sequence IDs to original headers
        """
        sequences = {}
        try:
            with open(fasta_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line.startswith('>'):
                        header = line
                        seq_id = self.sanitize_header(header)
                        sequences[seq_id] = header
        except Exception as e:
            logger.warning(f"Error reading FASTA file {fasta_file}: {str(e)}")
        
        return sequences
    
    def _sort_blast_output(self, output_file: Path) -> None:
        """
        Sort BLAST output file by pident in descending order.
        
        Args:
            output_file: Path to BLAST output TSV file
        """
        try:
            # Read the file
            with open(output_file, 'r') as f:
                lines = f.readlines()
            
            if not lines:
                return
            
            # Parse lines and sort by pident (column index 2, descending)
            data_lines = []
            for line in lines:
                if line.strip():
                    parts = line.strip().split('\t')
                    if len(parts) >= 3:
                        try:
                            pident = float(parts[2])
                            data_lines.append((pident, line))
                        except (ValueError, IndexError):
                            # Keep malformed lines at the end
                            data_lines.append((0.0, line))
            
            # Sort by pident (descending)
            data_lines.sort(key=lambda x: x[0], reverse=True)
            
            # Write back sorted data
            with open(output_file, 'w') as f:
                for _, line in data_lines:
                    f.write(line)
                    
        except Exception as e:
            logger.warning(f"Failed to sort BLAST output {output_file}: {str(e)}")
    
    def _run_blast_single(self, query_file: Path, output_file: Path, header: str) -> Tuple[bool, str]:
        """
        Run BLAST on a single sequence file.
        
        Args:
            query_file: Path to query FASTA file
            output_file: Path to output file
            header: Sequence header for logging
            
        Returns:
            Tuple of (success, message)
        """
        try:
            # Construct BLAST command
            cmd = [
                'blastn',
                '-query', str(query_file),
                '-db', str(self.db_path),
                '-out', str(output_file),
                '-outfmt', '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle'
            ] + self.blast_options.split()
            
            # Run BLAST
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            
            # Sort output by pident
            self._sort_blast_output(output_file)
            
            return True, f"✓ BLAST completed for {header}"
            
        except subprocess.CalledProcessError as e:
            error_msg = f"✗ BLAST failed for {header}: {e.stderr}"
            logger.error(error_msg)
            return False, error_msg
        except Exception as e:
            error_msg = f"✗ Unexpected error for {header}: {str(e)}"
            logger.error(error_msg)
            return False, error_msg
    
    def process_single_sequence(self, args: Tuple[Path, str, Path]) -> Tuple[bool, str, str]:
        """
        Process a single sequence (wrapper for multiprocessing).
        
        Args:
            args: Tuple of (fasta_file, base_output_prefix, output_dir)
            
        Returns:
            Tuple of (success, message, header)
        """
        fasta_file, base_output_prefix, output_dir = args
        
        try:
            # Read the first line to get header
            with open(fasta_file, 'r') as f:
                header = f.readline().strip()
            
            if not header.startswith('>'):
                return False, f"✗ Error: File doesn't appear to be in FASTA format", str(fasta_file)
            
            # Create output filename
            safe_header = self.sanitize_header(header)
            if base_output_prefix:
                output_file = output_dir / f"{base_output_prefix}_{safe_header}.tsv"
            else:
                output_file = output_dir / f"{safe_header}.tsv"
            
            # Track this sequence for CSV summary
            self.processed_sequences[safe_header] = header
            
            # Check if output already exists
            if output_file.exists():
                return True, f"○ Skipping {header} - output file already exists", header
            
            # Run BLAST
            success, message = self._run_blast_single(fasta_file, output_file, header)
            return success, message, header
            
        except Exception as e:
            return False, f"✗ Error processing {fasta_file}: {str(e)}", str(fasta_file)
    
    def split_multifasta(self, fasta_file: Path) -> Tuple[List[Path], Path]:
        """
        Split multi-FASTA file into individual sequence files.
        
        Args:
            fasta_file: Path to multi-FASTA file
            
        Returns:
            Tuple of (list of paths to individual sequence files, temp directory path)
        """
        temp_dir = Path(tempfile.mkdtemp())
        split_files = []
        
        current_file = None
        current_handle = None
        seq_count = 0
        
        try:
            with open(fasta_file, 'r') as f:
                for line in f:
                    if line.startswith('>'):
                        # Close previous file if open
                        if current_handle:
                            current_handle.close()
                        
                        # Track sequence for CSV summary
                        header = line.strip()
                        safe_header = self.sanitize_header(header)
                        self.processed_sequences[safe_header] = header
                        
                        # Start new file
                        seq_count += 1
                        current_file = temp_dir / f"seq_{seq_count}.fa"
                        current_handle = open(current_file, 'w')
                        split_files.append(current_file)
                    
                    if current_handle:
                        current_handle.write(line)
            
            # Close last file
            if current_handle:
                current_handle.close()
                
        except Exception as e:
            # Clean up on error
            if current_handle:
                current_handle.close()
            shutil.rmtree(temp_dir, ignore_errors=True)
            raise e
        
        return split_files, temp_dir
    
    def process_multifasta_parallel(self, fasta_file: Path, output_dir: Path) -> Tuple[int, int]:
        """
        Process multi-FASTA file in parallel.
        
        Args:
            fasta_file: Path to multi-FASTA file
            output_dir: Directory to save results
            
        Returns:
            Tuple of (processed_count, skipped_count)
        """        
        # Split multi-FASTA file
        split_files, temp_dir = self.split_multifasta(fasta_file)
        total_sequences = len(split_files)
        
        logger.info(f"Found {total_sequences} sequences in the multi-FASTA file")
        
        try:
            # Prepare arguments for parallel processing
            process_args = [(f, None, output_dir) for f in split_files]
            
            processed_count = 0
            skipped_count = 0
            
            # Process in parallel
            with ProcessPoolExecutor(max_workers=self.num_processes) as executor:
                # Submit all jobs
                future_to_file = {
                    executor.submit(self.process_single_sequence, args): args[0] 
                    for args in process_args
                }
                
                # Collect results
                for i, future in enumerate(as_completed(future_to_file), 1):
                    file_path = future_to_file[future]
                    try:
                        success, message, header = future.result()
                        
                        logger.info(f"[{i}/{total_sequences}] {message}")
                        
                        if success:
                            if "Skipping" in message:
                                skipped_count += 1
                            else:
                                processed_count += 1
                        
                    except Exception as e:
                        logger.error(f"[{i}/{total_sequences}] ✗ Error processing {file_path}: {str(e)}")
            
            logger.info(f"Completed multi-FASTA processing: {processed_count} processed, {skipped_count} skipped")
            return processed_count, skipped_count
            
        finally:
            # Clean up temporary files
            shutil.rmtree(temp_dir, ignore_errors=True)
    
    def process_directory_parallel(self, input_dir: Path) -> None:
        """
        Process all FASTA files in a directory in parallel.
        
        Args:
            input_dir: Directory containing FASTA files
        """
        # Find all FASTA files
        fasta_patterns = ['*.fasta', '*.fa']
        fasta_files = []
        for pattern in fasta_patterns:
            fasta_files.extend(input_dir.glob(pattern))
        
        if not fasta_files:
            raise ValueError(f"No FASTA files found in {input_dir}")
        
        # Extract all sequences from all FASTA files for tracking
        for fasta_file in fasta_files:
            sequences = self._extract_sequences_from_fasta(fasta_file)
            self.processed_sequences.update(sequences)
        
        total_files = len(fasta_files)
        logger.info(f"Found {total_files} FASTA files to process")
        logger.info("Starting BLAST searches...")
        
        dir_processed_count = 0
        dir_skipped_count = 0
        
        for i, fasta_file in enumerate(fasta_files, 1):
            filename = fasta_file.name
            base_name = fasta_file.stem
            
            logger.info(f"[{i}/{total_files}] Processing file: {filename}")
            
            # Count sequences in file
            with open(fasta_file, 'r') as f:
                seq_count = sum(1 for line in f if line.startswith('>'))
            
            if seq_count == 0:
                logger.error(f"  ✗ Error: No FASTA sequences found in {fasta_file}")
                continue
            elif seq_count == 1:
                logger.info("  Single sequence FASTA file")
                
                # Check if output already exists
                with open(fasta_file, 'r') as f:
                    header = f.readline().strip()
                safe_header = self.sanitize_header(header)
                output_file = self.output_dir / f"{base_name}_{safe_header}.tsv"
                
                if output_file.exists():
                    logger.info(f"  ○ Skipping {filename} - output file already exists")
                    dir_skipped_count += 1
                else:
                    success, message, _ = self.process_single_sequence((fasta_file, base_name, self.output_dir))
                    logger.info(f"  {message}")
                    if success:
                        dir_processed_count += 1
            else:
                logger.info(f"  Multi-FASTA file with {seq_count} sequences")
                
                # Create subdirectory for this file's results
                file_output_dir = self.output_dir / base_name
                file_output_dir.mkdir(exist_ok=True)
                
                processed, skipped = self.process_multifasta_parallel(fasta_file, file_output_dir)
                dir_processed_count += processed
                dir_skipped_count += skipped
        
        logger.info(f"Directory processing completed: {dir_processed_count} sequences processed, {dir_skipped_count} sequences skipped")
    
    def _generate_summary_csv(self) -> None:
        """
        Generate summary CSV file from all TSV files in the output directory.
        Now includes sequences with no BLAST hits.
        """
        if not self.output_csv:
            return
        
        logger.info("Generating summary CSV file...")
        
        # Find all TSV files recursively in output directory
        tsv_files = list(self.output_dir.rglob("*.tsv"))
        
        if not tsv_files:
            logger.warning("No TSV files found for CSV summary generation")
            return
        
        # Extract max_target_seqs from blast_options
        max_hits = 1000  # default
        if 'max_target_seqs' in self.blast_options:
            match = re.search(r'max_target_seqs\s+(\d+)', self.blast_options)
            if match:
                max_hits = int(match.group(1))
        
        # Create column headers
        headers = ['seq_id', 'original_header']
        for i in range(1, max_hits + 1):
            headers.extend([f'hit{i}', f'hit{i}_pident', f'hit{i}_length', f'hit{i}_mismatch', f'hit{i}_evalue'])
        
        # Dictionary to store results for each sequence
        sequence_results = {}
        
        # Initialize all sequences with empty results
        for seq_id, original_header in self.processed_sequences.items():
            sequence_results[seq_id] = {
                'original_header': original_header,
                'hits': []
            }
        
        # Process TSV files to populate hits
        sequences_with_hits = 0
        sequences_without_hits = 0
        
        for tsv_file in tsv_files:
            try:
                # Read TSV file
                with open(tsv_file, 'r') as f:
                    lines = f.readlines()
                
                if not lines:
                    # Empty file - find corresponding sequence ID from filename
                    tsv_filename = tsv_file.stem
                    # Try to match with processed sequences
                    for seq_id in self.processed_sequences:
                        if seq_id in tsv_filename or tsv_filename.endswith(seq_id):
                            logger.debug(f"No hits found for sequence: {seq_id}")
                            sequences_without_hits += 1
                            break
                    continue
                
                # Parse TSV data
                current_qseqid = None
                current_hits = []
                
                for line in lines:
                    line = line.strip()
                    if not line:
                        continue
                    
                    parts = line.split('\t')
                    if len(parts) < 13:  # Ensure we have all required columns
                        continue
                    
                    qseqid = parts[0]
                    sseqid = parts[1]
                    pident = parts[2]
                    length = parts[3]
                    mismatch = parts[4]
                    evalue = parts[10]
                    
                    # If this is a new query sequence, process the previous one
                    if current_qseqid is not None and qseqid != current_qseqid:
                        # Find matching sequence ID in processed_sequences
                        matching_seq_id = None
                        for seq_id in self.processed_sequences:
                            if seq_id == current_qseqid or current_qseqid in seq_id:
                                matching_seq_id = seq_id
                                break
                        
                        if matching_seq_id and matching_seq_id in sequence_results:
                            sequence_results[matching_seq_id]['hits'] = current_hits[:max_hits]
                            if current_hits:
                                sequences_with_hits += 1
                        
                        current_hits = []
                    
                    current_qseqid = qseqid
                    current_hits.append([sseqid, pident, length, mismatch, evalue])
                
                # Process the last query sequence
                if current_qseqid is not None:
                    # Find matching sequence ID in processed_sequences
                    matching_seq_id = None
                    for seq_id in self.processed_sequences:
                        if seq_id == current_qseqid or current_qseqid in seq_id:
                            matching_seq_id = seq_id
                            break
                    
                    if matching_seq_id and matching_seq_id in sequence_results:
                        sequence_results[matching_seq_id]['hits'] = current_hits[:max_hits]
                        if current_hits:
                            sequences_with_hits += 1
                
            except Exception as e:
                logger.warning(f"Error processing TSV file {tsv_file}: {str(e)}")
                continue
        
        # Count sequences without hits
        sequences_without_hits = len([seq for seq in sequence_results.values() if not seq['hits']])
        
        # Prepare CSV data
        csv_data = []
        for seq_id, result in sequence_results.items():
            row_data = [seq_id, result['original_header']]
            
            # Add hit data (up to max_hits)
            for i in range(max_hits):
                if i < len(result['hits']):
                    row_data.extend(result['hits'][i])
                else:
                    row_data.extend(['', '', '', '', ''])  # Empty values for missing hits
            
            csv_data.append(row_data)
        
        # Write CSV file
        try:
            csv_path = self.output_dir / self.output_csv
            with open(csv_path, 'w', newline='', encoding='utf-8') as csvfile:
                writer = csv.writer(csvfile)
                writer.writerow(headers)
                writer.writerows(csv_data)
            
            logger.info(f"✓ Summary CSV file created: {csv_path}")
            logger.info(f"  - Total sequences: {len(sequence_results)}")
            logger.info(f"  - Sequences with hits: {sequences_with_hits}")
            logger.info(f"  - Sequences without hits: {sequences_without_hits}")
            logger.info(f"  - Included up to {max_hits} hits per query")
            
        except Exception as e:
            logger.error(f"Failed to create summary CSV file: {str(e)}")
    
    def process_input(self, input_path: Path) -> None:
        """
        Process input (either directory or single file).
        
        Args:
            input_path: Path to input directory or file
        """
        if not input_path.exists():
            raise FileNotFoundError(f"Input path {input_path} does not exist")
        
        # Display configuration
        logger.info("===== BLAST Configuration =====")
        logger.info(f"Database: {self.db_path}")
        logger.info(f"Input: {input_path}")
        logger.info(f"Output directory: {self.output_dir}")
        logger.info(f"BLAST options: {self.blast_options}")
        logger.info(f"Parallel processes: {self.num_processes}")
        if self.output_csv:
            logger.info(f"Summary CSV: {self.output_csv}")
        logger.info("===============================")
        
        if input_path.is_dir():
            logger.info("Input is a directory. Processing all FASTA files...")
            self.process_directory_parallel(input_path)
        else:
            logger.info("Input is a single file...")
            
            # Validate file extension (with warning)
            if input_path.suffix not in ['.fa', '.fasta']:
                logger.warning("Input file does not have .fa or .fasta extension. Continuing anyway...")
            
            filename = input_path.name
            base_name = input_path.stem
            
            # Extract sequences for tracking
            sequences = self._extract_sequences_from_fasta(input_path)
            self.processed_sequences.update(sequences)
            
            # Count sequences
            with open(input_path, 'r') as f:
                seq_count = sum(1 for line in f if line.startswith('>'))
            
            if seq_count == 0:
                raise ValueError(f"No FASTA sequences found in {input_path}")
            elif seq_count == 1:
                logger.info(f"Processing single sequence FASTA file: {filename}")
                
                with open(input_path, 'r') as f:
                    header = f.readline().strip()
                safe_header = self.sanitize_header(header)
                output_file = self.output_dir / f"{base_name}_{safe_header}.tsv"
                
                if output_file.exists():
                    logger.info(f"  ○ Skipping {filename} - output file already exists")
                else:
                    success, message, _ = self.process_single_sequence((input_path, base_name, self.output_dir))
                    logger.info(f"  {message}")
            else:
                logger.info(f"Processing multi-FASTA file: {filename} with {seq_count} sequences")
                processed, skipped = self.process_multifasta_parallel(input_path, self.output_dir)
        
        # Generate summary CSV if requested
        if self.output_csv:
            self._generate_summary_csv()
        
        logger.info(f"BLAST searches completed. Results saved to {self.output_dir}")


def main():
    """Main function to run the BLAST script."""
    parser = argparse.ArgumentParser(
        description="Run BLASTn on FASTA files in parallel",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Auto-detect database in directory (works if only one database present)
  python blast_parallel.py -i sequences.fasta -d /path/to/databases/ -o results/
  
  # Specify exact database path
  python blast_parallel.py -i sequences.fasta -d /path/to/databases/nt -o results/
  
  # Create database from FASTA file (or use existing if present)
  python blast_parallel.py -i queries.fasta -d reference.fasta -o results/
  
  # Process directory with custom options and generate summary CSV
  python blast_parallel.py -i /fasta_dir/ -d /db_path/nr -o results/ -p 16 --output-csv summary.csv
  
  # Custom BLAST parameters
  python blast_parallel.py -i input.fasta -d db/mydb -o out/ --blast-opts "-evalue 1e-10 -max_target_seqs 10"
        """
    )
    
    # Required arguments
    parser.add_argument('-i', '--input', required=True, type=Path,
                        help='Input FASTA file or directory containing FASTA files')
    parser.add_argument('-d', '--database', required=True, type=Path,
                        help='Path to BLAST database directory (auto-detects databases), specific database path, or FASTA file to create database from')
    parser.add_argument('-o', '--output', required=True, type=Path,
                        help='Output directory for BLAST results')
    
    # Optional arguments
    parser.add_argument('-p', '--processes', type=int, default=None,
                        help=f'Number of parallel processes (default: {cpu_count()})')
    parser.add_argument('--blast-opts', default='-evalue 1e-5 -max_target_seqs 100 -num_threads 1',
                        help='Additional BLAST options (default: "-evalue 1e-5 -max_target_seqs 1000 -num_threads 1")')
    parser.add_argument('--output-csv', type=str, default=None,
                        help='Generate summary CSV file with specified filename (e.g. "summary.csv")')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Enable verbose logging')
    
    args = parser.parse_args()
    
    # Set logging level
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    try:
        # Create BLAST runner
        blast_runner = BLASTRunner(
            db_path=str(args.database),
            output_dir=str(args.output),
            blast_options=args.blast_opts,
            num_processes=args.processes,
            output_csv=args.output_csv
        )
        
        # Process input
        blast_runner.process_input(args.input)
        
    except Exception as e:
        logger.error(f"Error: {str(e)}")
        sys.exit(1)


if __name__ == "__main__":
    main()
