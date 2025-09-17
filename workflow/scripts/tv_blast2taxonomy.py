#!/usr/bin/env python3
"""
tv_blast2taxonomy.py

This script validates taxonomic assignments by comparing BLAST hits against expected taxonomy 
using hierarchical taxonomic matching with the BOLD (Barcode of Life Data System) taxonomy database. 
It processes DNA barcode sequences, their BLAST results, and provides comprehensive taxonomic 
validation with detailed reporting.

The script is designed for metabarcoding studies where sequences need to be taxonomically 
validated against known reference databases, particularly useful for biodiversity assessments 
and environmental DNA studies.

PROCESS
=============
1. INPUT PARSING AND VALIDATION:
   - Parses multi-FASTA sequences for downstream filtering and output
   - Parses BLAST results from standardised CSV format with hit columns (hit1, hit2, etc.)
   - Loads expected taxonomy via Process ID mapping
   - Imports BOLD taxonomy database from TSV format for reference lookups
   
2. SEQUENCE-BY-SEQUENCE PROCESSING:
   For each query sequence:
   a) Maps seq_id to Process ID from expected taxonomy
   b) Extracts expected taxonomic lineage (order → family → genus → species)
   c) Processes all available BLAST hits and extracts BOLD taxonomies
   d) Applies configurable minimum percent identity filtering
   e) Performs hierarchical taxonomic matching against expected lineage
   f) Selects optimal hit using rank-prioritized quality assessment
   g) Generates comprehensive taxonomic observation summaries
   
3. HIERARCHICAL TAXONOMIC MATCHING STRATEGY:
   The validation employs sophisticated multi-level matching:
   - Taxonomic Scope: Restricted to order, family, genus, and species ranks only
     (class and phylum matches are excluded to prevent overly broad assignments)
   - Matching Algorithm: Uses exact string comparison at each taxonomic level
     - Each BLAST hit taxonomy is compared against all levels of expected lineage
   - Rank Prioritisation: More specific matches take precedence
     - Species matches > genus matches > family matches > order matches
     - Records the actual taxonomic rank where matches occur
   - Match Validation: Only accepts matches where observed taxonomy exactly equals
     expected taxonomy at any rank within the allowed hierarchy
     
4. BOLD DATABASE INTEGRATION:
   - Extracts BOLD Barcode Index Numbers (BINs) from BLAST sseqids
   - Retrieves complete taxonomic lineages from BOLD TSV database
   
5. QUALITY-BASED HIT SELECTION:
   Multi-criteria optimisation for selecting best taxonomic assignments:
   Primary Criterion: Taxonomic Specificity
   - Species-level matches preferred over genus, family, or order
   - Ensures most precise taxonomic assignment possible
   Secondary Criteria (for hits at same taxonomic rank):
   - Highest percent identity (sequence similarity)
   - Longest alignment length (coverage)
   - Fewest mismatches (accuracy)
   - Lowest e-value (statistical significance)
   BIN-Based Optimisation:
   - For hits from same BOLD BIN with similar identity (±5% or ≥95%), prioritises longer alignments
   - Reduces redundancy while maintaining assignment quality
   
6. OUTPUT GENERATION:
   CSV Results File:
   - Validation outcomes (YES/NO matches)
   - Matched taxonomic ranks and values
   - Hit quality statistics (identity, length, mismatches, e-value)
   - Process ID linkages and expected taxonomy information
   Enhanced obs_taxonomy Field:
   - Shows taxonomies from top 50 most relevant hits in "Taxonomy (rank)" format
   - For matches: prioritises hits that match expected lineage
   - For non-matches: shows top 50 highest-quality hits overall
   - Provides comprehensive taxonomic context for manual review
   Filtered FASTA Output:
   - Contains only sequences with successful taxonomic matches
   - Maintains original sequence identifiers and formatting
   - Ready for downstream phylogenetic or ecological analyses

COMMAND LINE INTERFACE
=======================
REQUIRED ARGUMENTS:
  --input-blast FILE         BLAST results in CSV format with standardized hit columns
                            Expected columns: seq_id, hit1, hit1_pident, hit1_length,
                            hit1_mismatch, hit1_evalue, hit1_description, hit2, etc.               
  --input-bold FILE          BOLD taxonomy database in TSV format
                            Expected columns: bin, kingdom, phylum, class, order,
                            family, subfamily, tribe, genus, species, subspecies                 
  --input-taxonomy FILE      Expected taxonomy mappings in CSV format
                            Required columns: Process ID, and taxonomic rank columns
                            (phylum, class, order, family, genus, species)                        
  --input-fasta FILE         Multi-FASTA file with query sequences
                            Headers must contain identifiers linkable to Process IDs                          
  --output-csv FILE          Validation results in CSV format
                            Contains all validation outcomes and statistics
  --output-fasta FILE        Filtered FASTA file with only matched sequences
                            Subset of input sequences with successful taxonomy validation

OPTIONAL ARGUMENTS:
  --taxval-rank RANK         Primary taxonomic rank for validation assessment
                            Choices: order, family, genus, species (default: family)
                            Note: Matching occurs at ANY rank in expected lineage;
                            this parameter primarily affects expected taxonomy extraction                     
  --min-pident FLOAT         Minimum percent identity threshold for hit inclusion
                            Range: 0.0-100.0 (default: 0.0, no filtering)
                            Filters hits before taxonomic matching and analysis                      
  --log FILE                 Optional log file path for detailed processing information
                            Default: stdout only, includes debug information and statistics

VALIDATION EXAMPLES
===================
Example 1: Successful Family-Level Match
  Expected: Cosmopterigidae (family level)
  BLAST Hit: "Pyroderces sp." → BOLD taxonomy includes Cosmopterigidae family
  Expected Lineage: Arthropoda → Insecta → Lepidoptera → Cosmopterigidae → Pyroderces
  BLAST Taxonomy: Lepidoptera → Cosmopterigidae → Pyroderces → "Pyroderces sp."
  Result: MATCH at FAMILY level
  Output: match_taxonomy=YES, matched_rank=family, matched_taxonomy=Cosmopterigidae

Example 2: Order-Level Match Due to Incomplete Reference
  Expected: Cosmopterigidae (family level)
  BLAST Hit: "Lepidoptera environmental sample" → BOLD taxonomy only to order
  Expected Lineage: Arthropoda → Insecta → Lepidoptera → Cosmopterigidae → Pyroderces  
  BLAST Taxonomy: Lepidoptera (order only)
  Result: MATCH at ORDER level
  Output: match_taxonomy=YES, matched_rank=order, matched_taxonomy=Lepidoptera

Example 3: Hit Prioritization Scenario
  Hit A: Order match (Lepidoptera), 98.5% identity, 450bp alignment
  Hit B: Family match (Cosmopterigidae), 95.0% identity, 400bp alignment
  Hit C: Family match (Cosmopterigidae), 97.2% identity, 420bp alignment
  Selected: Hit C (family rank supersedes order; highest identity among family matches)
  
Example 4: No Match Scenario
  Expected: Cosmopterigidae (family)
  Top BLAST Hits: All from Diptera, Coleoptera, etc. (different orders)
  Result: NO taxonomic matches found at any allowed rank
  Output: match_taxonomy=NO, matched_rank=null, obs_taxonomy shows top 50 overall hits

DEPENDENCIES
============
- Python 3.6+

AUTHORS
=======
Dan Parsons @NHMUK
License: MIT
Version: 2.3





PROCESSING LOGIC:
================

1. INPUT PARSING:
   - Parses input FASTA sequences for later output filtering
   - Parses BLAST results from CSV summary file
   - Parses expected taxonomy CSV with Process ID mappings
   - Loads BOLD taxonomy database from TSV file

2. SEQUENCE PROCESSING:
   For each sequence in BLAST results:
   a) Find matching Process ID from sequence identifier
   b) Extract complete expected taxonomic lineage
   c) Process all BLAST hits to extract observed taxonomies
   d) Apply minimum percent identity filtering
   e) Compare observed taxonomies against expected lineage using hierarchical matching
   f) Select best matching hit using rank-based prioritization and hit quality
   g) For obs_taxonomy output: show top 50 matching hits if matches exist, otherwise top 50 overall hits

3. HIERARCHICAL TAXONOMIC MATCHING:
   The script employs a sophisticated matching strategy that:
   - Searches for matches at ANY taxonomic rank in the expected lineage
   - Prioritizes more specific matches (species > genus > family > order)
   - Records the actual rank where matches are found
   - DOES NOT ALLOW matches at class or phylum level (stops at order)
   - Uses exact string matching only (no suffix pattern matching)

4. BOLD TAXONOMY EXTRACTION:
   - Expects single CSV summary file for --input-blast
   - Extracts BOLD BIN from hit IDs (part after pipe separator)
   - Looks up BIN in BOLD taxonomy TSV for taxonomic lineage
   - Compares each taxonomic rank against expected lineage

5. RANK-BASED HIT SELECTION:
   When multiple hits match at different taxonomic ranks:
   - Primary criterion: Most specific taxonomic rank (species beats genus beats family, etc.)
   - Secondary criterion: Highest percent identity among hits at the best rank
   - Tertiary criterion: Longest alignment length
   - Quaternary criterion: Fewest mismatches and lowest e-value

6. OUTPUT GENERATION:
   - CSV results with validation outcomes, matched ranks, and hit statistics
   - Enhanced obs_taxonomy field showing "Taxonomy (rank)" format from top 50 relevant hits
   - FASTA file containing only sequences with taxonomic matches at any rank

COMMAND LINE OPTIONS:
====================

REQUIRED ARGUMENTS:
  --input-blast FILE         Input BLAST results CSV file
  --input-bold FILE          Input BOLD taxonomy TSV file
  --input-taxonomy FILE      Expected taxonomy CSV with Process ID mappings
  --input-fasta FILE         Input multi-FASTA sequences file
  --output-csv FILE          Output validation results CSV
  --output-fasta FILE        Output FASTA with matching sequences only

OPTIONAL ARGUMENTS:
  --taxval-rank RANK         Taxonomic rank for validation
                            Choices: order, family, genus, species
                            Default: family
                            Note: Used primarily for expected taxonomy extraction; 
                            actual matching occurs at any rank in the lineage
  --min-pident FLOAT         Minimum percent identity threshold for hits
                            Default: 0.0 (no filtering)
  --log FILE                 Log file path (optional, defaults to stdout)

VALIDATION EXAMPLES:
============================
Example Scenario:
- Expected: Cosmopterigidae (family level)
- BLAST hits find: "Lepidoptera sp." → extracts "Lepidoptera" 
- Expected lineage: Arthropoda > Insecta > Lepidoptera > Cosmopterigidae > Pyroderces
- Result: MATCH at ORDER level (Lepidoptera found in expected lineage)
- Output: match_taxonomy=YES, matched_rank=order

Prioritisation Example:
- Hit A: Matches at order level (Lepidoptera), 98% identity
- Hit B: Matches at family level (Cosmopterigidae), 95% identity  
- Hit C: Matches at family level (Cosmopterigidae), 97% identity
- Selected: Hit C (family rank beats order rank, highest identity within family matches)
"""

import argparse
import csv
import logging
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

def setup_logging(log_file: str = None):
    """Set up logging configuration."""
    handlers = [logging.StreamHandler(sys.stdout)]
    
    if log_file:
        handlers.append(logging.FileHandler(log_file, mode='w'))
    
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=handlers
    )
    
    return logging.getLogger(__name__)

def parse_fasta(fasta_file: str, logger) -> Dict[str, str]:
    """Parse FASTA file and return dictionary of seq_id -> sequence."""
    sequences = {}
    current_seq = ""
    current_id = ""
    
    try:
        with open(fasta_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if current_id:
                        sequences[current_id] = current_seq
                    current_id = line[1:]  # Remove '>' character
                    current_seq = ""
                else:
                    current_seq += line
            
            # Don't forget the last sequence
            if current_id:
                sequences[current_id] = current_seq
                
        logger.info(f"Parsed {len(sequences)} sequences from {fasta_file}")
        return sequences
        
    except FileNotFoundError:
        logger.error(f"FASTA file not found: {fasta_file}")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Error parsing FASTA file: {e}")
        sys.exit(1)

def parse_blast_csv(blast_file: str, logger) -> Dict[str, List[Dict]]:
    """Parse BLAST CSV file and return dictionary of seq_id -> list of ALL hits."""
    blast_data = {}
    
    try:
        with open(blast_file, 'r') as f:
            reader = csv.DictReader(f)
            
            for row in reader:
                seq_id = row['seq_id']
                hits = []
                
                # Extract ALL hit information (not just hit1-hit20)
                i = 1
                while True:
                    hit_col = f'hit{i}'
                    pident_col = f'hit{i}_pident'
                    length_col = f'hit{i}_length'
                    mismatch_col = f'hit{i}_mismatch'
                    evalue_col = f'hit{i}_evalue'
                    description_col = f'hit{i}_description'
                    
                    # Check if this hit column exists and has data
                    if hit_col not in row or not row[hit_col] or row[hit_col].lower() == 'null':
                        break
                    
                    # Helper function to safely convert values
                    def safe_float(value):
                        if not value or value.lower() == 'null':
                            return 0.0
                        try:
                            return float(value)
                        except (ValueError, TypeError):
                            return 0.0
                    
                    def safe_int(value):
                        if not value or value.lower() == 'null':
                            return 0
                        try:
                            return int(value)
                        except (ValueError, TypeError):
                            return 0
                    
                    hits.append({
                        'hit_id': row[hit_col],
                        'pident': safe_float(row.get(pident_col, '')),
                        'length': safe_int(row.get(length_col, '')),
                        'mismatch': safe_int(row.get(mismatch_col, '')),
                        'evalue': safe_float(row.get(evalue_col, '')),
                        'description': row.get(description_col, ''),
                        'hit_num': i
                    })
                    
                    i += 1
                
                blast_data[seq_id] = hits
                
        total_hits = sum(len(hits) for hits in blast_data.values())
        logger.info(f"Parsed {total_hits} BLAST hits for {len(blast_data)} sequences from {blast_file}")
        return blast_data
        
    except FileNotFoundError:
        logger.error(f"BLAST file not found: {blast_file}")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Error parsing BLAST file: {e}")
        sys.exit(1)

def parse_bold_tsv(bold_file: str, logger) -> Dict[str, Dict]:
    """Parse BOLD TSV file and return dictionary of BIN -> taxonomy info."""
    bold_data = {}
    
    try:
        with open(bold_file, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            
            for row in reader:
                bin_id = row['bin']
                bold_data[bin_id] = {
                    'kingdom': row.get('kingdom', ''),
                    'phylum': row.get('phylum', ''),
                    'class': row.get('class', ''),
                    'order': row.get('order', ''),
                    'family': row.get('family', ''),
                    'subfamily': row.get('subfamily', ''),
                    'tribe': row.get('tribe', ''),
                    'genus': row.get('genus', ''),
                    'species': row.get('species', ''),
                    'subspecies': row.get('subspecies', '')
                }
                
        logger.info(f"Parsed taxonomy data for {len(bold_data)} BINs from {bold_file}")
        return bold_data
        
    except FileNotFoundError:
        logger.error(f"BOLD file not found: {bold_file}")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Error parsing BOLD file: {e}")
        sys.exit(1)

def parse_taxonomy_csv(taxonomy_file: str, logger) -> Dict[str, Dict]:
    """Parse taxonomy CSV file and return dictionary of Process ID -> taxonomy info."""
    taxonomy_data = {}
    
    try:
        with open(taxonomy_file, 'r') as f:
            reader = csv.DictReader(f)
            
            for row in reader:
                process_id = row['Process ID']
                taxonomy_data[process_id] = {
                    'phylum': row.get('phylum', ''),
                    'class': row.get('class', ''),
                    'order': row.get('order', ''),
                    'family': row.get('family', ''),
                    'genus': row.get('genus', ''),
                    'species': row.get('species', '')
                }
                
        logger.info(f"Parsed expected taxonomy for {len(taxonomy_data)} Process IDs from {taxonomy_file}")
        return taxonomy_data
        
    except FileNotFoundError:
        logger.error(f"Taxonomy file not found: {taxonomy_file}")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Error parsing taxonomy file: {e}")
        sys.exit(1)

def extract_bold_bin(hit_id: str) -> Optional[str]:
    """Extract BOLD BIN from hit ID (part after the pipe)."""
    if '|' in hit_id:
        return hit_id.split('|')[-1]
    return None

def find_process_id(seq_id: str, taxonomy_data: Dict[str, Dict]) -> Optional[str]:
    """Find matching Process ID by partial match of seq_id."""
    # Try exact match first
    if seq_id in taxonomy_data:
        return seq_id
    
    # Try partial matches - look for Process IDs that are substrings of seq_id
    for process_id in taxonomy_data.keys():
        if process_id in seq_id:
            return process_id
    
    return None

def get_expected_lineage(expected_taxonomy_data: Dict[str, str]) -> Dict[str, str]:
    """Extract the complete expected taxonomic lineage (restricted to allowed ranks)."""
    # Only include allowed ranks (no class or phylum)
    return {
        'order': expected_taxonomy_data.get('order', '').strip(),
        'family': expected_taxonomy_data.get('family', '').strip(),
        'genus': expected_taxonomy_data.get('genus', '').strip(),
        'species': expected_taxonomy_data.get('species', '').strip()
    }

def find_expected_taxonomy(expected_taxonomy_data: Dict, taxval_rank: str, logger, seq_id: str = None, process_id: str = None) -> Tuple[str, str]:
    """Find expected taxonomy, traversing up ranks if necessary (restricted to allowed ranks)."""
    # Define taxonomic hierarchy (from most specific to least specific, no class/phylum)
    rank_hierarchy = ['species', 'genus', 'family', 'order']
    
    # Find the starting position in the hierarchy
    try:
        start_idx = rank_hierarchy.index(taxval_rank)
    except ValueError:
        logger.error(f"Invalid taxonomic rank: {taxval_rank}")
        return "", ""
    
    # Try to find taxonomy at the specified rank and progressively higher ranks
    for current_rank in rank_hierarchy[start_idx:]:
        expected_rank_value = expected_taxonomy_data.get(current_rank, '').strip()
        
        if expected_rank_value:
            if current_rank != taxval_rank:
                logger.info(f"Expected taxonomy found at higher rank '{current_rank}' for seq_id: {seq_id}, Process ID: {process_id}")
            return expected_rank_value, current_rank
    
    # No taxonomy found at any rank
    logger.warning(f"No expected taxonomy found at any rank for seq_id: {seq_id}, Process ID: {process_id}")
    return "", ""

def find_taxonomic_match_in_lineages(observed_taxonomy: Dict[str, str], expected_lineage: Dict[str, str], logger) -> Tuple[Optional[str], Optional[str]]:
    """
    Find taxonomic matches using exact string matching only.
    Returns (matched_taxonomy, matched_rank) or (None, None)
    """
    # Define rank hierarchy (most specific to least specific, no class/phylum)
    rank_hierarchy = ['species', 'genus', 'family', 'order']
    
    # Check each rank in observed taxonomy against expected lineage
    for rank in rank_hierarchy:
        observed_at_rank = observed_taxonomy.get(rank, '').strip()
        expected_at_rank = expected_lineage.get(rank, '').strip()
        
        if observed_at_rank and expected_at_rank:
            # Exact string match only
            if observed_at_rank == expected_at_rank:
                logger.debug(f"Exact match found: '{observed_at_rank}' at rank '{rank}'")
                return observed_at_rank, rank
    
    return None, None

def get_matching_hits_and_taxonomies(hits: List[Dict], bold_data: Dict[str, Dict], expected_lineage: Dict[str, str], min_pident: float, logger) -> Tuple[List[Tuple[str, str]], List[Dict]]:
    """Get all observed taxonomies and hits that match any level in the expected lineage."""
    observed_taxonomies = []  # List of (taxonomy, rank) tuples
    matching_hits = []
    filtered_count = 0
    
    # Define allowed ranks (no class or phylum)
    allowed_ranks = ['species', 'genus', 'family', 'order']
    
    for hit in hits:
        # Apply minimum percent identity filter
        if hit['pident'] < min_pident:
            filtered_count += 1
            continue
            
        bold_bin = extract_bold_bin(hit['hit_id'])
        if not bold_bin:
            continue
            
        if bold_bin not in bold_data:
            logger.warning(f"BOLD BIN {bold_bin} not found in BOLD database")
            continue
            
        hit_taxonomy = bold_data[bold_bin]
        
        # Extract taxonomic ranks for this hit (only allowed ranks)
        hit_taxonomies = []
        
        for rank in allowed_ranks:
            hit_rank_value = hit_taxonomy.get(rank, '').strip()
            if hit_rank_value:
                hit_taxonomies.append((hit_rank_value, rank))
                observed_taxonomies.append((hit_rank_value, rank))
        
        # Check for exact matches against expected lineage
        for taxonomy, rank in hit_taxonomies:
            expected_at_rank = expected_lineage.get(rank, '')
            if expected_at_rank and taxonomy == expected_at_rank:
                # Add the taxonomic info to the hit for debugging
                hit_with_taxonomy = hit.copy()
                hit_with_taxonomy['taxonomy_match'] = taxonomy
                hit_with_taxonomy['matched_rank'] = rank
                hit_with_taxonomy['bold_bin'] = bold_bin
                matching_hits.append(hit_with_taxonomy)
                break  # Only record the most specific match for this hit
    
    if filtered_count > 0:
        logger.info(f"Filtered out {filtered_count} hits below {min_pident}% identity threshold")
    
    return observed_taxonomies, matching_hits

def sort_hits_by_quality(hits: List[Dict]) -> List[Dict]:
    """Sort hits by quality (percent identity, length, mismatch, e-value)."""
    def hit_quality_key(hit):
        return (
            hit['pident'],           # Higher is better
            hit['length'],           # Higher is better  
            -hit['mismatch'],        # Lower is better (negative for sorting)
            -hit['evalue']           # Lower is better (negative for sorting)
        )
    
    return sorted(hits, key=hit_quality_key, reverse=True)

def get_top_50_taxonomies(hits_to_use: List[Dict], bold_data: Dict[str, Dict], logger) -> str:
    """Extract taxonomies from top 50 hits and format for obs_taxonomy field."""
    # Define allowed ranks (no class or phylum)
    allowed_ranks = ['species', 'genus', 'family', 'order']
    
    obs_taxonomy_formatted = []
    unique_obs = set()  # To avoid duplicates
    hits_processed = 0
    
    for hit in hits_to_use:
        if hits_processed >= 50:
            break
            
        bold_bin = extract_bold_bin(hit['hit_id'])
        if not bold_bin:
            continue
            
        if bold_bin not in bold_data:
            continue
            
        hit_taxonomy = bold_data[bold_bin]
        hits_processed += 1
        
        # Extract taxonomic ranks for this hit (only allowed ranks)
        for rank in allowed_ranks:
            hit_rank_value = hit_taxonomy.get(rank, '').strip()
            if hit_rank_value:
                formatted_entry = f"{hit_rank_value} ({rank})"
                if formatted_entry not in unique_obs:
                    unique_obs.add(formatted_entry)
                    obs_taxonomy_formatted.append(formatted_entry)
    
    return "; ".join(obs_taxonomy_formatted)

def find_best_matching_hit(matching_hits: List[Dict]) -> Optional[Dict]:
    """Find the best matching hit using rank-based prioritization followed by hit quality."""
    if not matching_hits:
        return None
    
    # Define rank hierarchy (most specific to least specific, no class/phylum)
    rank_priority = {'species': 1, 'genus': 2, 'family': 3, 'order': 4}
    
    # Group hits by matched rank
    hits_by_rank = {}
    for hit in matching_hits:
        matched_rank = hit.get('matched_rank', 'unknown')
        if matched_rank not in hits_by_rank:
            hits_by_rank[matched_rank] = []
        hits_by_rank[matched_rank].append(hit)
    
    # Find the most specific (lowest number) rank with matches
    best_rank = None
    best_rank_priority = float('inf')
    
    for rank, hits in hits_by_rank.items():
        if rank in rank_priority:
            priority = rank_priority[rank]
            if priority < best_rank_priority:
                best_rank_priority = priority
                best_rank = rank
    
    if best_rank is None:
        # Fallback to first hit if no recognizable ranks
        return matching_hits[0]
    
    # Among hits at the best rank, select by hit quality with BIN optimization
    best_rank_hits = hits_by_rank[best_rank]
    
    # Step 1: Find hit with highest percent identity
    def initial_selection_key(hit):
        return (
            hit['pident'],           # Higher is better
            hit['length'],           # Higher is better  
            -hit['mismatch'],        # Lower is better (negative for sorting)
            -hit['evalue']           # Lower is better (negative for sorting)
        )
    
    initial_best = max(best_rank_hits, key=initial_selection_key)
    
    # Step 2: BIN-based optimization
    # Look for hits with same BIN that have similar % identity but better length
    if 'bold_bin' not in initial_best:
        return initial_best
    
    same_bin_hits = [hit for hit in best_rank_hits 
                     if hit.get('bold_bin') == initial_best['bold_bin']]
    
    if len(same_bin_hits) <= 1:
        return initial_best
    
    # Define "similar" % identity (within 5% or 95% of the top hit's identity)
    identity_threshold = max(
        initial_best['pident'] - 5.0,  # Within 5% absolute
        initial_best['pident'] * 0.95  # Within 95% relative
    )
    
    # Find hits from same BIN with similar identity
    similar_identity_hits = [hit for hit in same_bin_hits 
                           if hit['pident'] >= identity_threshold]
    
    if not similar_identity_hits:
        return initial_best
    
    # Among hits with similar identity from same BIN, prioritise by length
    def bin_optimised_key(hit):
        return (
            hit['length'],           # Higher is better (prioritize length for same BIN)
            hit['pident'],           # Higher is better
            -hit['mismatch'],        # Lower is better
            -hit['evalue']           # Lower is better
        )
    
    bin_optimized_best = max(similar_identity_hits, key=bin_optimised_key)
    
    return bin_optimized_best

def write_output_csv(results: List[Dict], output_file: str, logger):
    """Write results to output CSV file."""
    try:
        with open(output_file, 'w', newline='') as f:
            fieldnames = ['seq_id', 'Process_ID', 'taxval_rank', 'expected_taxonomy', 'expected_taxonomy_rank', 
                         'obs_taxonomy', 'match_taxonomy', 'matched_rank', 'top_matching_hit', 
                         'pident', 'length', 'mismatch', 'evalue']
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            
            for result in results:
                writer.writerow(result)
                
        logger.info(f"Results written to {output_file}")
        
    except Exception as e:
        logger.error(f"Error writing output CSV: {e}")
        sys.exit(1)

def write_output_fasta(matching_sequences: Dict[str, str], output_file: str, logger):
    """Write matching sequences to output FASTA file."""
    try:
        with open(output_file, 'w') as f:
            for seq_id, sequence in matching_sequences.items():
                f.write(f">{seq_id}\n{sequence}\n")
                
        logger.info(f"Matching sequences written to {output_file} ({len(matching_sequences)} sequences)")
        
    except Exception as e:
        logger.error(f"Error writing output FASTA: {e}")
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(
        description="Validate taxonomic assignments from BLAST results using BOLD database"
    )
    parser.add_argument('--input-blast', required=True, 
                       help='Input BLAST results CSV file')
    parser.add_argument('--input-bold', required=True,
                       help='Input BOLD taxonomy TSV file')
    parser.add_argument('--input-taxonomy', required=True,
                       help='Input expected taxonomy CSV file')
    parser.add_argument('--input-fasta', required=True,
                       help='Input multi-FASTA file with sequences')
    parser.add_argument('--output-csv', required=True,
                       help='Output results CSV file')
    parser.add_argument('--output-fasta', required=True,
                       help='Output FASTA file with matching sequences')
    parser.add_argument('--taxval-rank', default='family',
                       choices=['order', 'family', 'genus', 'species'],
                       help='Taxonomic rank to validate at (default: family)')
    parser.add_argument('--min-pident', type=float, default=0.0,
                       help='Minimum percent identity threshold for hits to be considered (default: 0.0, no filtering)')
    parser.add_argument('--log', 
                       help='Log file path (optional, defaults to stdout only)')
    
    args = parser.parse_args()
    
    # Set up logging
    logger = setup_logging(args.log)
    
    logger.info("Starting BLAST taxonomy validation using BOLD database")
    logger.info(f"Validation rank: {args.taxval_rank}")
    logger.info(f"Minimum percent identity: {args.min_pident}%")
    logger.info("Rank restrictions: Only matches at order, family, genus, and species levels allowed")
    logger.info("obs_taxonomy output: Top 50 matching hits (if matches exist) or top 50 overall hits")
    logger.info(f"BOLD database: {args.input_bold}")
    logger.info(f"BLAST input (CSV file): {args.input_blast}")
    
    if args.log:
        logger.info(f"Log file: {args.log}")
    
    # Parse input files
    sequences = parse_fasta(args.input_fasta, logger)
    blast_data = parse_blast_csv(args.input_blast, logger)
    bold_data = parse_bold_tsv(args.input_bold, logger)
    taxonomy_data = parse_taxonomy_csv(args.input_taxonomy, logger)
    
    results = []
    matching_sequences = {}
    
    # Process each sequence
    for seq_id, hits in blast_data.items():
        logger.info(f"Processing sequence: {seq_id} ({len(hits)} hits)")
        
        # Find corresponding Process ID
        process_id = find_process_id(seq_id, taxonomy_data)
        if not process_id:
            logger.error(f"No matching Process ID found for sequence: {seq_id}")
            sys.exit(1)
        
        # Get complete expected taxonomic lineage (restricted to allowed ranks)
        expected_lineage = get_expected_lineage(taxonomy_data[process_id])
        logger.info(f"Expected lineage (restricted): {expected_lineage}")
        
        # Find expected taxonomy at specified rank (with traversal if needed)
        expected_taxonomy, expected_taxonomy_rank = find_expected_taxonomy(
            taxonomy_data[process_id], args.taxval_rank, logger, seq_id, process_id)
        
        # Get observed taxonomies and matching hits
        obs_taxonomies, matching_hits = get_matching_hits_and_taxonomies(
            hits, bold_data, expected_lineage, args.min_pident, logger)
        
        # Determine obs_taxonomy string from top 50 relevant hits
        if matching_hits:
            # Use top 50 matching hits, sorted by quality (same as overall hits)
            sorted_matching_hits = sort_hits_by_quality(matching_hits)
            obs_taxonomy_str = get_top_50_taxonomies(sorted_matching_hits, bold_data, logger)
            logger.info(f"Using top 50 from {len(matching_hits)} matching hits for obs_taxonomy")
        else:
            # Use top 50 overall hits, sorted by quality
            # Filter hits by min_pident first
            filtered_hits = [hit for hit in hits if hit['pident'] >= args.min_pident]
            sorted_overall_hits = sort_hits_by_quality(filtered_hits)
            obs_taxonomy_str = get_top_50_taxonomies(sorted_overall_hits, bold_data, logger)
            logger.info(f"No matches found. Using top 50 from {len(filtered_hits)} overall hits for obs_taxonomy")
        
        # Log matching hits for debugging
        if matching_hits:
            logger.info(f"Found {len(matching_hits)} hits matching expected lineage:")
            for hit in matching_hits[:5]:  # Show first 5 matching hits
                matched_rank = hit.get('matched_rank', 'unknown')
                matched_taxonomy = hit.get('taxonomy_match', 'unknown')
                logger.info(f"  - {hit['hit_id']}: {matched_taxonomy} at {matched_rank} rank, {hit['pident']}% identity, {hit['length']}bp, BIN: {hit.get('bold_bin', 'unknown')}")
        
        # Determine if there's a match and find the best matching hit
        if matching_hits:
            match_taxonomy = "YES"
            best_hit = find_best_matching_hit(matching_hits)
            matched_rank = best_hit.get('matched_rank', 'unknown')
            matched_taxonomy = best_hit.get('taxonomy_match', 'unknown')
            logger.info(f"Selected best matching hit: {best_hit['hit_id']} (taxonomy: {matched_taxonomy} at {matched_rank} rank, pident: {best_hit['pident']}%, length: {best_hit['length']}bp)")
        else:
            match_taxonomy = "NO"
            matched_rank = "null"
            # For non-matching cases, still show the top overall hit
            filtered_hits = [hit for hit in hits if hit['pident'] >= args.min_pident]
            best_hit = filtered_hits[0] if filtered_hits else None
            if best_hit:
                logger.info(f"No taxonomic matches found. Showing top overall hit: {best_hit['hit_id']}")
        
        # Log the processing results for this sample
        logger.info(f"Sample processing results:")
        logger.info(f"  taxval_rank: {args.taxval_rank}")
        logger.info(f"  expected_taxonomy: {expected_taxonomy}")
        logger.info(f"  expected_taxonomy_rank: {expected_taxonomy_rank}")
        logger.info(f"  obs_taxonomy: {obs_taxonomy_str[:200]}..." if len(obs_taxonomy_str) > 200 else f"  obs_taxonomy: {obs_taxonomy_str}")
        logger.info(f"  match_taxonomy: {match_taxonomy}")
        logger.info(f"  matched_rank: {matched_rank}")
        
        # Create result dictionary
        result = {
            'seq_id': seq_id,
            'Process_ID': process_id,
            'taxval_rank': args.taxval_rank,
            'expected_taxonomy': expected_taxonomy,
            'expected_taxonomy_rank': expected_taxonomy_rank,
            'obs_taxonomy': obs_taxonomy_str,
            'match_taxonomy': match_taxonomy,
            'matched_rank': matched_rank,
            'top_matching_hit': best_hit['hit_id'] if best_hit else 'No hits',
            'pident': best_hit['pident'] if best_hit else '',
            'length': best_hit['length'] if best_hit else '',
            'mismatch': best_hit['mismatch'] if best_hit else '',
            'evalue': best_hit['evalue'] if best_hit else ''
        }
        
        # Add to matching sequences for FASTA output if there was a match
        if match_taxonomy == "YES":
            if seq_id in sequences:
                matching_sequences[seq_id] = sequences[seq_id]
            else:
                logger.warning(f"Sequence {seq_id} not found in input FASTA")
        
        results.append(result)
    
    # Write output files
    write_output_csv(results, args.output_csv, logger)
    write_output_fasta(matching_sequences, args.output_fasta, logger)
    
    # Summary statistics
    total_sequences = len(results)
    matched_sequences = len(matching_sequences)
    logger.info(f"Summary: {matched_sequences}/{total_sequences} sequences had matches")
    logger.info("BLAST taxonomy validation completed successfully")

if __name__ == "__main__":
    main()
