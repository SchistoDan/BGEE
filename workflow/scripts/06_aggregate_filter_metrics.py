#!/usr/bin/env python3
"""
Aggregate Filter Metrics
Combines metrics from all filtering steps and generates final statistics
Updated to source data directly from specific columns in each input CSV
"""

import os
import sys
import csv
import argparse
import pandas as pd
from datetime import datetime
from collections import defaultdict
from typing import Dict, List

def log_message(message: str, log_file=None, stdout=False):
    """Log message to file and optionally stdout"""
    timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    formatted_msg = f"[{timestamp}] {message}"
    
    if log_file:
        log_file.write(formatted_msg + '\n')
        log_file.flush()
    
    if stdout:
        print(formatted_msg, flush=True)

def normalise_base_name(name: str) -> tuple:
    """Normalise base names by removing common suffixes to get core base name
    Returns: (normalised_name, removed_suffix)"""
    # Remove suffixes in order of priority (longest first to avoid partial matches)
    suffixes_to_remove = [
        '_fcleaner_merge',
        '_fcleaner',
        '_human_filtered',
        '_at_filtered', 
        '_outlier_filtered',
        '_reference_filtered'
    ]
    
    normalised = name.strip()
    removed_suffix = ''
    
    for suffix in suffixes_to_remove:
        if normalised.endswith(suffix):
            normalised = normalised[:-len(suffix)]
            removed_suffix = suffix
            break  # Only remove one suffix to avoid over-trimming
    
    return normalised, removed_suffix

def parse_human_metrics(csv_file: str) -> Dict[str, Dict]:
    """Parse human COX1 filter metrics CSV"""
    metrics = {}
    try:
        with open(csv_file, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                if row.get('sequence_id') == 'FILE_SUMMARY':
                    base_name = row.get('base_name', '')
                    if base_name:
                        normalised_name, removed_suffix = normalise_base_name(base_name)
                        metrics[normalised_name] = {
                            'input_reads': int(row.get('input_count', 0) or 0),
                            'removed_human': int(row.get('removed_count', 0) or 0),
                            'original_name': base_name,
                            'removed_suffix': removed_suffix
                        }
    except Exception as e:
        print(f"Warning: Could not read human metrics {csv_file}: {str(e)}")
    return metrics

def parse_at_metrics(csv_file: str) -> Dict[str, Dict]:
    """Parse AT content filter metrics CSV"""
    metrics = {}
    try:
        with open(csv_file, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                sample_name = row.get('sample_name', '')
                if sample_name:
                    normalised_name, removed_suffix = normalise_base_name(sample_name)
                    metrics[normalised_name] = {
                        'removed_at_distance': int(row.get('removed_sequences', 0) or 0),
                        'original_name': sample_name,
                        'removed_suffix': removed_suffix
                    }
    except Exception as e:
        print(f"Warning: Could not read AT metrics {csv_file}: {str(e)}")
    return metrics

def parse_outlier_metrics(csv_file: str) -> Dict[str, Dict]:
    """Parse statistical outlier filter metrics CSV"""
    metrics = {}
    try:
        with open(csv_file, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                base_name = row.get('base_name', '')
                if base_name:
                    normalised_name, removed_suffix = normalise_base_name(base_name)
                    metrics[normalised_name] = {
                        'removed_outliers': int(row.get('removed_count', 0) or 0),
                        'original_name': base_name,
                        'removed_suffix': removed_suffix
                    }
    except Exception as e:
        print(f"Warning: Could not read outlier metrics {csv_file}: {str(e)}")
    return metrics

def parse_reference_metrics(csv_file: str) -> Dict[str, Dict]:
    """Parse reference filter metrics CSV"""
    metrics = {}
    try:
        with open(csv_file, 'r') as f:
            reader = csv.DictReader(f)
            
            for row in reader:
                # Only process FILE_SUMMARY records
                if row.get('sequence_id') == 'FILE_SUMMARY':
                    # Extract sample identifier from file_path instead of using base_name
                    file_path = row.get('file_path', '')
                    if file_path:
                        # Extract filename from path
                        filename = os.path.basename(file_path)
                        
                        # Remove the _outlier_filtered.fasta suffix to get the sample identifier
                        if filename.endswith('_outlier_filtered.fasta'):
                            sample_identifier = filename[:-len('_outlier_filtered.fasta')]
                        else:
                            # Fallback: just remove .fasta extension
                            sample_identifier = filename.replace('.fasta', '')
                        
                        # Use the full sample_identifier as the key
                        metrics[sample_identifier] = {
                            'removed_reference': int(row.get('removed_count', 0) or 0),
                            'original_name': sample_identifier,
                            'removed_suffix': ''
                        }
                        
    except Exception as e:
        print(f"Warning: Could not read reference metrics {csv_file}: {str(e)}")
    return metrics

def parse_consensus_metrics(csv_file: str) -> Dict[str, Dict]:
    """Parse consensus generation metrics CSV"""
    metrics = {}
    try:
        with open(csv_file, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                sample_name = row.get('sample_name', '')
                if sample_name:
                    # Keep the full sample_name with _fcleaner suffix for output
                    # But also create normalised version for matching
                    normalised_name, removed_suffix = normalise_base_name(sample_name)
                    metrics[normalised_name] = {
                        'sample_name': sample_name,  # Keep original with _fcleaner suffix
                        'cleaned_reads': int(row.get('cleaned_reads', 0) or 0),
                        'final_ambig_bases': int(row.get('ambiguous_bases', 0) or 0),
                        'cov_percent': float(row.get('cov_percent', 0) or 0),
                        'cov_avg': float(row.get('cov_avg', 0) or 0),
                        'cov_med': float(row.get('cov_med', 0) or 0),
                        'cov_max': int(row.get('cov_max', 0) or 0),
                        'cov_min': int(row.get('cov_min', 0) or 0),
                        'original_name': sample_name,
                        'removed_suffix': removed_suffix
                    }
    except Exception as e:
        print(f"Warning: Could not read consensus metrics {csv_file}: {str(e)}")
    return metrics

def aggregate_metrics(human_metrics: Dict, at_metrics: Dict, outlier_metrics: Dict,
                     reference_metrics: Dict, consensus_metrics: Dict) -> Dict[str, Dict]:
    """Aggregate metrics from all filtering steps"""
    
    # Get all unique normalised base names
    all_base_names = set()
    all_base_names.update(human_metrics.keys())
    all_base_names.update(at_metrics.keys())
    all_base_names.update(outlier_metrics.keys())
    all_base_names.update(reference_metrics.keys())
    all_base_names.update(consensus_metrics.keys())
    
    combined_metrics = {}
    
    for normalised_name in all_base_names:
        # Get data from each step (use 0 for missing data)
        human_data = human_metrics.get(normalised_name, {})
        at_data = at_metrics.get(normalised_name, {})
        outlier_data = outlier_metrics.get(normalised_name, {})
        reference_data = reference_metrics.get(normalised_name, {})
        consensus_data = consensus_metrics.get(normalised_name, {})
        
        # Determine the appropriate sample_name for output
        # Priority: 1) Use consensus sample_name if available, 2) Reconstruct with appropriate suffix
        if 'sample_name' in consensus_data:
            sample_name = consensus_data['sample_name']
        else:
            # Find any _fcleaner suffix from any of the sources
            fcleaner_suffix = ''
            for data in [consensus_data, human_data, at_data, outlier_data, reference_data]:
                if data.get('removed_suffix') in ['_fcleaner', '_fcleaner_merge']:
                    fcleaner_suffix = data['removed_suffix']
                    break
            
            # If no _fcleaner suffix found, default to _fcleaner
            if not fcleaner_suffix:
                fcleaner_suffix = '_fcleaner'
                
            sample_name = f"{normalised_name}{fcleaner_suffix}"
        
        combined_metrics[normalised_name] = {
            'sample_name': sample_name,
            'input_reads': human_data.get('input_reads', 0),
            'removed_human': human_data.get('removed_human', 0),
            'removed_at_distance': at_data.get('removed_at_distance', 0),
            'removed_outliers': outlier_data.get('removed_outliers', 0),
            'removed_reference': reference_data.get('removed_reference', 0),
            'cleaned_reads': consensus_data.get('cleaned_reads', 0),
            'final_ambig_bases': consensus_data.get('final_ambig_bases', 0),
            'cov_percent': consensus_data.get('cov_percent', 0.0),
            'cov_avg': consensus_data.get('cov_avg', 0.0),
            'cov_med': consensus_data.get('cov_med', 0.0),
            'cov_max': consensus_data.get('cov_max', 0),
            'cov_min': consensus_data.get('cov_min', 0)
        }
    
    return combined_metrics

def main():
    parser = argparse.ArgumentParser(description='Aggregate metrics from all filtering steps')
    parser.add_argument('--human-metrics', help='Human COX1 filter metrics CSV (optional)')
    parser.add_argument('--at-metrics', help='AT content filter metrics CSV (optional)')
    parser.add_argument('--outlier-metrics', help='Statistical outlier filter metrics CSV (optional)')
    parser.add_argument('--reference-metrics', help='Reference filter metrics CSV (optional)')
    parser.add_argument('--consensus-metrics', help='Consensus generation metrics CSV (optional)')
    parser.add_argument('--output-dir', required=True, help='Output directory')
    parser.add_argument('--combined-statistics', required=True, help='Output combined statistics CSV')
    parser.add_argument('--threads', type=int, default=1, help='Number of threads (for compatibility)')
    
    args = parser.parse_args()
    
    print("Aggregating filter metrics from all processing steps...")
    
    # Parse all metrics files with new specific parsers (only if provided and exist)
    print("Parsing metrics files...")
    
    human_metrics = {}
    if args.human_metrics and os.path.exists(args.human_metrics):
        human_metrics = parse_human_metrics(args.human_metrics)
        print(f"  Human filter: {len(human_metrics)} samples")
    else:
        print(f"  Human filter: SKIPPED (file not provided or doesn't exist)")
    
    at_metrics = {}
    if args.at_metrics and os.path.exists(args.at_metrics):
        at_metrics = parse_at_metrics(args.at_metrics)
        print(f"  AT content filter: {len(at_metrics)} samples")
    else:
        print(f"  AT content filter: SKIPPED (file not provided or doesn't exist)")
    
    outlier_metrics = {}
    if args.outlier_metrics and os.path.exists(args.outlier_metrics):
        outlier_metrics = parse_outlier_metrics(args.outlier_metrics)
        print(f"  Outlier filter: {len(outlier_metrics)} samples")
    else:
        print(f"  Outlier filter: SKIPPED (file not provided or doesn't exist)")
    
    reference_metrics = {}
    if args.reference_metrics and os.path.exists(args.reference_metrics):
        reference_metrics = parse_reference_metrics(args.reference_metrics)
        print(f"  Reference filter: {len(reference_metrics)} samples")
    else:
        print(f"  Reference filter: SKIPPED (file not provided or doesn't exist)")
    
    consensus_metrics = {}
    if args.consensus_metrics and os.path.exists(args.consensus_metrics):
        consensus_metrics = parse_consensus_metrics(args.consensus_metrics)
        print(f"  Consensus generation: {len(consensus_metrics)} samples")
    else:
        print(f"  Consensus generation: SKIPPED (file not provided or doesn't exist)")
    
    # Check if we have any data at all
    total_files_found = len([x for x in [human_metrics, at_metrics, outlier_metrics, reference_metrics, consensus_metrics] if x])
    if total_files_found == 0:
        print("\nError: No valid metrics files found. At least one metrics file must be provided.")
        sys.exit(1)
    
    print(f"\nFound data from {total_files_found} metrics files.")
    
    # Aggregate all metrics
    combined_metrics = aggregate_metrics(
        human_metrics, at_metrics, outlier_metrics, reference_metrics, consensus_metrics
    )
    
    print(f"Combined metrics for {len(combined_metrics)} samples")
    
    # Write combined statistics CSV with specified headers
    fieldnames = [
        'sample_name', 'input_reads', 'removed_human', 'removed_at_distance',
        'removed_outliers', 'removed_reference', 'cleaned_reads', 'final_ambig_bases',
        'cov_percent', 'cov_avg', 'cov_med', 'cov_max', 'cov_min'
    ]
    
    with open(args.combined_statistics, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        
        for normalised_name, metrics in combined_metrics.items():
            # Ensure all fields are present
            row = {field: metrics.get(field, 0) for field in fieldnames}
            writer.writerow(row)
    
    # Generate summary statistics
    total_samples = len(combined_metrics)
    successful_samples = len([m for m in combined_metrics.values() if m['cleaned_reads'] > 0])
    total_input_reads = sum(m['input_reads'] for m in combined_metrics.values())
    total_cleaned_reads = sum(m['cleaned_reads'] for m in combined_metrics.values())
    total_removed_human = sum(m['removed_human'] for m in combined_metrics.values())
    total_removed_at = sum(m['removed_at_distance'] for m in combined_metrics.values())
    total_removed_outliers = sum(m['removed_outliers'] for m in combined_metrics.values())
    total_removed_reference = sum(m['removed_reference'] for m in combined_metrics.values())
    total_ambiguous_bases = sum(m['final_ambig_bases'] for m in combined_metrics.values())
    
    print(f"\nFilter Aggregation Summary:")
    print(f"Total samples processed: {total_samples}")
    print(f"Samples with final sequences: {successful_samples}")
    print(f"Total input sequences: {total_input_reads}")
    print(f"Total cleaned sequences: {total_cleaned_reads}")
    print(f"")
    print(f"Sequences removed by filter:")
    print(f"  Human COX1 similarity: {total_removed_human}")
    print(f"  AT content difference: {total_removed_at}")
    print(f"  Statistical outliers: {total_removed_outliers}")
    print(f"  Reference outliers: {total_removed_reference}")
    print(f"")
    print(f"Total ambiguous bases in consensus: {total_ambiguous_bases}")
    print(f"")
    print(f"Combined statistics written to: {args.combined_statistics}")

if __name__ == "__main__":
    main()
