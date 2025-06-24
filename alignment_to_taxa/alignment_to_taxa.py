#!/usr/bin/env python3
"""
Taxonomy Assignment Tool

This script processes BAM alignments and taxonomy files to assign taxonomy to sequences
and create abundance tables with hash-based identifiers.

Usage:
    ./taxonomy_processor.py --bam INPUT.bam --abundance ABUNDANCE.tsv --taxonomy TAXONOMY.tsv
                           [--output-dir OUTPUT_DIR] [--score-threshold 0.97] [--lca-threshold 0.8]
"""

import pysam
import pandas as pd
from collections import defaultdict, Counter
import os
import hashlib
import argparse
import sys
from datetime import datetime

def parse_arguments():
    parser = argparse.ArgumentParser(description='Assign taxonomy to sequences and create abundance tables')
    
    # Required arguments
    parser.add_argument('--bam', required=True, help='Input BAM file with alignments')
    parser.add_argument('--abundance', required=True, help='Tab-separated abundance file')
    parser.add_argument('--taxonomy', required=True, help='Tab-separated taxonomy file')
    
    # Optional arguments
    parser.add_argument('--output-dir', default='.', help='Output directory (default: current directory)')
    parser.add_argument('--score-threshold', type=float, default=0.97, 
                        help='Identity threshold for keeping multiple matches (default: 0.97)')
    parser.add_argument('--lca-threshold', type=float, default=0.8, 
                        help='Agreement required at each rank to retain it (default: 0.8)')
    parser.add_argument('--prefix', default='', 
                        help='Prefix for output files (default: none)')
    parser.add_argument('--verbose', action='store_true', help='Print detailed progress')
    
    return parser.parse_args()

def log_message(message, verbose=True):
    """Print log message with timestamp"""
    timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    print(f"[{timestamp}] {message}", flush=True)

def hash_taxonomy(taxonomy):
    """Generate hash ID for a taxonomy string"""
    return hashlib.md5(taxonomy.encode()).hexdigest()[:10]

def main():
    # Parse command line arguments
    args = parse_arguments()
    
    # Define taxonomy ranks
    RANKS = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
    
    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Define output file paths
    prefix = args.prefix + "_" if args.prefix else ""
    abundance_output = os.path.join(args.output_dir, f"{prefix}abundance_table.tsv")
    taxonomy_output = os.path.join(args.output_dir, f"{prefix}taxonomy_table.tsv")
    mapping_output = os.path.join(args.output_dir, f"{prefix}sequence_hash_mapping.tsv")
    log_output = os.path.join(args.output_dir, f"{prefix}taxonomy_assignment_log.txt")
    
    # Initialize log file
    with open(log_output, 'w') as log_file:
        log_file.write(f"Taxonomy Assignment Tool - Run on {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        log_file.write(f"BAM file: {args.bam}\n")
        log_file.write(f"Abundance file: {args.abundance}\n")
        log_file.write(f"Taxonomy file: {args.taxonomy}\n")
        log_file.write(f"Score threshold: {args.score_threshold}\n")
        log_file.write(f"LCA threshold: {args.lca_threshold}\n\n")
    
    # Load taxonomy into a dictionary
    log_message(f"Loading taxonomy from {args.taxonomy}", args.verbose)
    taxonomy_dict = {}
    with open(args.taxonomy) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) == 2:
                taxonomy_dict[parts[0]] = parts[1]
    log_message(f"Loaded {len(taxonomy_dict)} taxonomy entries", args.verbose)
    
    # Load abundance data
    log_message(f"Loading abundance data from {args.abundance}", args.verbose)
    abundance = defaultdict(lambda: defaultdict(int))
    sample_header = []
    with open(args.abundance) as f:
        header = f.readline().strip().split("\t")[1:]
        sample_header = header
        for line in f:
            parts = line.strip().split("\t")
            seq_id = parts[0]
            counts = list(map(int, parts[1:]))
            for sample, count in zip(header, counts):
                abundance[seq_id][sample] = count
    log_message(f"Loaded abundance data for {len(abundance)} sequences", args.verbose)
    
    # Process alignments
    log_message(f"Processing alignments from {args.bam}", args.verbose)
    bamfile = pysam.AlignmentFile(args.bam, "rb")
    query_hits = defaultdict(list)
    all_refs = set()
    
    for read in bamfile:
        if not read.is_unmapped:
            query = read.query_name
            ref = bamfile.get_reference_name(read.reference_id)
            all_refs.add(ref)
            score = read.get_tag("AS") if read.has_tag("AS") else read.mapping_quality / 100.0
            query_hits[query].append((ref, score))
    bamfile.close()
    log_message(f"Processed {len(query_hits)} queries with alignments", args.verbose)
    log_message(f"Found {len(all_refs)} unique reference sequences", args.verbose)
    
    # Check if refs match taxonomy keys directly
    if args.verbose:
        test_refs = list(all_refs)[:min(100, len(all_refs))]
        ref_in_tax = sum(1 for ref in test_refs if ref in taxonomy_dict)
        log_message(f"Direct ref matches in taxonomy: {ref_in_tax}/{len(test_refs)}", args.verbose)
        
        if ref_in_tax < len(test_refs) * 0.5:
            log_message("Warning: Many references don't match taxonomy keys directly", args.verbose)
    
    def get_lca(taxa_list):
        """Compute the lowest common ancestor for a list of taxonomies"""
        # Parse taxonomies into lists of ranks
        parsed_taxa = []
        for tax in taxa_list:
            ranks = tax.split(";")
            while len(ranks) < len(RANKS):  # Ensure we have entries for all ranks
                ranks.append("")
            parsed_taxa.append(ranks)
        
        # Compute LCA rank by rank
        lca_ranks = []
        for i in range(len(RANKS)):
            # Get all values at this rank
            rank_values = [taxa[i] for taxa in parsed_taxa if taxa[i]]
            
            # Skip if no values at this rank
            if not rank_values:
                lca_ranks.append("")
                continue
                
            # Count occurrences of each value
            counts = Counter(rank_values)
            most_common, count = counts.most_common(1)[0]
            
            # Check if we meet the threshold
            if count / len(rank_values) >= args.lca_threshold:
                lca_ranks.append(most_common)
            else:
                # Add empty ranks for the rest
                lca_ranks.extend([""] * (len(RANKS) - i))
                break
        
        return lca_ranks
    
    # Compute taxonomy assignments
    log_message("Computing taxonomy assignments", args.verbose)
    query_taxonomy = {}
    n_single = 0
    n_lca = 0
    n_unmatched = 0
    total_refs_checked = 0
    tax_lookup_failures = 0
    
    # Process in smaller batches to provide progress feedback
    queries = list(query_hits.keys())
    total_queries = len(queries)
    batch_size = min(100000, max(1, total_queries // 10))  # 10 progress updates or 100k batch size
    n_batches = (total_queries + batch_size - 1) // batch_size
    
    for batch_idx in range(n_batches):
        start_idx = batch_idx * batch_size
        end_idx = min(start_idx + batch_size, total_queries)
        batch_queries = queries[start_idx:end_idx]
        
        for query in batch_queries:
            hits = query_hits[query]
            max_score = max(score for _, score in hits)
            filtered = [(ref, score) for ref, score in hits if score >= args.score_threshold * max_score]
            
            # Get taxonomies for all references that passed the score filter
            taxonomies = []
            for ref, _ in filtered:
                total_refs_checked += 1
                
                if ref in taxonomy_dict:
                    tax = taxonomy_dict[ref]
                    if tax:  # Make sure taxonomy is not empty
                        taxonomies.append(tax)
                else:
                    tax_lookup_failures += 1
            
            if len(taxonomies) == 0:
                n_unmatched += 1
                continue
            elif len(taxonomies) == 1:
                query_taxonomy[query] = taxonomies[0].split(";")
                # Pad with empty strings if needed
                while len(query_taxonomy[query]) < len(RANKS):
                    query_taxonomy[query].append("")
                n_single += 1
            else:
                lca_ranks = get_lca(taxonomies)
                if any(lca_ranks):  # Check if we have at least one non-empty rank
                    query_taxonomy[query] = lca_ranks
                    n_lca += 1
                else:
                    n_unmatched += 1
        
        if args.verbose:
            log_message(f"Processed {end_idx}/{total_queries} queries ({end_idx/total_queries:.1%})", args.verbose)
    
    # Print and log taxonomy lookup stats
    lookup_failure_rate = tax_lookup_failures/total_refs_checked if total_refs_checked > 0 else 0
    log_message(f"Total reference lookups: {total_refs_checked}", args.verbose)
    log_message(f"Taxonomy lookup failures: {tax_lookup_failures} ({lookup_failure_rate:.1%})", args.verbose)
    
    with open(log_output, 'a') as log_file:
        log_file.write(f"Total reference lookups: {total_refs_checked}\n")
        log_file.write(f"Taxonomy lookup failures: {tax_lookup_failures} ({lookup_failure_rate:.1%})\n")
    
    # Create a mapping from taxonomy string to hash ID
    log_message("Creating hash IDs for taxonomies", args.verbose)
    tax_to_hash = {}
    for query, tax_ranks in query_taxonomy.items():
        tax_string = ";".join(tax_ranks)
        if tax_string not in tax_to_hash:
            tax_to_hash[tax_string] = f"{hash_taxonomy(tax_string)}"
    
    # Map queries to hash IDs
    query_to_hash = {query: tax_to_hash[";".join(tax_ranks)] for query, tax_ranks in query_taxonomy.items()}
    
    # Collect unique taxonomies and their hash IDs
    unique_taxonomies = {}
    for tax_string, hash_id in tax_to_hash.items():
        tax_ranks = tax_string.split(";")
        while len(tax_ranks) < len(RANKS):
            tax_ranks.append("")
        unique_taxonomies[hash_id] = tax_ranks
    
    # Merge counts by hash ID
    log_message("Merging counts by taxonomy", args.verbose)
    merged_counts = defaultdict(lambda: defaultdict(int))
    for seq_id, hash_id in query_to_hash.items():
        for sample, count in abundance[seq_id].items():
            merged_counts[hash_id][sample] += count
    
    # Write hash-based abundance table
    log_message(f"Writing abundance table to {abundance_output}", args.verbose)
    with open(abundance_output, "w") as f:
        f.write("TaxID\t" + "\t".join(sample_header) + "\n")
        for hash_id in sorted(merged_counts):
            counts = [str(merged_counts[hash_id][s]) for s in sample_header]
            f.write(hash_id + "\t" + "\t".join(counts) + "\n")
    
    # Write tabular taxonomy file with hash IDs
    log_message(f"Writing taxonomy table to {taxonomy_output}", args.verbose)
    with open(taxonomy_output, "w") as f:
        f.write("TaxID\t" + "\t".join(RANKS) + "\n")
        for hash_id, tax_ranks in sorted(unique_taxonomies.items()):
            f.write(hash_id + "\t" + "\t".join([rank if rank else "NA" for rank in tax_ranks]) + "\n")
    
    # Write sequence to hash ID mapping for reference
    log_message(f"Writing sequence to hash ID mapping to {mapping_output}", args.verbose)
    with open(mapping_output, "w") as f:
        f.write("Sequence\tTaxID\n")
        for seq_id, hash_id in sorted(query_to_hash.items()):
            f.write(f"{seq_id}\t{hash_id}\n")
    
    # Final stats
    assignment_rate = len(query_taxonomy)/len(query_hits) if query_hits else 0
    log_message(f"Total reads processed: {len(query_hits)}", True)
    log_message(f"Assigned directly (1 hit): {n_single}", True)
    log_message(f"Assigned via LCA (multiple hits): {n_lca}", True)
    log_message(f"No taxonomy assignments: {n_unmatched}", True)
    log_message(f"Total assignments: {len(query_taxonomy)}", True)
    log_message(f"Assignment rate: {assignment_rate:.1%}", True)
    log_message(f"Final taxonomies: {len(unique_taxonomies)}", True)
    
    # Write final stats to log file
    with open(log_output, 'a') as log_file:
        log_file.write("\nFinal Statistics:\n")
        log_file.write(f"Total reads processed: {len(query_hits)}\n")
        log_file.write(f"Assigned directly (1 hit): {n_single}\n")
        log_file.write(f"Assigned via LCA (multiple hits): {n_lca}\n")
        log_file.write(f"No taxonomy assignments: {n_unmatched}\n")
        log_file.write(f"Total assignments: {len(query_taxonomy)}\n")
        log_file.write(f"Assignment rate: {assignment_rate:.1%}\n")
        log_file.write(f"Final taxonomies: {len(unique_taxonomies)}\n")
    
    log_message("Processing complete! Output files:", True)
    log_message(f"  - Abundance table: {abundance_output}", True)
    log_message(f"  - Taxonomy table: {taxonomy_output}", True)
    log_message(f"  - Sequence mapping: {mapping_output}", True)
    log_message(f"  - Log file: {log_output}", True)

if __name__ == "__main__":
    main()
