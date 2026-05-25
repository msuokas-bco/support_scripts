#!/usr/bin/env python3
# alignment_to_taxa.py v1.3.1
# Changelog:
#   v1.3.1 - BAM opened as context manager so the handle is closed safely on exceptions
#            Sort-order guard: raises ValueError if BAM header declares SO:coordinate
#            Score filter guard: when max_score == 0, all hits are passed to LCA
#              instead of being filtered by a threshold of zero (which admits everything)
#   v1.3.0 - Streaming BAM mode: query_hits dict eliminated entirely.
#            BAM scan and taxonomy assignment merged into one pass using
#              itertools.groupby over a name-sorted BAM. Peak RAM is now
#              taxonomy_dict + abundance + one query's hits instead of all hits.
#            Requires name-sorted BAM (samtools sort -n).
#   v1.2.1 - Use dict.get() for single-lookup taxonomy fetch in BAM scan
#            Remove redundant sample_header init and header alias in abundance loading
#   v1.2.0 - RAM and temporary-object optimizations:
#              Inline taxonomy lookup during BAM scan (stored strings are shared
#                references to taxonomy_dict values, not new per-hit allocations)
#              all_refs set only built when --verbose is active
#              abundance stored as list-per-row dict instead of nested defaultdict
#              query_taxonomy stores semicolon-joined strings; dedup via set()
#              merged_counts inner structure changed from defaultdict to list
#              query_hits, abundance, query_taxonomy freed with del after last use
#              Full keys-list copy eliminated; islice-based batching over dict iter
#   v1.1.0 - Fixed log_message to respect --verbose flag (was always printing)
#            Supplementary alignments (0x800, chimeric reads) now discarded to
#              prevent false multi-mapping signals corrupting LCA assignment
#            LCA denominator changed to total reference count so references with
#              no annotation at a rank count against the threshold, preventing
#              over-assignment at deeper ranks
#            MD5 replaced with SHA256 for taxonomy hash IDs
#            Batch progress uses ceiling division for consistent update count
"""
Taxonomy Assignment Tool

This script processes BAM alignments and taxonomy files to assign taxonomy to sequences
and create abundance tables with hash-based identifiers.

Usage:
    ./alignment_to_taxa.py --bam INPUT.bam --abundance ABUNDANCE.tsv --taxonomy TAXONOMY.tsv
                           [--output-dir OUTPUT_DIR] [--score-threshold 0.97] [--lca-threshold 0.8]
"""

import pysam
from collections import defaultdict, Counter
import itertools
import os
import hashlib
import argparse
from datetime import datetime

def parse_arguments():
    parser = argparse.ArgumentParser(description='Assign taxonomy to sequences and create abundance tables')

    # Required arguments
    parser.add_argument('--bam', required=True, help='Input BAM file with alignments (must be name-sorted: samtools sort -n)')
    parser.add_argument('--abundance', required=True, help='Tab-separated abundance file')
    parser.add_argument('--taxonomy', required=True, help='Tab-separated taxonomy file')

    # Optional arguments
    parser.add_argument('--output-dir', default='.', help='Output directory (default: current directory)')
    parser.add_argument('--score-threshold', type=float, default=0.97,
                        help='Minimum fraction of best alignment score for retaining secondary hits. '
                             'Hits below this fraction are discarded; if only one hit remains it is '
                             'assigned directly, otherwise LCA is used (default: 0.97)')
    parser.add_argument('--lca-threshold', type=float, default=0.8,
                        help='Agreement required at each rank to retain it (default: 0.8)')
    parser.add_argument('--prefix', default='',
                        help='Prefix for output files (default: none)')
    parser.add_argument('--verbose', action='store_true', help='Print detailed progress')

    return parser.parse_args()

def log_message(message, verbose=False):
    """Print log message with timestamp. Only prints if verbose is True."""
    if verbose:
        timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        print(f"[{timestamp}] {message}", flush=True)

def hash_taxonomy(taxonomy):
    """Generate hash ID for a taxonomy string"""
    return hashlib.sha256(taxonomy.encode()).hexdigest()[:10]

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
        next(f)  # skip header line (ReferenceID\tTaxonomy)
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) == 2:
                taxonomy_dict[parts[0]] = parts[1]
    log_message(f"Loaded {len(taxonomy_dict)} taxonomy entries", args.verbose)

    # Load abundance data as list-per-row: seq_id -> list[int] parallel to sample_header.
    # A plain dict of lists avoids the per-row nested-dict overhead of defaultdict(defaultdict(int)).
    log_message(f"Loading abundance data from {args.abundance}", args.verbose)
    abundance = {}
    with open(args.abundance) as f:
        sample_header = f.readline().strip().split("\t")[1:]
        for line in f:
            parts = line.strip().split("\t")
            seq_id = parts[0]
            counts = list(map(int, parts[1:]))
            if len(counts) != len(sample_header):
                raise ValueError(f"Abundance file row '{seq_id}' has {len(counts)} values, expected {len(sample_header)}")
            abundance[seq_id] = counts
    log_message(f"Loaded abundance data for {len(abundance)} sequences", args.verbose)

    def get_lca(taxa_list):
        """Compute the lowest common ancestor for a list of taxonomies.
        Returns a semicolon-joined string with exactly len(RANKS) fields."""
        parsed_taxa = []
        for tax in taxa_list:
            ranks = tax.split(";")
            while len(ranks) < len(RANKS):
                ranks.append("")
            parsed_taxa.append(ranks)

        lca_ranks = []
        for i in range(len(RANKS)):
            rank_values = [taxa[i] for taxa in parsed_taxa if taxa[i]]
            if not rank_values:
                lca_ranks.append("")
                continue
            # Denominator is total number of references (including those with no
            # annotation at this rank) so that missing annotations count against
            # the threshold and prevent over-assignment at deeper ranks.
            counts = Counter(rank_values)
            most_common, count = counts.most_common(1)[0]
            if count / len(parsed_taxa) >= args.lca_threshold:
                lca_ranks.append(most_common)
            else:
                lca_ranks.extend([""] * (len(RANKS) - i))
                break

        return ";".join(lca_ranks)

    # Stream BAM grouped by query name.
    # Requires a name-sorted BAM so all hits for a query are consecutive.
    # query_hits is never accumulated — each query's hits are processed and
    # discarded immediately, keeping peak RAM to taxonomy_dict + abundance +
    # one query's hits at a time.
    log_message(f"Processing alignments and assigning taxonomy from {args.bam}", args.verbose)

    query_taxonomy = {}
    n_single = 0
    n_lca = 0
    n_unmatched = 0
    total_queries = 0
    total_refs_checked = 0
    tax_lookup_failures = 0
    n_supplementary = 0
    all_refs = set() if args.verbose else None

    def iter_mapped(bam):
        nonlocal n_supplementary
        for read in bam:
            if read.is_supplementary:
                n_supplementary += 1
                continue
            if not read.is_unmapped:
                yield read

    with pysam.AlignmentFile(args.bam, "rb") as bamfile:
        so = dict(bamfile.header.get("HD", {})).get("SO", "")
        if so == "coordinate":
            raise ValueError(
                f"BAM is coordinate-sorted (SO:coordinate). "
                f"Name-sort it first with: samtools sort -n -o sorted.bam {args.bam}"
            )

        for query_name, group in itertools.groupby(iter_mapped(bamfile), key=lambda r: r.query_name):
            hits = []
            for read in group:
                ref = bamfile.get_reference_name(read.reference_id)
                score = read.get_tag("AS") if read.has_tag("AS") else read.mapping_quality / 100.0
                tax = taxonomy_dict.get(ref)
                if all_refs is not None:
                    all_refs.add(ref)
                hits.append((tax, score))

            total_queries += 1
            max_score = max(score for _, score in hits)
            if max_score > 0:
                filtered = [(tax, score) for tax, score in hits if score >= args.score_threshold * max_score]
            else:
                # All hits tied at zero; pass all to LCA rather than filtering by a degenerate threshold.
                filtered = hits

            taxonomies = []
            for tax, _ in filtered:
                total_refs_checked += 1
                if tax is None:
                    tax_lookup_failures += 1
                elif tax:
                    taxonomies.append(tax)

            if not taxonomies:
                n_unmatched += 1
            elif len(taxonomies) == 1:
                parts = taxonomies[0].split(";")
                while len(parts) < len(RANKS):
                    parts.append("")
                query_taxonomy[query_name] = ";".join(parts)
                n_single += 1
            else:
                lca_str = get_lca(taxonomies)
                if any(lca_str.split(";")):
                    query_taxonomy[query_name] = lca_str
                    n_lca += 1
                else:
                    n_unmatched += 1

            if args.verbose and total_queries % 1_000_000 == 0:
                log_message(f"Processed {total_queries:,} queries", args.verbose)

    log_message(f"Processed {total_queries:,} queries with alignments", args.verbose)
    log_message(f"Chimeric reads discarded (supplementary alignments): {n_supplementary}", args.verbose)

    if args.verbose:
        log_message(f"Found {len(all_refs)} unique reference sequences", args.verbose)
        test_refs = list(all_refs)[:min(100, len(all_refs))]
        ref_in_tax = sum(1 for ref in test_refs if ref in taxonomy_dict)
        log_message(f"Direct ref matches in taxonomy: {ref_in_tax}/{len(test_refs)}", args.verbose)
        if ref_in_tax < len(test_refs) * 0.5:
            log_message("Warning: Many references don't match taxonomy keys directly", args.verbose)
    all_refs = None

    # Print and log taxonomy lookup stats
    lookup_failure_rate = tax_lookup_failures / total_refs_checked if total_refs_checked > 0 else 0
    log_message(f"Total reference lookups: {total_refs_checked}", args.verbose)
    log_message(f"Taxonomy lookup failures: {tax_lookup_failures} ({lookup_failure_rate:.1%})", args.verbose)

    with open(log_output, 'a') as log_file:
        log_file.write(f"Chimeric reads discarded (supplementary alignments): {n_supplementary}\n")
        log_file.write(f"Total reference lookups: {total_refs_checked}\n")
        log_file.write(f"Taxonomy lookup failures: {tax_lookup_failures} ({lookup_failure_rate:.1%})\n")

    # Build hash IDs from the deduplicated set of taxonomy strings — one hash
    # computation per unique taxonomy instead of one per query.
    log_message("Creating hash IDs for taxonomies", args.verbose)
    tax_to_hash = {ts: hash_taxonomy(ts) for ts in set(query_taxonomy.values())}

    # Map queries to hash IDs, then free query_taxonomy
    query_to_hash = {query: tax_to_hash[ts] for query, ts in query_taxonomy.items()}
    del query_taxonomy

    # Build unique taxonomy rank lists from the already-deduplicated tax_to_hash
    unique_taxonomies = {}
    for tax_string, hash_id in tax_to_hash.items():
        tax_ranks = tax_string.split(";")
        while len(tax_ranks) < len(RANKS):
            tax_ranks.append("")
        unique_taxonomies[hash_id] = tax_ranks

    # Merge counts by hash ID.
    # Inner structure is a pre-allocated list (parallel to sample_header) rather
    # than a nested defaultdict to avoid per-taxonomy dict overhead.
    log_message("Merging counts by taxonomy", args.verbose)
    n_samples = len(sample_header)
    merged_counts = defaultdict(lambda: [0] * n_samples)
    for seq_id, hash_id in query_to_hash.items():
        if seq_id in abundance:
            row = abundance[seq_id]
            dest = merged_counts[hash_id]
            for i in range(n_samples):
                dest[i] += row[i]
    del abundance  # no longer needed

    # Write hash-based abundance table
    log_message(f"Writing abundance table to {abundance_output}", args.verbose)
    with open(abundance_output, "w") as f:
        f.write("TaxID\t" + "\t".join(sample_header) + "\n")
        for hash_id in sorted(merged_counts):
            counts = [str(merged_counts[hash_id][i]) for i in range(n_samples)]
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
    assignment_rate = len(query_to_hash) / total_queries if total_queries > 0 else 0
    log_message(f"Total reads processed: {total_queries}", True)
    log_message(f"Assigned directly (1 hit): {n_single}", True)
    log_message(f"Assigned via LCA (multiple hits): {n_lca}", True)
    log_message(f"No taxonomy assignments: {n_unmatched}", True)
    log_message(f"Total assignments: {len(query_to_hash)}", True)
    log_message(f"Assignment rate: {assignment_rate:.1%}", True)
    log_message(f"Final taxonomies: {len(unique_taxonomies)}", True)

    # Write final stats to log file
    with open(log_output, 'a') as log_file:
        log_file.write("\nFinal Statistics:\n")
        log_file.write(f"Total reads processed: {total_queries}\n")
        log_file.write(f"Assigned directly (1 hit): {n_single}\n")
        log_file.write(f"Assigned via LCA (multiple hits): {n_lca}\n")
        log_file.write(f"No taxonomy assignments: {n_unmatched}\n")
        log_file.write(f"Total assignments: {len(query_to_hash)}\n")
        log_file.write(f"Assignment rate: {assignment_rate:.1%}\n")
        log_file.write(f"Final taxonomies: {len(unique_taxonomies)}\n")

    log_message("Processing complete! Output files:", True)
    log_message(f"  - Abundance table: {abundance_output}", True)
    log_message(f"  - Taxonomy table: {taxonomy_output}", True)
    log_message(f"  - Sequence mapping: {mapping_output}", True)
    log_message(f"  - Log file: {log_output}", True)

if __name__ == "__main__":
    main()
