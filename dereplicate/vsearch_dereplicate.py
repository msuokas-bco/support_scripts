#!/usr/bin/env python3
# vsearch_dereplicate.py v1.1.0
# Changelog:
#   v1.1.0 - Robust sample ID extraction using regex (supports dots in filenames)
#            Added --strip_suffix option to remove trailing suffixes (e.g. _trimmed)
#              from sample IDs, so they match metadata without post-processing
#            SHA1 replaced with SHA256 for sequence hashing
#            vsearch stderr suppressed during normal runs, shown only on failure
#            Output parent directories created automatically if they do not exist
#            Python 3.9+ version guard added at startup

import os
import re
import sys
import subprocess
import tempfile
import argparse
import hashlib
import pandas as pd
from collections import defaultdict
import gzip
from Bio import SeqIO
from concurrent.futures import ProcessPoolExecutor

if sys.version_info < (3, 9):
    sys.exit("Python 3.9 or later is required (uses str.removesuffix)")

def parse_args():
    parser = argparse.ArgumentParser(description='Dereplicate sequences from multiple FASTQ files')
    parser.add_argument('--input_dir', required=True, help='Directory containing fastq.gz files')
    parser.add_argument('--output_fasta', required=True, help='Output FASTA file with dereplicated sequences')
    parser.add_argument('--output_table', required=True, help='Output TSV table with sequence abundances')
    parser.add_argument('--derep_prefix', action='store_true', help='Use prefix dereplication instead of full-length')
    parser.add_argument('--min_seq_length', type=int, default=1, help='Minimum sequence length')
    parser.add_argument('--min_unique_size', type=int, default=1, help='Minimum unique sequence abundance')
    parser.add_argument('--threads', type=int, default=1, help='Number of parallel samples to process')
    parser.add_argument('--strip_suffix', default='', help='Suffix to strip from sample IDs after extension removal (e.g. _trimmed)')
    return parser.parse_args()

def extract_sample_id(fastq_file, strip_suffix=''):
    """Extract sample ID from filename.

    Uses regex to strip .fastq or .fastq.gz extension, preserving any dots
    in the sample name itself (e.g. sample.1_trimmed.fastq.gz -> sample.1_trimmed).
    If strip_suffix is set (e.g. '_trimmed'), it is removed from the result
    (e.g. sample.1_trimmed -> sample.1).
    """
    sample_id = re.sub(r'\.fastq(\.gz)?$', '', os.path.basename(fastq_file))
    if strip_suffix:
        sample_id = sample_id.removesuffix(strip_suffix)
    return sample_id

def convert_fastq_to_fasta(fastq_file, output_dir, strip_suffix=''):
    """Convert a fastq.gz file to fasta format with sample ID in the header"""
    sample_id = extract_sample_id(fastq_file, strip_suffix)
    output_fasta = os.path.join(output_dir, f"{sample_id}.fasta")
    
    with gzip.open(fastq_file, "rt") as handle, open(output_fasta, "w") as output:
        for record in SeqIO.parse(handle, "fastq"):
            output.write(f">{sample_id}_{record.id}\n{str(record.seq)}\n")
    
    return output_fasta

def process_sample(fastq_file, temp_dir, args):
    """Process a single FASTQ file"""
    fasta_file = convert_fastq_to_fasta(fastq_file, temp_dir, args.strip_suffix)
    sample_id = extract_sample_id(fastq_file, args.strip_suffix)
    
    derep_fasta = os.path.join(temp_dir, f"{sample_id}_derep.fasta")
    uc_file = os.path.join(temp_dir, f"{sample_id}_derep.uc")
    
    cmd = [
        'vsearch',
        '--derep_prefix' if args.derep_prefix else '--derep_fulllength',
        fasta_file,
        '--output', derep_fasta,
        '--uc', uc_file,
        '--xsize',
        '--minseqlength', str(args.min_seq_length),
        '--minuniquesize', str(args.min_unique_size),
        '--fasta_width', '0'
    ]
    
    try:
        subprocess.run(cmd, check=True, stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        print(e.stderr.decode() if e.stderr else '', file=sys.stderr)
        raise
    
    return sample_id, derep_fasta, uc_file

def parse_uc_file(uc_file):
    """Parse UC file to get sequence clusters and sizes"""
    seq_counts = defaultdict(int)
    
    with open(uc_file, 'r') as f:
        for line in f:
            if line.startswith('H'):
                parts = line.strip().split('\t')
                target_seq = parts[9].split(';')[0]
                seq_counts[target_seq] += 1
            elif line.startswith('S'):
                parts = line.strip().split('\t')
                seq_id = parts[8].split(';')[0]
                seq_counts[seq_id] += 1
    
    return seq_counts

def generate_sha256(sequence):
    """Generate SHA256 hash for a sequence"""
    return hashlib.sha256(sequence.encode()).hexdigest()

def main():
    args = parse_args()

    os.makedirs(os.path.dirname(os.path.abspath(args.output_fasta)), exist_ok=True)
    os.makedirs(os.path.dirname(os.path.abspath(args.output_table)), exist_ok=True)

    with tempfile.TemporaryDirectory() as temp_dir:
        fastq_files = [os.path.join(args.input_dir, f) for f in os.listdir(args.input_dir) 
                      if f.endswith('.fastq.gz')]
        
        if not fastq_files:
            print("No .fastq.gz files found in the input directory!")
            return
        
        print(f"Found {len(fastq_files)} FASTQ files for processing")
        
        sample_results = []
        with ProcessPoolExecutor(max_workers=args.threads) as executor:
            futures = [executor.submit(process_sample, f, temp_dir, args) for f in fastq_files]
            for future in futures:
                sample_results.append(future.result())
        
        abundance_data = defaultdict(lambda: defaultdict(int))
        sequence_data = {}
        
        for sample_id, derep_fasta, uc_file in sample_results:
            seq_counts = parse_uc_file(uc_file)
            
            with open(derep_fasta, 'r') as f:
                for record in SeqIO.parse(f, 'fasta'):
                    seq_str = str(record.seq)
                    seq_id = generate_sha256(seq_str)
                    
                    if seq_id not in sequence_data:
                        sequence_data[seq_id] = seq_str
                    
                    original_id = record.id
                    
                    # Always add count (zero if not present)
                    abundance_data[seq_id][sample_id] += seq_counts.get(original_id, 0)
        
        # Write dereplicated sequences to FASTA
        with open(args.output_fasta, 'w') as f:
            for seq_id, seq in sequence_data.items():
                f.write(f">{seq_id}\n{seq}\n")
        
        sample_ids = sorted([extract_sample_id(f, args.strip_suffix) for f in fastq_files])
        table_data = []
        
        for seq_id in sequence_data:
            row = {'sequence_id': seq_id}
            for sample in sample_ids:
                row[sample] = abundance_data[seq_id].get(sample, 0)
            table_data.append(row)
        
        df = pd.DataFrame(table_data)
        df.to_csv(args.output_table, sep='\t', index=False)
        
        print(f"Dereplication complete!")
        print(f"Found {len(sequence_data)} unique sequences across {len(sample_ids)} samples")
        print(f"Results written to:")
        print(f"  - Sequences: {args.output_fasta}")
        print(f"  - Abundance table: {args.output_table}")

if __name__ == "__main__":
    main()

