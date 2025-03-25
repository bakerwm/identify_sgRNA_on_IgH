#!/usr/bin/env python3

"""
Custom k-mer based repeat finder script.
This script analyzes a FASTA file to find repeated sequences using a k-mer approach.
"""

import argparse
import os
from collections import defaultdict
from Bio import SeqIO

def find_repeats(seq, min_length, min_copies):
    """
    Find repeats in a sequence using a k-mer approach.
    
    Args:
        seq: The DNA sequence to analyze
        min_length: Minimum length of repeats to find
        min_copies: Minimum number of copies required
        
    Returns:
        List of tuples (repeat_seq, positions, copy_count)
    """
    repeats = []
    
    # Start with the minimum length and search for repeats
    for k in range(min_length, min_length + 30):  # Check repeats up to min_length + 30 bp
        kmers = defaultdict(list)
        
        # Generate all k-mers and their positions
        for i in range(len(seq) - k + 1):
            kmer = seq[i:i+k]
            kmers[kmer].append(i)
        
        # Filter to keep only k-mers with sufficient copies
        for kmer, positions in kmers.items():
            if len(positions) >= min_copies:
                repeats.append((kmer, positions, len(positions)))
    
    return repeats

def process_fasta(fasta_file, min_length, min_copies, output_file):
    """Process a FASTA file to find repeats in each sequence"""
    all_repeats = []
    
    # Parse the FASTA file
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq_id = record.id
        sequence = str(record.seq).upper()
        
        print(f"Analyzing sequence: {seq_id} (length: {len(sequence)})")
        
        # Find repeats in this sequence
        seq_repeats = find_repeats(sequence, min_length, min_copies)
        
        # Add sequence identifier to each repeat
        for repeat_seq, positions, copy_count in seq_repeats:
            for pos in positions:
                all_repeats.append({
                    'seq_id': seq_id,
                    'start': pos,
                    'end': pos + len(repeat_seq),
                    'repeat_seq': repeat_seq,
                    'copy_count': copy_count
                })
    
    # Write results to BED file
    with open(output_file, 'w') as out:
        for repeat in all_repeats:
            # BED format: chrom start end name score strand
            out.write(f"{repeat['seq_id']}\t{repeat['start']}\t{repeat['end']}\t"
                     f"repeat_{len(repeat['repeat_seq'])}bp_x{repeat['copy_count']}\t"
                     f"{repeat['copy_count']}\t+\n")
    
    print(f"Found {len(all_repeats)} repeats")

def main():
    parser = argparse.ArgumentParser(description="Find repeats in DNA sequences using a k-mer approach")
    parser.add_argument('--input', required=True, help="Input FASTA file")
    parser.add_argument('--min_length', type=int, default=13, help="Minimum length of repeats to find")
    parser.add_argument('--min_copies', type=int, default=5, help="Minimum number of copies required")
    parser.add_argument('--output', required=True, help="Output BED file")
    
    args = parser.parse_args()
    
    process_fasta(args.input, args.min_length, args.min_copies, args.output)
    
    print(f"Results written to {args.output}")

if __name__ == "__main__":
    main() 