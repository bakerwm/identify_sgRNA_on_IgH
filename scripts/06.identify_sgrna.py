#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = "Ming Wang"
__email__ = "wangm08@hotmail.com"
__version__ = "1.0"
__date__ = "2025-03-25"

"""
Script to identify sgRNA sequences from specified sequence
Criteria:
- 23 bp in length
- Ending with GG (for Cas9/Sp sgRNA)
- Consider both forward and reverse-complement strands

Input:
- Filtered repeats in BED format

Example:
12      111410728       111411078       TRF_repeat_23bp_x15.6   15      +       GGTTAGTCCTCAGTTAGTCCTCA

columns:
1. chromosome
2. start
3. end
4. repeat_name
5. copy_number
6. strand
7. consensus_sequence

Output:
- sgRNA candidates in BED format

Example:
12      111410730       111411079       TRF_repeat_34bp_x10.4   10      +       TTAGTCCTCAGGTTAGTCCTCAGTTAGTCCTCAG      rev     sgRNA_000001    ACTAACTGAGGACTAACCTGAGG

columns:
1. chromosome
2. start
3. end
4. repeat_name
5. copy_number
6. strand
7. consensus_sequence
8. sgRNA_direction
9. sgRNA_id
10. sgRNA_sequence
"""

import argparse
import os
import sys
import pathlib
from Bio import Seq
from Bio.Seq import Seq

def identify_sgRNA(x, strand='fwd', size=23, suffix='GG'):
    """
    Identify sgRNA from DNA sequence

    Args:
        x: str, Consensus sequence
        strand: str, Specify strand to extract sgRNA from consensus sequence, default: 'fwd'
        size: int, Length of the sgRNA, default: 23
        suffix: str, Suffix of the sgRNA, default: 'GG' for SpCas9/Cas9
    """
    # extract sgRNA seq
    sgRNA_list = []
    if not isinstance(x, str):
        return sgRNA_list
    if len(x) < size:
        return sgRNA_list
    if strand == 'rev':
        x = str(Seq(x).reverse_complement())
    # extract sgRNA from consensus sequence
    for i in range(len(x) - size + 1):
        window = x[i:i+size]
        if window.endswith(suffix):
            sgRNA_list.append(window)

    return sgRNA_list

def identify_sgrna_candidates(input_file, output_file, size=23, suffix='GG'):
    """
    Identify sgRNA candidates from filtered repeats

    input_file:
        - Consensus sequence in column 7
    
    Args:
        input_file: Path to the input BED file with filtered repeats
        output_file: Path to the output file containing sgRNA candidates
    """
    total_sgRNAs = 0
    total_repeats = 0
    output_dir = pathlib.Path(output_file).parent
    output_dir.mkdir(parents=True, exist_ok=True)
    
    with open(input_file, 'r') as f, open(output_file, 'w') as out:
        for line in f:
            total_repeats += 1
            fields = line.strip().split('\t')
            
            if len(fields) < 7:
                continue
            # Extract repeat information: fwd (forward strand)
            sgRNA_list = identify_sgRNA(fields[6], 'fwd', size, suffix)
            if len(sgRNA_list) > 0:
                for sgRNA in sgRNA_list:
                    total_sgRNAs += 1
                    out.write('\t'.join(fields + ['fwd', f'sgRNA_{total_sgRNAs:06d}', sgRNA]) + '\n')

            # Extract repeat information: rev (reverse strand)
            sgRNA_list = identify_sgRNA(fields[6], 'rev', size, suffix)
            if len(sgRNA_list) > 0:
                for sgRNA in sgRNA_list:
                    total_sgRNAs += 1
                    out.write('\t'.join(fields + ['rev', f'sgRNA_{total_sgRNAs:06d}', sgRNA]) + '\n')
    
    print(f"Total repeats: {total_repeats}")
    print(f"Total sgRNA: {total_sgRNAs}")

def main():
    parser = argparse.ArgumentParser(description="Extract sgRNA candidates from filtered repeats")
    parser.add_argument('--input', required=True, help="Input BED file with filtered repeats")
    parser.add_argument('--output', required=True, help="Output file for sgRNA candidates")
    parser.add_argument('--size', type=int, default=23, help="Length of the sgRNA, default: 23")
    parser.add_argument('--suffix', type=str, default='GG', help="Suffix of the sgRNA, default: 'GG'")
    args = parser.parse_args()
    
    # Check if input file exists
    if not os.path.exists(args.input):
        print(f"Error: Input file {args.input} does not exist", file=sys.stderr)
        sys.exit(1)
    
    # Identify sgRNA candidates
    identify_sgrna_candidates(args.input, args.output, args.size, args.suffix)

if __name__ == "__main__":
    main() 
