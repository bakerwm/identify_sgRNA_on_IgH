#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = "Ming Wang"
__email__ = "wangm08@hotmail.com"
__version__ = "1.0"
__date__ = "2025-03-25"

"""
Script to combine repeat results from different tools into a unified format.

Format of TRF output:
see documentation: https://tandem.bu.edu/trf/definitions#table

Example:
Sequence: 12:111254829-116010093
584 644 20 3.3 20 70 25 71 11 31 39 18 1.86 GCCTGGGTCAGCATGGCCCT GCCTGGGTCAGCATGGCCCTGCCTGGGAAGCATGGCCTGGGTCAGCATGGCCCTGCCTGGG

column-1: start (584)
column-2: end (644)
column-3: repeat size (20)
column-4: copy number (3.3)
column-5: consensus size (20)
column-6: percent of matches (70)
column-7: percent of indels (25)
column-8: alignment score (71)
column-9-12: percent composition for each of the four nucleotides, A, C, G, T (11 31 39 18)
column-13: entropy measure based on percent composition (1.86)
column-14: consensus sequence (GCCTGGGTCAGCATGGCCCT)
column-15: repeat sequence (GCCTGGGTCAGCATGGCCCTGCCTGGGAAGCATGGCCTGGGTCAGCATGGCCCTGCCTGGG)

Source:
Table Explanation:
The summary table includes the following information:

Indices of the repeat relative to the start of the sequence.
Period size of the repeat.
Number of copies aligned with the consensus pattern.
Size of consensus pattern (may differ slightly from the period size).
Percent of matches between adjacent copies overall.
Percent of indels between adjacent copies overall.
Alignment score.
Percent composition for each of the four nucleotides.
Entropy measure based on percent composition.
"""

import argparse
import os
import re
import glob

def parse_trf_dat(trf_file, min_length=13, min_copy_number=10):
    """Parse Tandem Repeats Finder .dat output file"""
    repeats = []
    
    with open(trf_file, 'r') as f:
        current_seq = None
        
        for line in f:
            line = line.strip()
            
            # Extract sequence name
            if line.startswith('Sequence:'):
                current_seq = line.split()[1]
                continue
            
            # Skip header lines
            if not line or line.startswith('Parameters'):
                continue
                
            # Process repeat data lines
            if current_seq and re.match(r'\d+', line):
                fields = line.split()
                if len(fields) >= 15:
                    start = int(fields[0]) - 1 # Convert to 0-based for BED
                    end = int(fields[1])
                    repeat_length = int(fields[2])
                    copy_number = float(fields[3])
                    consensus_size = int(fields[4])
                    consensus_seq = ""
                    if re.match(r'^[ATGCN]+$', fields[13]) and len(fields[13]) == consensus_size:
                        consensus_seq = fields[13]
                    
                    # Only include repeats that meet our criteria
                    if repeat_length >= min_length and copy_number >= min_copy_number:
                        repeats.append({
                            'chrom': current_seq,
                            'start': start,
                            'end': end,
                            'name': f"TRF_repeat_{consensus_size}bp_x{copy_number:.1f}",
                            'score': int(copy_number),
                            'strand': '+',
                            'consensus_seq': consensus_seq
                        })
    
    return repeats

def parse_repeatmasker(rm_file, min_length=13, min_copy_number=10):
    """Parse RepeatMasker .out output file"""
    repeats = []
    
    with open(rm_file, 'r') as f:
        # Skip header lines
        for _ in range(3):
            next(f, None)
            
        for line in f:
            fields = line.strip().split()
            if len(fields) >= 14:
                try:
                    score = int(fields[0])
                    sequence_name = fields[4]
                    start = int(fields[5])
                    end = int(fields[6])
                    repeat_type = fields[9]
                    repeat_class = fields[10] if len(fields) > 10 else "Unknown"
                    
                    # Convert RepeatMasker coordinates to BED format
                    # RepeatMasker output is 1-based, BED is 0-based
                    repeats.append({
                        'chrom': sequence_name,
                        'start': start - 1,
                        'end': end,
                        'name': f"RM_{repeat_type}_{repeat_class}",
                        'score': score,
                        'strand': '+'
                    })
                except (ValueError, IndexError):
                    continue
    
    return repeats

def parse_kmer_bed(kmer_file, min_length=13, min_copy_number=10):
    """Parse BED file from k-mer repeat finder"""
    repeats = []
    
    with open(kmer_file, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) >= 6:
                # BED format: chrom start end name score strand
                repeats.append({
                    'chrom': fields[0],
                    'start': int(fields[1]),
                    'end': int(fields[2]),
                    'name': f"KMER_{fields[3]}",
                    'score': int(fields[4]),
                    'strand': fields[5]
                })
    
    return repeats

def parse_misa(misa_file, min_length=13, min_copy_number=10):
    """Parse MISA output file for microsatellites"""
    repeats = []
    
    with open(misa_file, 'r') as f:
        next(f)  # Skip header
        
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) >= 6:
                sequence_id = fields[0]
                ssr_nr = fields[1]
                ssr_type = fields[2]
                ssr_size = fields[3]
                start = int(fields[5])
                end = start + int(ssr_size) - 1
                
                # Estimate copy number from SSR type and size
                unit_size = int(ssr_type.split('-')[0])
                copy_number = int(ssr_size) // unit_size
                
                if int(ssr_size) >= min_length and copy_number >= min_copy_number:
                    repeats.append({
                        'chrom': sequence_id,
                        'start': start - 1,  # Convert to 0-based for BED
                        'end': end,
                        'name': f"MISA_{ssr_type}_x{copy_number}",
                        'score': copy_number,
                        'strand': '+'
                    })
    
    return repeats

def combine_and_filter_repeats(repeats, min_length=13, min_copy_number=10):
    """Combine repeats and filter based on criteria"""
    filtered = []
    
    for repeat in repeats:
        length = repeat['end'] - repeat['start']
        copy_number = repeat['score']
        
        if length >= min_length and copy_number >= min_copy_number:
            filtered.append(repeat)
    
    return filtered

def write_bed(repeats, output_file):
    """Write repeats to BED file"""
    with open(output_file, 'w') as out:
        for repeat in repeats:
            consensus = repeat.get('consensus_seq', '')
            out.write(f"{repeat['chrom']}\t{repeat['start']}\t{repeat['end']}\t"
                     f"{repeat['name']}\t{repeat['score']}\t{repeat['strand']}\t{consensus}\n")

def main():
    parser = argparse.ArgumentParser(description="Combine repeat results from different tools")
    parser.add_argument('--trf', help="TRF output .dat file")
    parser.add_argument('--repeatmasker', help="RepeatMasker output .out file")
    parser.add_argument('--kmer', help="K-mer repeat finder BED output", required=False)
    parser.add_argument('--misa', help="MISA output file", required=False)
    parser.add_argument('--min_length', type=int, default=13, help="Minimum repeat length")
    parser.add_argument('--min_copy_number', type=int, default=10, help="Minimum copy number")
    parser.add_argument('--output_bed', required=True, help="Output BED file")
    
    args = parser.parse_args()
    
    all_repeats = []
    
    # Parse results from each tool
    if args.trf:
        trf_repeats = parse_trf_dat(args.trf)
        print(f"Found {len(trf_repeats)} repeats from TRF")
        all_repeats.extend(trf_repeats)
    
    if args.repeatmasker:
        # Handle both single file and glob pattern
        rm_files = glob.glob(args.repeatmasker)
        for rm_file in rm_files:
            rm_repeats = parse_repeatmasker(rm_file)
            print(f"Found {len(rm_repeats)} repeats from RepeatMasker in {rm_file}")
            all_repeats.extend(rm_repeats)
    
    if args.kmer:
        kmer_repeats = parse_kmer_bed(args.kmer)
        print(f"Found {len(kmer_repeats)} repeats from k-mer analysis")
        all_repeats.extend(kmer_repeats)
    
    if args.misa:
        # Handle both single file and glob pattern
        misa_files = glob.glob(args.misa)
        for misa_file in misa_files:
            misa_repeats = parse_misa(misa_file)
            print(f"Found {len(misa_repeats)} repeats from MISA in {misa_file}")
            all_repeats.extend(misa_repeats)
    
    # Filter based on criteria
    # print(">>>", f'{all_repeats[0]}')
    filtered_repeats = combine_and_filter_repeats(
        all_repeats, args.min_length, args.min_copy_number)
    
    print(f"Total: {len(all_repeats)} repeats, filtered to {len(filtered_repeats)} repeats")
    
    # Write to output file
    write_bed(filtered_repeats, args.output_bed)
    print(f"Results written to {args.output_bed}")

if __name__ == "__main__":
    main() 