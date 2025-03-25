#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = "Ming Wang"
__email__ = "wangm08@hotmail.com"
__version__ = "1.0"
__date__ = "2025-03-25"

"""
Script to filter repeat sequences based on specified criteria:
- Minimum length > 13 bp
- Minimum copy number > 10

Output:
- Filtered BED file with repeat sequences that meet the criteria

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
"""

import argparse
import re
import os
import sys

def filter_repeats(input_file, min_length, min_copy_number, output_file):
    """
    Filter repeats based on criteria
    
    Args:
        input_file: Path to the input BED file
        min_length: Minimum length of repeats to keep
        min_copy_number: Minimum copy number of repeats to keep
        output_file: Path to the output filtered BED file
    """
    filtered_repeats = []
    total_count = 0
    
    with open(input_file, 'r') as f:
        for line in f:
            total_count += 1
            fields = line.strip().split('\t')
            
            if len(fields) < 6:
                continue
                
            # Extract repeat information
            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            name = fields[3]
            
            # Parse score field as copy number
            try:
                copy_number = int(fields[4])
            except ValueError:
                # Try to extract copy number from the name field if it's not in the score
                copy_match = re.search(r'x(\d+)', name)
                if copy_match:
                    copy_number = int(copy_match.group(1))
                else:
                    copy_number = 0
            
            strand = fields[5]
            
            # Calculate length
            length = end - start
            
            # Apply filters
            if length >= min_length and copy_number >= min_copy_number:
                filtered_repeats.append(fields)
    
    # Write filtered repeats to output file
    with open(output_file, 'w') as out:
        for fields in filtered_repeats:
            out.write('\t'.join(fields) + '\n')
    
    print(f"Total repeats: {total_count}")
    print(f"Filtered repeats (length >= {min_length}bp, copy number >= {min_copy_number}): {len(filtered_repeats)}")

def extract_chromosome_wide_coordinates(bed_record):
    """
    Extract chromosome-wide coordinates for a repeat sequence.
    
    Args:
        bed_record (list): A BED file record as a list of fields
    
    Returns:
        list: Modified BED record with chromosome-wide coordinates
    """
    # check if chr:start-end format
    pattern = r'.*:\d+-\d+$'
    if re.match(pattern, bed_record[0]):
        # Parse the chromosome region from the first field
        chrom_region = bed_record[0].split(':')
        chromosome = chrom_region[0]
        region_start, region_end = map(int, chrom_region[1].split('-'))
        
        # Parse the local repeat coordinates 
        local_start = int(bed_record[1])
        local_end = int(bed_record[2])
        
        # Calculate chromosome-wide coordinates
        chrom_wide_start = region_start + local_start
        chrom_wide_end = region_start + local_end
        
        # Create new record with chromosome-wide coordinates
        new_record = bed_record.copy()
        new_record[0] = chromosome
        new_record[1] = str(chrom_wide_start)
        new_record[2] = str(chrom_wide_end)
    else:
        new_record = bed_record.copy() # not changed
    
    return new_record

def fix_repeats_coordinates(input_file, output_file):
    """
    Fix the coordinates of the repeat sequences to be chromosome-wide.
    """
    with open(input_file, 'r') as f, open(output_file, 'w') as out:
        for line in f:
            fields = line.strip().split('\t') 
            fixed_record = extract_chromosome_wide_coordinates(fields)
            out.write('\t'.join(fixed_record) + '\n')

def main():
    parser = argparse.ArgumentParser(description="Filter repeat sequences based on criteria")
    parser.add_argument('--input', required=True, help="Input BED file with repeats")
    parser.add_argument('--min_length', type=int, default=13, help="Minimum length of repeats to keep")
    parser.add_argument('--min_copy_number', type=int, default=10, help="Minimum copy number of repeats to keep")
    parser.add_argument('--output', required=True, help="Output filtered BED file")
    
    args = parser.parse_args()
    
    # Check if input file exists
    if not os.path.exists(args.input):
        print(f"Error: Input file {args.input} does not exist", file=sys.stderr)
        sys.exit(1)
    
    # Filter the repeats by length and copy number
    filtered_repeats = args.output.replace('.filtered.bed', '.by_len_and_copy.bed')
    filter_repeats(args.input, args.min_length, args.min_copy_number, filtered_repeats)
    print(f"Filtered results written to {filtered_repeats}")

    # Fix the coordinates of the repeat sequences to be chromosome-wide
    fix_repeats_coordinates(filtered_repeats, args.output)
    print(f"Fixed coordinates written to {args.output}")

if __name__ == "__main__":
    import re  # Import here to avoid potential circular import
    main() 