#!/usr/bin/env python3

"""
Script to extract sequences from regions around IgH genes with specified flanking distance.
"""

import argparse
import os
import subprocess
import tempfile

def merge_overlapping_regions(regions):
    """Merge overlapping genomic regions"""
    if not regions:
        return []
        
    # Sort regions by chromosome and start position
    sorted_regions = sorted(regions, key=lambda r: (r[0], r[1]))
    
    merged = [sorted_regions[0]]
    
    for region in sorted_regions[1:]:
        prev_region = merged[-1]
        
        # Check if current region overlaps with previous one (same chromosome)
        if region[0] == prev_region[0] and region[1] <= prev_region[2]:
            # Extend previous region if needed
            merged[-1] = (prev_region[0], prev_region[1], max(prev_region[2], region[2]))
        else:
            # Add as a new region
            merged.append(region)
    
    return merged

def main():
    parser = argparse.ArgumentParser(description="Extract sequences from regions around IgH genes")
    parser.add_argument('--igh_genes', required=True, help="BED file with IgH gene locations")
    parser.add_argument('--genome', required=True, help="Path to the genome FASTA file")
    parser.add_argument('--flanking', type=int, required=True, help="Flanking distance in bp")
    parser.add_argument('--output', required=True, help="Output FASTA file with extracted sequences")
    
    args = parser.parse_args()
    
    print(f"Reading IgH gene locations from {args.igh_genes}")
    igh_regions = []
    
    with open(args.igh_genes, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) < 3:
                continue
                
            chrom = fields[0]
            start = max(0, int(fields[1]) - args.flanking)
            end = int(fields[2]) + args.flanking
            
            igh_regions.append((chrom, start, end))
    
    # Merge overlapping regions
    merged_regions = merge_overlapping_regions(igh_regions)
    print(f"Found {len(igh_regions)} IgH regions, merged into {len(merged_regions)} non-overlapping regions")
    
    # Create a temporary BED file with the regions to extract
    with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as tmp:
        tmp_bed = tmp.name
        for region in merged_regions:
            tmp.write(f"{region[0]}\t{region[1]}\t{region[2]}\n")
    
    # Use bedtools to extract the sequences
    print(f"Extracting sequences...")
    cmd = f"bedtools getfasta -fi {args.genome} -bed {tmp_bed} -fo {args.output}"
    subprocess.run(cmd, shell=True, check=True)
    
    # Cleanup
    os.unlink(tmp_bed)
    
    print(f"Extracted sequences written to {args.output}")

if __name__ == "__main__":
    main() 