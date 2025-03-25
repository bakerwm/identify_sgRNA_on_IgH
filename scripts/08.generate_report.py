#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = "Ming Wang"
__email__ = "wangm08@hotmail.com"
__version__ = "1.0"
__date__ = "2025-03-25"

"""
Script to generate a summary report of the identified repeat sequences.
"""

import argparse
import os
import sys
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict, Counter

def load_bed_file(bed_file):
    """Load a BED file and return a list of entries"""
    entries = []
    
    with open(bed_file, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) >= 6:
                entries.append({
                    'chrom': fields[0],
                    'start': int(fields[1]),
                    'end': int(fields[2]),
                    'name': fields[3],
                    'score': int(fields[4]),
                    'strand': fields[5]
                })
    
    return entries

def calculate_distance_to_igh(repeats, igh_genes):
    """Calculate distance from each repeat to the nearest IgH gene"""
    # Process IgH genes
    igh_regions = []
    for gene in igh_genes:
        igh_regions.append((gene['chrom'], gene['start'], gene['end']))
    
    # Calculate distances
    for repeat in repeats:
        min_distance = float('inf')
        
        for chrom, start, end in igh_regions:
            if repeat['chrom'] == chrom:
                # Calculate distance to this IgH gene
                if repeat['end'] < start:
                    # Repeat is upstream of the gene
                    distance = start - repeat['end']
                elif repeat['start'] > end:
                    # Repeat is downstream of the gene
                    distance = repeat['start'] - end
                else:
                    # Repeat overlaps with the gene
                    distance = 0
                
                min_distance = min(min_distance, distance)
        
        # Add distance to repeat info
        repeat['distance_to_igh'] = min_distance if min_distance != float('inf') else None

def generate_statistics(repeats):
    """Generate statistics on the repeats"""
    stats = {
        'total_count': len(repeats),
        'length_distribution': [],
        'copy_number_distribution': [],
        'distance_distribution': [],
        'repeat_types': Counter(),
        'chromosomes': Counter()
    }
    
    for repeat in repeats:
        length = repeat['end'] - repeat['start']
        copy_number = repeat['score']
        distance = repeat['distance_to_igh']
        
        stats['length_distribution'].append(length)
        stats['copy_number_distribution'].append(copy_number)
        
        if distance is not None:
            stats['distance_distribution'].append(distance)
        
        # Extract repeat type from name
        repeat_type = repeat['name'].split('_')[0]
        stats['repeat_types'][repeat_type] += 1
        
        # Count chromosomes
        stats['chromosomes'][repeat['chrom']] += 1
    
    return stats

def plot_statistics(stats, output_dir):
    """Create plots visualizing the statistics"""
    # Length distribution
    plt.figure(figsize=(10, 6))
    plt.hist(stats['length_distribution'], bins=30)
    plt.xlabel('Repeat Length (bp)')
    plt.ylabel('Count')
    plt.title('Distribution of Repeat Lengths')
    plt.savefig(f"{output_dir}/length_distribution.png")
    
    # Copy number distribution
    plt.figure(figsize=(10, 6))
    plt.hist(stats['copy_number_distribution'], bins=30)
    plt.xlabel('Copy Number')
    plt.ylabel('Count')
    plt.title('Distribution of Copy Numbers')
    plt.savefig(f"{output_dir}/copy_number_distribution.png")
    
    # Distance to IgH genes
    if stats['distance_distribution']:
        plt.figure(figsize=(10, 6))
        plt.hist(stats['distance_distribution'], bins=30)
        plt.xlabel('Distance to nearest IgH gene (bp)')
        plt.ylabel('Count')
        plt.title('Distribution of Distances to IgH Genes')
        plt.savefig(f"{output_dir}/distance_distribution.png")
    
    # Repeat types pie chart
    plt.figure(figsize=(10, 6))
    labels = list(stats['repeat_types'].keys())
    values = list(stats['repeat_types'].values())
    plt.pie(values, labels=labels, autopct='%1.1f%%')
    plt.title('Repeat Types')
    plt.savefig(f"{output_dir}/repeat_types.png")
    
    # Chromosome distribution
    plt.figure(figsize=(12, 6))
    chrom_keys = list(stats['chromosomes'].keys())
    chrom_values = list(stats['chromosomes'].values())
    plt.bar(chrom_keys, chrom_values)
    plt.xlabel('Chromosome')
    plt.ylabel('Number of Repeats')
    plt.title('Repeat Distribution by Chromosome')
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.savefig(f"{output_dir}/chromosome_distribution.png")

def write_report(stats, repeats, output_file):
    """Write a text report summarizing the findings"""
    with open(output_file, 'w') as f:
        f.write("# Repeat Sequence Analysis Report\n\n")
        
        f.write(f"## Summary Statistics\n\n")
        f.write(f"Total repeats found: {stats['total_count']}\n")
        f.write(f"Average repeat length: {np.mean(stats['length_distribution']):.2f} bp\n")
        f.write(f"Average copy number: {np.mean(stats['copy_number_distribution']):.2f}\n")
        
        f.write("\n## Repeat Types\n\n")
        for repeat_type, count in stats['repeat_types'].most_common():
            f.write(f"- {repeat_type}: {count} ({count/stats['total_count']*100:.1f}%)\n")
        
        f.write("\n## Top 10 Longest Repeats\n\n")
        sorted_by_length = sorted(repeats, key=lambda r: r['end'] - r['start'], reverse=True)
        f.write("| Chromosome | Start | End | Length | Copy Number | Name |\n")
        f.write("|------------|-------|-----|--------|-------------|------|\n")
        for repeat in sorted_by_length[:10]:
            length = repeat['end'] - repeat['start']
            f.write(f"| {repeat['chrom']} | {repeat['start']} | {repeat['end']} | {length} | {repeat['score']} | {repeat['name']} |\n")
        
        f.write("\n## Top 10 Highest Copy Number Repeats\n\n")
        sorted_by_copies = sorted(repeats, key=lambda r: r['score'], reverse=True)
        f.write("| Chromosome | Start | End | Length | Copy Number | Name |\n")
        f.write("|------------|-------|-----|--------|-------------|------|\n")
        for repeat in sorted_by_copies[:10]:
            length = repeat['end'] - repeat['start']
            f.write(f"| {repeat['chrom']} | {repeat['start']} | {repeat['end']} | {length} | {repeat['score']} | {repeat['name']} |\n")
        
        # f.write("\n## Top 10 Closest Repeats to IgH Genes\n\n")
        # closest_repeats = [r for r in repeats if r.get('distance_to_igh', float('inf')) != float('inf')]
        # sorted_by_distance = sorted(closest_repeats, key=lambda r: r.get('distance_to_igh', float('inf')))
        
        # if sorted_by_distance:
        #     f.write("| Chromosome | Start | End | Length | Copy Number | Distance to IgH (bp) | Name |\n")
        #     f.write("|------------|-------|-----|--------|-------------|---------------------|------|\n")
        #     for repeat in sorted_by_distance[:10]:
        #         length = repeat['end'] - repeat['start']
        #         f.write(f"| {repeat['chrom']} | {repeat['start']} | {repeat['end']} | {length} | "
        #                f"{repeat['score']} | {repeat['distance_to_igh']} | {repeat['name']} |\n")
        # else:
        #     f.write("No repeats with distance information available.\n")

def main():
    parser = argparse.ArgumentParser(description="Generate a report on identified repeat sequences")
    parser.add_argument('--repeats', required=True, help="BED file with filtered repeats")
    parser.add_argument('--igh_genes', required=True, help="BED file with IgH gene locations")
    parser.add_argument('--output_dir', required=True, help="Output directory for report and plots")
    
    args = parser.parse_args()
    
    # Check if input files exist
    for file_path in [args.repeats, args.igh_genes]:
        if not os.path.exists(file_path):
            print(f"Error: Input file {file_path} does not exist", file=sys.stderr)
            sys.exit(1)
    
    # Create output directory if it doesn't exist
    if args.output_dir and not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
    
    # Load data
    print("Loading repeat data...")
    repeats = load_bed_file(args.repeats)
    print(f"Loaded {len(repeats)} repeats")
    
    print("Loading IgH gene data...")
    igh_genes = load_bed_file(args.igh_genes)
    print(f"Loaded {len(igh_genes)} IgH genes")
    
    # Calculate distances to IgH genes
    print("Calculating distances to IgH genes...")
    calculate_distance_to_igh(repeats, igh_genes)
    
    # Generate statistics
    print("Generating statistics...")
    stats = generate_statistics(repeats)
    
    # Create visualizations
    print("Creating plots...")
    plot_statistics(stats, args.output_dir)
    
    # Write text report
    print("Writing report...")
    write_report(stats, repeats, f"{args.output_dir}/report.md")
    
    print(f"Report and visualizations generated with prefix {args.output_dir}")

if __name__ == "__main__":
    main() 