#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Generate a comprehensive report for the sgRNA design pipeline.

This script generates HTML and text reports summarizing the results of the sgRNA
design pipeline, including statistics about the repeats, sgRNAs, and their distribution.
It also creates visualizations to help interpret the results.

Author: Ming Wang
Version: 1.1
Date: 2025-03-26
"""

import os
import sys
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
from collections import Counter
import json

def parse_bed6plus(filename, header=None):
    """Parse a BED6+ file into a pandas DataFrame."""
    if not os.path.exists(filename):
        print(f"Error: File not found: {filename}")
        return None
    
    try:
        if header is None:
            df = pd.read_csv(filename, sep='\t', header=None)
        else:
            df = pd.read_csv(filename, sep='\t', names=header)
        return df
    except Exception as e:
        print(f"Error parsing file {filename}: {e}")
        return None

def parse_sgrna_file(filename):
    """Parse the sgRNA output file with specific column names."""
    header = ['chrom', 'start', 'end', 'repeat_name', 'copy_number', 'strand', 
              'consensus_seq', 'direction', 'sgrna_id', 'sgrna_seq']
    return parse_bed6plus(filename, header=header)

def parse_repeats_file(filename):
    """Parse the repeats BED file."""
    header = ['chrom', 'start', 'end', 'repeat_name', 'copy_number', 'strand', 
              'consensus_seq']
    return parse_bed6plus(filename, header=header)

def parse_igh_genes(filename):
    """Parse the IgH genes BED file."""
    header = ['chrom', 'start', 'end', 'name', 'score', 'strand']
    return parse_bed6plus(filename, header=header)

def calculate_statistics(sgrna_df, repeats_df):
    """Calculate statistics about the sgRNAs and repeats."""
    stats = {}
    
    # Basic counts
    stats['total_sgrna'] = len(sgrna_df) if sgrna_df is not None else 0
    stats['total_repeats'] = len(repeats_df) if repeats_df is not None else 0
    
    if sgrna_df is not None and len(sgrna_df) > 0:
        # sgRNA direction distribution
        direction_counts = Counter(sgrna_df['direction'])
        stats['direction_fwd'] = direction_counts.get('fwd', 0)
        stats['direction_rev'] = direction_counts.get('rev', 0)
        
        # sgRNA length distribution
        stats['sgrna_length'] = int(sgrna_df['sgrna_seq'].str.len().mean())
        
        # GC content calculation
        sgrna_df['gc_content'] = sgrna_df['sgrna_seq'].apply(
            lambda x: (x.count('G') + x.count('C')) / len(x) * 100
        )
        stats['gc_content_mean'] = round(sgrna_df['gc_content'].mean(), 2)
        stats['gc_content_min'] = round(sgrna_df['gc_content'].min(), 2)
        stats['gc_content_max'] = round(sgrna_df['gc_content'].max(), 2)
        
        # Copy number statistics
        stats['copy_number_mean'] = round(sgrna_df['copy_number'].mean(), 2)
        stats['copy_number_min'] = int(sgrna_df['copy_number'].min())
        stats['copy_number_max'] = int(sgrna_df['copy_number'].max())
    
    if repeats_df is not None and len(repeats_df) > 0:
        # Repeat length statistics
        repeats_df['repeat_length'] = repeats_df['end'] - repeats_df['start']
        stats['repeat_length_mean'] = round(repeats_df['repeat_length'].mean(), 2)
        stats['repeat_length_min'] = int(repeats_df['repeat_length'].min())
        stats['repeat_length_max'] = int(repeats_df['repeat_length'].max())
    
    return stats

def create_visualizations(sgrna_df, repeats_df, output_dir):
    """Create visualizations of the sgRNA and repeat data."""
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    if sgrna_df is not None and len(sgrna_df) > 0:
        # Figure 1: sgRNA direction distribution
        plt.figure(figsize=(8, 6))
        direction_counts = Counter(sgrna_df['direction'])
        plt.bar(['Forward', 'Reverse'], 
                [direction_counts.get('fwd', 0), direction_counts.get('rev', 0)])
        plt.title('sgRNA Direction Distribution')
        plt.ylabel('Count')
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, 'sgrna_direction_dist.png'), dpi=300)
        plt.close()
        
        # Figure 2: GC content distribution
        plt.figure(figsize=(10, 6))
        plt.hist(sgrna_df['gc_content'], bins=20, alpha=0.7, color='skyblue')
        plt.axvline(sgrna_df['gc_content'].mean(), color='red', linestyle='dashed', 
                    linewidth=1, label=f'Mean: {sgrna_df["gc_content"].mean():.2f}%')
        plt.title('sgRNA GC Content Distribution')
        plt.xlabel('GC Content (%)')
        plt.ylabel('Frequency')
        plt.legend()
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, 'sgrna_gc_content.png'), dpi=300)
        plt.close()
        
        # Figure 3: Copy number distribution
        plt.figure(figsize=(10, 6))
        plt.hist(sgrna_df['copy_number'], bins=20, alpha=0.7, color='lightgreen')
        plt.title('Repeat Copy Number Distribution')
        plt.xlabel('Copy Number')
        plt.ylabel('Frequency')
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, 'repeat_copy_number.png'), dpi=300)
        plt.close()
        
        # Figure 4: sgRNA position distribution
        if 'start' in sgrna_df.columns:
            plt.figure(figsize=(12, 6))
            plt.hist(sgrna_df['start'], bins=50, alpha=0.7, color='salmon')
            plt.title('sgRNA Genomic Position Distribution')
            plt.xlabel('Genomic Position')
            plt.ylabel('Frequency')
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, 'sgrna_position_dist.png'), dpi=300)
            plt.close()
    
    if repeats_df is not None and len(repeats_df) > 0:
        # Figure 5: Repeat length distribution
        repeats_df['repeat_length'] = repeats_df['end'] - repeats_df['start']
        plt.figure(figsize=(10, 6))
        plt.hist(repeats_df['repeat_length'], bins=20, alpha=0.7, color='plum')
        plt.title('Repeat Length Distribution')
        plt.xlabel('Length (bp)')
        plt.ylabel('Frequency')
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, 'repeat_length_dist.png'), dpi=300)
        plt.close()

def generate_html_report(stats, output_dir):
    """Generate an HTML report with the statistics and visualizations."""
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    
    html = f"""<!DOCTYPE html>
<html>
<head>
    <title>sgRNA Design Report</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 20px; line-height: 1.6; }}
        h1, h2, h3 {{ color: #333366; }}
        .container {{ max-width: 1200px; margin: 0 auto; }}
        .stats-container {{ display: flex; flex-wrap: wrap; gap: 20px; margin-bottom: 30px; }}
        .stat-box {{ background: #f5f5f5; border-radius: 5px; padding: 15px; flex: 1; min-width: 200px; }}
        .image-container {{ display: flex; flex-wrap: wrap; gap: 20px; justify-content: center; }}
        .image-box {{ margin-bottom: 20px; text-align: center; }}
        .image-box img {{ max-width: 100%; height: auto; border: 1px solid #ddd; }}
        table {{ width: 100%; border-collapse: collapse; margin-bottom: 20px; }}
        th, td {{ padding: 8px; text-align: left; border-bottom: 1px solid #ddd; }}
        th {{ background-color: #f2f2f2; }}
    </style>
</head>
<body>
    <div class="container">
        <h1>sgRNA Design Pipeline Report</h1>
        <p>Report generated on {timestamp}</p>
        
        <h2>Summary Statistics</h2>
        <div class="stats-container">
            <div class="stat-box">
                <h3>sgRNA Statistics</h3>
                <p>Total sgRNAs: <strong>{stats.get('total_sgrna', 0)}</strong></p>
                <p>Forward strand: <strong>{stats.get('direction_fwd', 0)}</strong></p>
                <p>Reverse strand: <strong>{stats.get('direction_rev', 0)}</strong></p>
                <p>Average length: <strong>{stats.get('sgrna_length', 'N/A')}</strong> bp</p>
            </div>
            
            <div class="stat-box">
                <h3>GC Content</h3>
                <p>Mean: <strong>{stats.get('gc_content_mean', 'N/A')}%</strong></p>
                <p>Min: <strong>{stats.get('gc_content_min', 'N/A')}%</strong></p>
                <p>Max: <strong>{stats.get('gc_content_max', 'N/A')}%</strong></p>
            </div>
            
            <div class="stat-box">
                <h3>Repeat Statistics</h3>
                <p>Total repeats: <strong>{stats.get('total_repeats', 0)}</strong></p>
                <p>Mean length: <strong>{stats.get('repeat_length_mean', 'N/A')}</strong> bp</p>
                <p>Copy number range: <strong>{stats.get('copy_number_min', 'N/A')} - {stats.get('copy_number_max', 'N/A')}</strong></p>
            </div>
        </div>
        
        <h2>Visualizations</h2>
        <div class="image-container">
            <div class="image-box">
                <img src="sgrna_direction_dist.png" alt="sgRNA Direction Distribution">
                <p>Figure 1: Distribution of sgRNA directions</p>
            </div>
            
            <div class="image-box">
                <img src="sgrna_gc_content.png" alt="sgRNA GC Content Distribution">
                <p>Figure 2: Distribution of GC content in sgRNAs</p>
            </div>
            
            <div class="image-box">
                <img src="repeat_copy_number.png" alt="Repeat Copy Number Distribution">
                <p>Figure 3: Distribution of repeat copy numbers</p>
            </div>
            
            <div class="image-box">
                <img src="sgrna_position_dist.png" alt="sgRNA Position Distribution">
                <p>Figure 4: Genomic position distribution of sgRNAs</p>
            </div>
            
            <div class="image-box">
                <img src="repeat_length_dist.png" alt="Repeat Length Distribution">
                <p>Figure 5: Distribution of repeat lengths</p>
            </div>
        </div>
        
        <h2>Methods</h2>
        <p>
            This report was generated using the sgRNA design pipeline, which identifies potential 
            CRISPR-Cas9 sgRNA targets in regions with tandem repeats. The pipeline:
        </p>
        <ol>
            <li>Identifies tandem repeats in the target region</li>
            <li>Extracts candidate sgRNAs (23bp, ending with NGG PAM)</li>
            <li>Filters sgRNAs for specificity and on-target scoring</li>
        </ol>
        
        <h2>Recommendations</h2>
        <p>
            Based on the analysis, we recommend selecting sgRNAs with the following characteristics:
        </p>
        <ul>
            <li>GC content between 40-60%</li>
            <li>Located in regions with copy number ≥ 10</li>
            <li>Minimal predicted off-target effects</li>
        </ul>
        
        <p>
            <em>Note: This report is automatically generated. Please review the results 
            carefully before selecting sgRNAs for experimental validation.</em>
        </p>
    </div>
</body>
</html>
"""
    
    # Write HTML report
    with open(os.path.join(output_dir, 'report.html'), 'w') as f:
        f.write(html)
    
    # Save statistics as JSON
    with open(os.path.join(output_dir, 'statistics.json'), 'w') as f:
        json.dump(stats, f, indent=2)
    
    # Generate text report
    with open(os.path.join(output_dir, 'report.txt'), 'w') as f:
        f.write("sgRNA Design Pipeline Report\n")
        f.write("===========================\n\n")
        f.write(f"Report generated on {timestamp}\n\n")
        
        f.write("Summary Statistics\n")
        f.write("------------------\n")
        f.write(f"Total sgRNAs: {stats.get('total_sgrna', 0)}\n")
        f.write(f"Forward strand: {stats.get('direction_fwd', 0)}\n")
        f.write(f"Reverse strand: {stats.get('direction_rev', 0)}\n")
        f.write(f"Average length: {stats.get('sgrna_length', 'N/A')} bp\n\n")
        
        f.write(f"Mean GC content: {stats.get('gc_content_mean', 'N/A')}%\n")
        f.write(f"GC content range: {stats.get('gc_content_min', 'N/A')}% - {stats.get('gc_content_max', 'N/A')}%\n\n")
        
        f.write(f"Total repeats: {stats.get('total_repeats', 0)}\n")
        f.write(f"Mean repeat length: {stats.get('repeat_length_mean', 'N/A')} bp\n")
        f.write(f"Copy number range: {stats.get('copy_number_min', 'N/A')} - {stats.get('copy_number_max', 'N/A')}\n\n")
        
        f.write("Visualizations are available in the HTML report.\n")

def main():
    parser = argparse.ArgumentParser(description='Generate a comprehensive report for the sgRNA design pipeline')
    parser.add_argument('--repeats', required=True, help='Path to the filtered repeats BED file')
    parser.add_argument('--sgrna', required=True, help='Path to the on-target sgRNA file')
    parser.add_argument('--igh_genes', required=True, help='Path to the IgH genes BED file')
    parser.add_argument('--output', required=True, help='Directory to save the report')
    
    args = parser.parse_args()
    
    # Create output directory if it doesn't exist
    if not os.path.exists(args.output):
        os.makedirs(args.output)
    
    # Parse input files
    sgrna_df = parse_sgrna_file(args.sgrna)
    repeats_df = parse_repeats_file(args.repeats)
    igh_genes_df = parse_igh_genes(args.igh_genes)
    
    if sgrna_df is None or repeats_df is None:
        print("Error: Failed to parse input files. Check that they exist and have the correct format.")
        sys.exit(1)
    
    # Calculate statistics
    stats = calculate_statistics(sgrna_df, repeats_df)
    
    # Create visualizations
    create_visualizations(sgrna_df, repeats_df, args.output)
    
    # Generate report
    generate_html_report(stats, args.output)
    
    print(f"Report generated successfully. View it at {os.path.join(args.output, 'report.html')}")

if __name__ == "__main__":
    main() 