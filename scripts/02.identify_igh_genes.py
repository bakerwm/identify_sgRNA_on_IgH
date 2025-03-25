#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = "Ming Wang"
__email__ = "wangm08@hotmail.com"
__version__ = "1.0"
__date__ = "2025-03-25"

"""
Script to identify IgH genes in the human genome by searching gene_name for:
IGH*, Igh*, immunoglobulin, heavy

IgH genes are typically located on chromosome 12.
with the following coordinates:
chr14:105,583,730-106,879,812 (-) (GRCh38/hg38)

Output:
  - BED file with IgH genes
"""

import argparse
import pybedtools

def identify_igh_genes(gtf_file, output_file, chr_list=['14']):
    """
    Parse GTF file to identify IgH genes.
    IgH gene segments typically have IGHV, IGHD, IGHJ, or IGHC prefixes.
    """
    print("Searching for IgH genes...")
    igh_genes = []
    
    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
                
            fields = line.strip().split('\t')
            if len(fields) < 9 or fields[2] != 'gene':
                continue
            
            # Check if the chromosome is in the list
            if fields[0] not in chr_list:
                continue
                
            # Extract gene information from the attributes field
            attr_dict = {}
            for attr in fields[8].split(';'):
                attr = attr.strip()
                if not attr:
                    continue
                key, value = attr.split(' ', 1)
                attr_dict[key] = value.strip('"')
            
            gene_name = attr_dict.get('gene_name', '')
            gene_id = attr_dict.get('gene_id', '')
            
            # Check if this is an IgH gene
            if (gene_name.startswith('IGH') or 
                gene_name.startswith('Igh') or  # Added to catch human gene naming
                (gene_id.startswith('ENSMUSG') and 
                ('immunoglobulin' in attr_dict.get('description', '').lower() and 
                 'heavy' in attr_dict.get('description', '').lower()))):
                
                chrom = fields[0]
                start = int(fields[3]) - 1  # Convert to 0-based for BED format
                end = int(fields[4])
                
                igh_genes.append({
                    'chrom': chrom,
                    'start': start,
                    'end': end,
                    'name': gene_name if gene_name else gene_id,
                    'score': 0,
                    'strand': fields[6]
                })
    
    print(f"Found {len(igh_genes)} IgH genes")
    
    # Write results to BED file
    with open(output_file, 'w') as out:
        for gene in igh_genes:
            out.write(f"{gene['chrom']}\t{gene['start']}\t{gene['end']}\t{gene['name']}\t{gene['score']}\t{gene['strand']}\n")
    
    # If we didn't find any IgH genes through annotation, use known coordinates
    if len(igh_genes) == 0:
        print("No IgH genes found in annotation, using known coordinates")
        # human IgH locus is located on chromosome 12, approximately at 113,258,768-116,009,954 in mm10
        with open(output_file, 'w') as out:
            out.write("chr12\t113258768\t116009954\tIgH_locus\t0\t+\n")

# Merge the regions into a single cluster region
def merge_regions(regions_file, output_file, max_distance = 1000000):
    """Merge overlapping genomic regions
    Using command: bedtools merge -i <regions_file> -d <max_distance> > <output_file>

    Using pybedtools instead.
    import pybedtools
    pybedtools.BedTool(regions_file).merge(d=max_distance).saveas(output_file)
    
    Args:
        regions_file: BED file with genomic regions
        output_file: BED file to write merged regions
        max_distance: Maximum distance between regions to merge
    """
    # sort the bed file by chr and start
    sorted_bed = pybedtools.BedTool(regions_file).sort()
    # Using pybedtools
    merged_bed = sorted_bed.merge(d=max_distance, c='4,5,6', o='first,first,first')
    # save to file
    merged_bed.saveas(output_file)   

def main():
    parser = argparse.ArgumentParser(description="Identify IgH genes in human genome")
    parser.add_argument('--gtf', required=True, help="Path to the genome GTF file")
    parser.add_argument('--output', required=True, help="Output BED file with IgH gene locations")
    parser.add_argument('--chr_list', required=False, default=['14'], help="List of chromosomes to search for IgH gene, default: 14")
    args = parser.parse_args()

    # Identify IgH genes
    igh_genes_bed = args.output.replace('_regions.bed', '_genes.bed')
    identify_igh_genes(args.gtf, igh_genes_bed, args.chr_list)

    # Merge regions
    merge_regions(igh_genes_bed, args.output, 1000000)
    
    print(f"IgH gene locations written to {args.output}")

if __name__ == "__main__":
    main() 
