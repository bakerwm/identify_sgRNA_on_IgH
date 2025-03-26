# Identify sgRNAs tageting IgH cluster in human genome

## Criteria
+ sgRNA is 23mer and ending with GG
+ sgRNA targets should be in tandom repeat (copy number > 10)
+ sgRNA targeting IgH and its **upstream** 2 Mbp region
+ sgRNA MUST not targeting other regions (off-target)

## Pipeline Overview

This pipeline performs the following steps:

1. Download and prepare the human reference genome (fasta, gtf) (GRCh38)
2. Identify IgH gene locations in the human genome (gene_name: IGH, Igh)
3. Detect tandom repeats using Tandem Repeats Finder (TRF)
4. Filter repeats based on criteria (length > 13bp, copy number > 10)
5. Extract candidate sgRNA loci, 23 bp in length and ending with GG (SpCas9)
6. Remove off-targets by genome-wide alignment

## Requirements

The pipeline requires the following tools and packages:

- Python 3.6+
  - pybedtools
  - biopython
  - NumPy
- bowtie2
- samtools
- bedtools
- trf (Tandem Repeats Finder)

## Usage

To run the complete pipeline:

```bash
bash workflow.sh
```

### Individual Components

The pipeline is composed of multiple scripts that can be run independently:

- `01.download_genome.sh`: Download and prepare the human reference genome
- `02.identify_igh_genes.py`: Identify IgH gene locations in the human genome
- `03.find_repeats.sh`: Find repeat sequences using various tools
- `04.combind_repeats.py`: Combine repeats generated by all tools
- `05.filter_repeats.py`: Filter repeats based on criteria
- `06.identify_sgrna.py`: Identify all possible sgRNAs within consensus sequence
- `07.remove_off_targets.sh`: Remove off-targets sgRNAs, that could map to other regions
- `08.generate_report.py`: Generate a summary report of the results [optional]

## Output

The pipeline generates the following outputs in the `results` directory:

- `igh_genes.bed`: BED file of IgH gene locations
- `igh_regions.bed`: BED file of IgH cluster location
- `igh_regions_flank_left_2Mb.bed`: BED file of IgH cluster and downstream 2 Mbp
- `igh_regions_flank_left_2Mb.fa`: FASTA file of IgH cluster and dwonstream 2 Mbp
- `repeats/`: Directory containing raw outputs from repeat-finding tools
- `repeats/all_repeats.filtered.bed`: BED6+1 file of tandom repeats
- `sgRNA/`: Directory containing sgRNA sequences and files
- `sgRNA/sgrna_raw.on_target.txt`: BED6+4 file of on target sgRNAs

**sgRNA on target**

file: `results/sgRNA/sgrna_raw.on_target.txt`
columns:
1. chromosome  
2. start  
3. end  
4. tandom_repeat_name  
5. copy_number  
6. strand  
7. consensus_sequence  
8. direction_of_sgRNA:fwd/rev  
9. sgRNA_id  
10. sgRNA_sequence

## Customization

You can modify the following parameters in the workflow script:

- `GENOME_BUILD`: GRCh38
- `RELEASE`: Ensembl release, eg: 102
- `BOWTIE2_IDX`: Bowtie2 index, if not specified, will build in `${GENOME_DIR}`
- `FLANKING_DISTANCE`: Flanking distance around IgH genes (default: 2 Mbp)
- `MIN_LENGTH`: Minimum length of repeats (default: 13 bp)
- `MIN_COPY_NUMBER`: Minimum copy number (default: 10)
- `N_CPU`: Number of CPUs to run program (default: 12)
- `GENOME_DIR`: Directory for genome files (default: genome)
- `OUTPUT_DIR`: Directory for output files (default: results)

## License

This software is provided under the MIT License.

## Citation

If you use this pipeline in your research, please cite:

[Citation information] 