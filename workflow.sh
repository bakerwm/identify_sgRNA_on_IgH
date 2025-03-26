#!/bin/bash

# Date: 2025-03-25
# Author: Ming Wang
# Email: wangm08@hotmail.com
# Version: 1.0

# Description:
#
# Design sgRNAs (SpCas9, NGG) targeting IgH cluster and its upstream region, 2 Mbp.
#
# How-To:
# 1. Identify IgH cluster region and its upstream region in human (GRCh38, Ensembl release 102)
#    - chr14:105,583,730-106,879,812 (-) (GRCh38/hg38)
#    - IgH cluster and upstream 2 Mbp: chr14:105,583,730-108,879,812 (-) (GRCh38/hg38)
# 2. Identify tandom repeat sequences in the IgH cluster and its upstream region
# 3. Design sgRNAs targeting the random repeat sequences, on consensus sequence
# 4. Remove off-target sgRNAs (bowtie2 alignment)
# 5. Generate report
#
# Input:
#    - Human reference genome (GRCh38/hg38)
#    - Human annotation (GRCh38/hg38)
#    - IgH cluster and its upstream region (GRCh38/hg38)
#
# Output:
#    - results/sgrna/sgrna_candidates.on_target.txt
#
#    columns:
#      1. chromosome
#      2. start
#      3. end
#      4. tandom_repeat_name
#      5. copy_number
#      6. strand
#      7. consensus_sequence
#      8. direction_of_sgRNA:fwd/rev
#      9. sgRNA_id
#      10. sgRNA_sequence
#
# Gettig started:
# $ bash workflow.sh

################################################################################
# BEGIN: GLOBAL VARIABLES                                                      #
export GENOME_BUILD="GRCh38"
export RELEASE=102
export FLANKING_DISTANCE=2000000  # 2 Mbp around IgH genes
export IgH_REGION="14"  # Human IgH is on chromosome 14
export MIN_COPY=10      # Minimum copy number of repeat sequences
export MIN_LENGTH=13    # Minimum length of repeat sequences
export BOWTIE2_IDX="/data/biodata/genome_db/GRCh38/Ensembl/bowtie2_index/GRCh38"
# export BOWTIE2_IDX=""   # will build bowtie2 index
# END: GLOBAL VARIABLES                                                        #
################################################################################

set -e  # Exit on error

export N_CPU=12

################################################################################
# Configuration                                                                #
GENOME_DIR="genome"
OUTPUT_DIR="results"
SCRIPTS_DIR="scripts"
GENOME_FASTA="${GENOME_DIR}/${GENOME_BUILD}.fa"
GENOME_FAI="${GENOME_DIR}/${GENOME_BUILD}.fa.fai"
GENOME_GTF="${GENOME_DIR}/${GENOME_BUILD}.gtf"
# Create necessary directories
mkdir -p $GENOME_DIR $OUTPUT_DIR $SCRIPTS_DIR
################################################################################
# Step 1: Download human reference genome if not already present
echo "Step 1: Preparing human reference genome, annotation-GTF..."
bash ${SCRIPTS_DIR}/01.download_genome.sh ${GENOME_BUILD} ${GENOME_DIR} ${RELEASE}

# Step 2: Identify IgH gene locations
echo "Step 2: Identifying IgH gene locations..."
igh_regions_bed="${OUTPUT_DIR}/igh_regions.bed"
if [ ! -f ${igh_regions_bed} ]; then
    python ${SCRIPTS_DIR}/02.identify_igh_genes.py \
        --gtf ${GENOME_GTF} \
        --output ${igh_regions_bed}
fi

# Step 3: Extract sequences from regions around IgH genes
echo "Step 3: Extracting sequences from IgH regions and flanking regions..."
# human, upstream, right flanking region: 2 Mbp
flanking_bed="${OUTPUT_DIR}/igh_regions_flank_right_2Mb.bed"
flanking_fa="${OUTPUT_DIR}/igh_regions_flank_right_2Mb.fa"
if [ ! -f ${flanking_bed} ]; then
    bedtools slop -l 0 -r ${FLANKING_DISTANCE} -g ${GENOME_FAI} \
        -i ${OUTPUT_DIR}/igh_regions.bed  \
        > ${flanking_bed}
    bedtools getfasta -fi ${GENOME_FASTA} -bed ${flanking_bed} -fo ${flanking_fa}
    # alternative: use samtools faidx <genome.fa> <region> > <output.fa>
fi

# Step 4: Find repeat sequences in the extracted regions
# Using only TRF and RepeatMasker (skipping k-mer analysis and MISA)
echo "Step 4: Finding repeat sequences..."
raw_repeats_bed="${OUTPUT_DIR}/repeats/all_repeats.raw.bed"
if [ ! -f ${raw_repeats_bed} ]; then
    bash ${SCRIPTS_DIR}/03.find_repeats.sh \
        ${flanking_fa} \
        ${raw_repeats_bed} \
        ${MIN_LENGTH} \
        ${MIN_COPY}
fi

# Step 5: Filter repeats based on criteria (length > 13bp, copy number > 10)
echo "Step 5: Filtering repeats based on criteria..."
filtered_repeats_bed="${OUTPUT_DIR}/repeats/all_repeats.filtered.bed"
if [ ! -f ${filtered_repeats_bed} ]; then
    python ${SCRIPTS_DIR}/05.filter_repeats.py \
        --input ${raw_repeats_bed} \
        --min_length ${MIN_LENGTH} \
        --min_copy_number ${MIN_COPY} \
        --output ${filtered_repeats_bed}
fi

# Step 6: Extract sgRNA candidate loci (23bp, ending with GG)
echo "Step 6: Extracting sgRNA candidates..."
sgrna_txt="${OUTPUT_DIR}/sgRNA/sgrna_raw.txt"
sgrna_fa="${sgrna_txt%.txt}.fa"
if [ ! -f ${sgrna_fa} ]; then
    python ${SCRIPTS_DIR}/06.identify_sgrna.py \
        --input ${filtered_repeats_bed} \
        --output ${sgrna_txt} --size 23 --suffix GG
fi

# Step 7: Map sgRNA sequences to reference genome
# Args: <output_dir> <sgrna.txt> <bowtie2_idx> [n_cpu] [target_regions_bed]
echo "Step 7: remove off-target sgRNAs..."
bash ${SCRIPTS_DIR}/07.remove_off_targets.sh \
    ${OUTPUT_DIR}/sgRNA \
    ${sgrna_txt} \
    ${BOWTIE2_IDX} \
    ${N_CPU} \
    ${OUTPUT_DIR}/igh_regions.bed

# Step 8: Generate summary and report
echo "Step 8: Generating report..."
# python ${SCRIPTS_DIR}/08.generate_report.py \
#     --repeats ${filtered_repeats_bed} \
#     --igh_genes ${OUTPUT_DIR}/igh_regions.bed \
#     --output ${OUTPUT_DIR}/report

echo "Pipeline completed successfully!"
echo "Results are saved in: ${OUTPUT_DIR}/sgRNA/sgrna_raw.on_target.txt" 
