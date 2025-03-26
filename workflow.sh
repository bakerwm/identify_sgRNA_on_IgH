#!/bin/bash

# Date: 2025-03-25
# Author: Ming Wang
# Email: wangm08@hotmail.com
# Version: 1.0

# Description:
#
# Design sgRNAs (SpCas9, NGG) targeting specific region on genome
#
# How-To:
# 1. Identify target region on genome (BED)
# 2. Identify tandom repeat sequences in the target region
# 3. Design sgRNAs targeting the random repeat sequences, on consensus sequence
# 4. Remove off-target sgRNAs (bowtie2 alignment)
# 5. (to-do) Filter sgRNA for better performance
# 6. Generate report
#
# Input:
#    1. Genome_build (eg: GRCm38, GRCh38)
#    2. Target region (BED3)
#    3. Flanking distance (int, bp)
#
# Output:
#    - results/sgrna/sgrna_raw.on_target.txt
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
#
# Gettig started:
# $ bash workflow.sh config.sh

set -e  # Exit on error

if [[ $# -ne 1 ]]; then
    echo "Usage: $0 <config_file>"
    exit 1
fi
config_file=$1

if [[ ! -f ${config_file} ]] ; then
    echo "Error: config file not exists: ${config_file}"
    exit 1
fi
source ${config_file} # load configuration

# ################################################################################
# # BEGIN: GLOBAL VARIABLES                                                      #
# export N_CPU=12
# export GENOME_BUILD="GRCm38"
# export RELEASE=102
# export TARGET_BED="igh_regions.bed"  # Required, BED3 format
# export FLANKING_LEFT=2000000  # 2 Mbp around IgH genes
# export FLANKING_RIGHT=0  # 2 Mbp around IgH genes
# export MIN_COPY=10      # Minimum copy number of repeat sequences
# export MIN_LENGTH=13    # Minimum length of repeat sequences
# export BOWTIE2_IDX="/data/biodata/genome_db/GRCm38/Ensembl/bowtie2_index/GRCm38"   # will build bowtie2 index
# # END: GLOBAL VARIABLES                                                        #
# ################################################################################

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

# Step 0: Check required file
if [[ ! -f ${TARGET_BED} ]]; then
    echo "Error: Target bed file not exists: ${TARGET_BED}"
    exit 1
fi

# Step 1: Download reference genome if not already present
echo "Step 1: Preparing reference genome, build bowtie2 index..."
bash ${SCRIPTS_DIR}/01.download_genome.sh ${GENOME_BUILD} ${GENOME_DIR} ${RELEASE}

# Step 2: Prepare target bed file
echo "Step 2: Preparing target fasta file..."
FNAME=$(basename ${TARGET_BED} .bed)
target_extended_bed="${OUTPUT_DIR}/${FNAME}.extended.bed"
target_extended_fa="${OUTPUT_DIR}/${FNAME}.extended.fa"
if [ ! -f ${target_extended_bed} ]; then
    bash ${SCRIPTS_DIR}/02.prepare_target_bed.sh \
        ${target_extended_bed} \
        ${TARGET_BED} \
        ${GENOME_FASTA} \
        ${FLANKING_LEFT} \
        ${FLANKING_RIGHT}
fi

# Step 3: Find repeat sequences in the target region
echo "Step 3: Finding repeat sequences in the target region..."
repeats_raw_bed="${OUTPUT_DIR}/repeats/all_repeats.raw.bed"
if [ ! -f ${repeats_raw_bed} ]; then
    bash ${SCRIPTS_DIR}/03.find_repeats.sh \
        ${target_extended_fa} \
        ${repeats_raw_bed} \
        ${MIN_LENGTH} \
        ${MIN_COPY}
fi

# Step 4: Filter repeats based on criteria (length > 13bp, copy number > 10)
echo "Step 4: Filtering repeats based on criteria..."
filtered_repeats_bed="${OUTPUT_DIR}/repeats/all_repeats.filtered.bed"
if [ ! -f ${filtered_repeats_bed} ]; then
    python ${SCRIPTS_DIR}/05.filter_repeats.py \
        --input ${repeats_raw_bed} \
        --min_length ${MIN_LENGTH} \
        --min_copy_number ${MIN_COPY} \
        --output ${filtered_repeats_bed}
fi

# Step 5: Extract sgRNA candidate loci (23bp, ending with GG)
echo "Step 5: Extracting sgRNA candidates..."
sgrna_txt="${OUTPUT_DIR}/sgRNA/sgrna_raw.txt"
sgrna_fa="${sgrna_txt%.txt}.fa"
if [ ! -f ${sgrna_fa} ]; then
    python ${SCRIPTS_DIR}/06.identify_sgrna.py \
        --input ${filtered_repeats_bed} \
        --output ${sgrna_txt} --size 23 --suffix GG
fi

# Step 6: Remove off-target sgRNAs
# Args: <output_dir> <sgrna.txt> <bowtie2_idx> [n_cpu] [target_regions_bed]
echo "Step 6: remove off-target sgRNAs..."
bash ${SCRIPTS_DIR}/07.remove_off_targets.sh \
    ${OUTPUT_DIR}/sgRNA \
    ${sgrna_txt} \
    ${BOWTIE2_IDX} \
    ${N_CPU} \
    ${target_extended_bed}

# Step 7: Filter sgRNAs for better performance
# TODO:

# # Step 8: Generate summary and report
# echo "Step 8: Generating report..."
# python ${SCRIPTS_DIR}/08.generate_report.py \
#     --repeats ${filtered_repeats_bed} \
#     --igh_genes ${OUTPUT_DIR}/igh_regions.bed \
#     --output ${OUTPUT_DIR}/report

# Create symlink to the final sgRNA file
sgrna_txt="${OUTPUT_DIR}/sgrna_raw.on_target.txt"
[[ ! -f ${sgrna_txt} ]] && ln -s sgRNA/sgrna_raw.on_target.txt ${OUTPUT_DIR}/

echo "Pipeline completed successfully!"
echo "Results are saved in: ${OUTPUT_DIR}/sgrna_raw.on_target.txt" 
