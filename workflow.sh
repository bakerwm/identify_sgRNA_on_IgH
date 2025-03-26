#!/bin/bash

# Date: 2025-03-25
# Author: Ming Wang
# Email: wangm08@hotmail.com
# Version: 1.1

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

# Function to handle errors
handle_error() {
    local line=$1
    local cmd=$2
    echo "Error at line $line: Command '$cmd' failed with exit code $?"
    exit 1
}

# Set error trap
trap 'handle_error ${LINENO} "$BASH_COMMAND"' ERR

# Function to log progress with timestamps
log_progress() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

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
LOG_FILE="${OUTPUT_DIR}/workflow_$(date '+%Y%m%d_%H%M%S').log"

# Create necessary directories
mkdir -p $GENOME_DIR $OUTPUT_DIR $SCRIPTS_DIR

# Redirect all output to both terminal and log file
exec > >(tee -a "${LOG_FILE}") 2>&1

log_progress "Starting sgRNA design pipeline for ${GENOME_BUILD}"
log_progress "Configuration: ${config_file}"
################################################################################

# Step 0: Check required file
log_progress "Step 0: Checking required files and dependencies"
if [[ ! -f ${TARGET_BED} ]]; then
    echo "Error: Target bed file not exists: ${TARGET_BED}"
    exit 1
fi

# Validate required tools
for cmd in bowtie2 samtools bedtools trf python; do
    if ! command -v $cmd &> /dev/null; then
        echo "Error: Required tool '$cmd' not found in PATH"
        echo "Please install required dependencies or activate conda environment"
        exit 1
    fi
done

# Step 1: Download reference genome if not already present
log_progress "Step 1: Preparing reference genome, build bowtie2 index..."
bash ${SCRIPTS_DIR}/01.download_genome.sh ${GENOME_BUILD} ${GENOME_DIR} ${RELEASE}

# Step 2: Prepare target bed file
log_progress "Step 2: Preparing target fasta file..."
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
else
    log_progress "  > Skipping: File already exists: ${target_extended_bed}"
fi

# Step 3: Find repeat sequences in the target region
log_progress "Step 3: Finding repeat sequences in the target region..."
repeats_raw_bed="${OUTPUT_DIR}/repeats/all_repeats.raw.bed"
if [ ! -f ${repeats_raw_bed} ]; then
    bash ${SCRIPTS_DIR}/03.find_repeats.sh \
        ${target_extended_fa} \
        ${repeats_raw_bed} \
        ${MIN_LENGTH} \
        ${MIN_COPY}
else
    log_progress "  > Skipping: File already exists: ${repeats_raw_bed}"
fi

# Step 4: Filter repeats based on criteria (length > 13bp, copy number > 10)
log_progress "Step 4: Filtering repeats based on criteria..."
filtered_repeats_bed="${OUTPUT_DIR}/repeats/all_repeats.filtered.bed"
if [ ! -f ${filtered_repeats_bed} ]; then
    python ${SCRIPTS_DIR}/05.filter_repeats.py \
        --input ${repeats_raw_bed} \
        --min_length ${MIN_LENGTH} \
        --min_copy_number ${MIN_COPY} \
        --output ${filtered_repeats_bed}
else
    log_progress "  > Skipping: File already exists: ${filtered_repeats_bed}"
fi

# Step 5: Extract sgRNA candidate loci (23bp, ending with GG)
log_progress "Step 5: Extracting sgRNA candidates..."
sgrna_txt="${OUTPUT_DIR}/sgRNA/sgrna_raw.txt"
sgrna_fa="${sgrna_txt%.txt}.fa"
if [ ! -f ${sgrna_fa} ]; then
    python ${SCRIPTS_DIR}/06.identify_sgrna.py \
        --input ${filtered_repeats_bed} \
        --output ${sgrna_txt} --size 23 --suffix GG
else
    log_progress "  > Skipping: File already exists: ${sgrna_fa}"
fi

# Step 6: Remove off-target sgRNAs
# Args: <output_dir> <sgrna.txt> <bowtie2_idx> [n_cpu] [target_regions_bed]
log_progress "Step 6: remove off-target sgRNAs..."
sgrna_on_target="${OUTPUT_DIR}/sgRNA/sgrna_raw.on_target.txt"
if [ ! -f ${sgrna_on_target} ]; then
    bash ${SCRIPTS_DIR}/07.remove_off_targets.sh \
        ${OUTPUT_DIR}/sgRNA \
        ${sgrna_txt} \
        ${BOWTIE2_IDX} \
        ${N_CPU} \
        ${target_extended_bed}
else
    log_progress "  > Skipping: File already exists: ${sgrna_on_target}"
fi

# Step 7: Filter sgRNAs for better performance
log_progress "Step 7: Filter sgRNAs for better performance..."
# TODO: Implement this step in future versions

# Step 8: Generate summary and report
log_progress "Step 8: Generating report..."
report_dir="${OUTPUT_DIR}/report"
if [ ! -d ${report_dir} ]; then
    mkdir -p ${report_dir}
    # Uncomment when the report script is ready
    # python ${SCRIPTS_DIR}/08.generate_report.py \
    #     --repeats ${filtered_repeats_bed} \
    #     --igh_genes ${TARGET_BED} \
    #     --sgrna ${sgrna_on_target} \
    #     --output ${report_dir}
    echo "Report generation is not implemented yet, skipping..."
else
    log_progress "  > Skipping: Report directory already exists: ${report_dir}"
fi

# Create symlink to the final sgRNA file
sgrna_txt="${OUTPUT_DIR}/sgrna_raw.on_target.txt"
[[ ! -f ${sgrna_txt} ]] && ln -s sgRNA/sgrna_raw.on_target.txt ${OUTPUT_DIR}/

log_progress "Pipeline completed successfully!"
log_progress "Results are saved in: ${OUTPUT_DIR}/sgrna_raw.on_target.txt" 
