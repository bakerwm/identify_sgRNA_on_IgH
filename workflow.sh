#!/bin/bash

# Date: 2025-03-25
# Author: Ming Wang
# Email: wangm08@hotmail.com
# Version: 1.1

# Description:
#
# Design sgRNAs (SpCas9, NGG) targeting chromosome 21 of human genome
#
# How-To:
# 1. Identify chromosome 21 of human genome (GRCh38, Ensembl release 102)
#    - chr21 (GRCh38/hg38)
# 2. Identify tandom repeat sequences
# 3. Design sgRNAs targeting the random repeat sequences, on consensus sequence
# 4. Remove off-target sgRNAs (bowtie2 alignment to MOUSE genome)
# 5. Generate report
#
# Input:
#    - Human reference genome (GRCh38/hg38)
#    - MOUSE genome (GRCm38/mm10)
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

# Setup logging
LOG_FILE="${OUTPUT_DIR}/workflow_$(date '+%Y%m%d_%H%M%S').log"
mkdir -p ${OUTPUT_DIR}

# Redirect all output to both terminal and log file
exec > >(tee -a "${LOG_FILE}") 2>&1

log_progress "Starting sgRNA design pipeline for human chromosome 21"
log_progress "Configuration: ${config_file}"

# Step 0: Validate required tools
log_progress "Step 0: Checking required tools and dependencies"
for cmd in bowtie2 samtools bedtools trf python; do
    if ! command -v $cmd &> /dev/null; then
        echo "Error: Required tool '$cmd' not found in PATH"
        echo "Please install required dependencies or activate conda environment"
        exit 1
    fi
done

# Step 1: Download human reference genome if not already present
log_progress "Step 1: Preparing human reference genome, annotation-GTF..."
bash ${SCRIPTS_DIR}/01.download_genome.sh ${MOUSE_BUILD} ${HUMAN_BUILD} ${GENOME_DIR} ${RELEASE}

# Step 2: Find repeat sequences in the extracted regions
# Using only TRF and RepeatMasker (skipping k-mer analysis and MISA)
log_progress "Step 2: Identifying repeat sequences..."
raw_repeats_bed="${OUTPUT_DIR}/repeats/all_repeats.raw.bed"
if [ ! -f ${raw_repeats_bed} ]; then
    bash ${SCRIPTS_DIR}/03.find_repeats.sh \
        ${HUMAN_CHR21_FA} \
        ${raw_repeats_bed} \
        ${MIN_LENGTH} \
        ${MIN_COPY}
else
    log_progress "  > Skipping: File already exists: ${raw_repeats_bed}"
fi

# Step 3: Filter repeats based on criteria (length > 13bp, copy number > 10)
log_progress "Step 3: Filtering repeats based on criteria..."
filtered_repeats_bed="${OUTPUT_DIR}/repeats/all_repeats.filtered.bed"
if [ ! -f ${filtered_repeats_bed} ]; then
    python ${SCRIPTS_DIR}/05.filter_repeats.py \
        --input ${raw_repeats_bed} \
        --min_length ${MIN_LENGTH} \
        --min_copy_number ${MIN_COPY} \
        --output ${filtered_repeats_bed}
else
    log_progress "  > Skipping: File already exists: ${filtered_repeats_bed}"
fi

# Step 4: Extract sgRNA candidate loci (23bp, ending with GG)
log_progress "Step 4: Extracting sgRNA candidates..."
sgrna_txt="${OUTPUT_DIR}/sgRNA/sgrna_raw.txt"
sgrna_fa="${sgrna_txt%.txt}.fa"
if [ ! -f ${sgrna_fa} ]; then
    python ${SCRIPTS_DIR}/06.identify_sgrna.py \
        --input ${filtered_repeats_bed} \
        --output ${sgrna_txt} --size 23 --suffix GG
else
    log_progress "  > Skipping: File already exists: ${sgrna_fa}"
fi

# Step 5: Map sgRNA sequences to reference genome
# Args: <output_dir> <sgrna.txt> <bowtie2_idx> [n_cpu] [target_regions_bed]
log_progress "Step 5: Removing off-target sgRNAs..."
sgrna_on_target="${OUTPUT_DIR}/sgRNA/sgrna_raw.on_target.txt"
if [ ! -f ${sgrna_on_target} ]; then
    bash ${SCRIPTS_DIR}/07.remove_off_targets.sh \
        ${OUTPUT_DIR}/sgRNA \
        ${sgrna_txt} \
        ${MOUSE_BOWTIE2_IDX} \
        ${HUMAN_CHR21_BOWTIE2_IDX} \
        ${N_CPU}
else
    log_progress "  > Skipping: File already exists: ${sgrna_on_target}"
fi

# Step 6: Generate summary and report
log_progress "Step 6: Generating report..."
report_dir="${OUTPUT_DIR}/report"
if [ ! -d ${report_dir} ]; then
    mkdir -p ${report_dir}
    # Uncomment when the report script is ready
    # python ${SCRIPTS_DIR}/08.generate_report.py \
    #     --repeats ${filtered_repeats_bed} \
    #     --sgrna ${sgrna_on_target} \
    #     --output ${report_dir}
    log_progress "  > Report generation is not implemented yet, skipping..."
else
    log_progress "  > Skipping: Report directory already exists: ${report_dir}"
fi

# Create symlink to the final sgRNA file
sgrna_txt="${OUTPUT_DIR}/sgrna_raw.on_target.txt"
[[ ! -f ${sgrna_txt} ]] && ln -s sgRNA/sgrna_raw.on_target.txt ${OUTPUT_DIR}/

log_progress "Pipeline completed successfully!"
log_progress "Results are available in ${OUTPUT_DIR}/sgrna_raw.on_target.txt" 
