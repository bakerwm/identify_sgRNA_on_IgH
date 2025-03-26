#!/bin/bash

# Date: 2025-03-25
# Author: Ming Wang
# Email: wangm08@hotmail.com
# Version: 1.1

# Description:
# on-target sgRNAs are mapped to human chr21
# off-target sgRNAs, mapped to mouse genome
#
# usage: ./map_sgrna_to_genome.sh <sgrna_fasta> <target_regions_bed> <output_prefix>

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

# Function to log messages with timestamps
log_progress() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

# Map sgRNA sequences to reference genome
# Input: sgRNA sequences in fasta format
# Output: alignment in bed format
# Args: <sgrna.fa> <aln.bam> <bowtie2_idx> [n_cpu]
map_to_genome() {
    [[ $# -lt 3 ]] && {
        log_progress "ERROR: Insufficient arguments for map_to_genome function"
        echo "Usage: ./map_sgrna_to_genome.sh <sgrna.fa> <aln.bam> <bowtie2_idx> [n_cpu]"
        exit 1
    }
    local input_fa=$1
    local output_bam=$2
    local bowtie2_idx=$3
    local n_cpu=${4:-8} # default 8 threads
    # see environment variable N_CPU and BOWTIE2_IDX in workflow.sh

    # Check if input file exists
    if [[ ! -f ${input_fa} ]]; then
        log_progress "ERROR: Input FASTA file not found: ${input_fa}"
        exit 1
    fi

    # Check if bowtie2 index exists
    if [[ ! -f ${bowtie2_idx}.1.bt2 && ! -f ${bowtie2_idx}.1.bt2l ]]; then
        log_progress "ERROR: Bowtie2 index not found: ${bowtie2_idx}"
        exit 1
    fi

    # report all alignments
    local log_file="${output_bam%.bam}.bowtie2.log"
    if [ ! -f ${output_bam} ]; then
        log_progress "Mapping sgRNAs to reference using bowtie2: ${bowtie2_idx}"
        # -a: output all alignments
        bowtie2 -a -f -p ${n_cpu} --end-to-end --no-unal \
            -x ${bowtie2_idx} -U ${input_fa} 2> ${log_file} | \
            samtools view -bhS - | \
            samtools sort -o ${output_bam} -
        samtools index ${output_bam}
        log_progress "Alignment complete: ${output_bam}"
    else
        log_progress "Using existing alignment file: ${output_bam}"
    fi
}

# Find off-target sgRNAs
# if target_regions_bed is provided, only off-target sgRNAs mapped to target regions are reported
# otherwise, all off-target sgRNAs are reported
# Args: <off_target_file> <aln.bam> [target_regions_bed]
find_off_targets() {
    [[ $# -lt 2 ]] && {
        log_progress "ERROR: Insufficient arguments for find_off_targets function"
        echo "Usage: ./find_off_targets <off_target_file> <aln.bam> [target_regions_bed]"
        exit 1
    }
    local off_target_file=$1
    local input_bam=$2
    local target_bed=${3:-""}
    
    # Check if input file exists
    if [[ ! -f ${input_bam} ]]; then
        log_progress "ERROR: Input BAM file not found: ${input_bam}"
        exit 1
    fi
    
    local prefix="$(basename ${input_bam} .bam)"
    local output_dir="$(dirname ${off_target_file})"
    mkdir -p ${output_dir}

    if [ -n "${target_bed}" ]; then
        # Check if target bed exists
        if [[ ! -f ${target_bed} ]]; then
            log_progress "ERROR: Target BED file not found: ${target_bed}"
            exit 1
        fi
        
        log_progress "Extracting sgRNAs that map outside target regions"
        # Extract off-target sgRNAs, mapped to outside target regions
        local ext_bed="${output_dir}/${prefix}.off_target.bed"
        bedtools intersect -v -bed -a ${input_bam} -b ${target_bed} | \
            awk '{print $4}' | sort | uniq > ${off_target_file}
    else
        log_progress "Extracting all mapped sgRNAs"
        # Extract all alignments as off-target sgRNAs
        samtools view ${input_bam} | awk '{print $1}' | sort | uniq > ${off_target_file}
    fi
    log_progress "Found $(wc -l < ${off_target_file}) sgRNAs in ${off_target_file}"
}

# Args: <output_dir> <sgrna.txt> <mouse_bowtie2_idx> <human_chr21_bowtie2_idx> [n_cpu]
main() {
    [[ $# -lt 4 ]] && {
        log_progress "ERROR: Insufficient arguments for main function"
        echo "Usage: ./remove_off_targets.sh <output_dir> <sgrna.txt> <mouse_bowtie2_idx> <human_chr21_bowtie2_idx> [n_cpu]"
        exit 1
    }
    local output_dir=$1
    local sgrna_txt=$2  # including id and seq in column-9 and column-10
    local mouse_bowtie2_idx=$3
    local human_chr21_bowtie2_idx=$4
    local n_cpu=${5:-12} # default 12 threads

    # Validate input
    if [[ ! -f ${sgrna_txt} ]]; then
        log_progress "ERROR: sgRNA file not found: ${sgrna_txt}"
        exit 1
    fi

    local prefix="$(basename ${sgrna_txt} .txt)"
    mkdir -p ${output_dir}

    log_progress "Starting off-target removal for sgRNAs"
    log_progress "Input: ${sgrna_txt}"
    
    # Step1. convert sgrna.txt to fasta format
    # column-9: id, column-10: seq in sgrna_txt file
    local sgrna_fa="${output_dir}/${prefix}.fa"
    log_progress "Converting sgRNA file to FASTA format: ${sgrna_fa}"
    awk '{print ">" $9 "\n" $10}' ${sgrna_txt} > ${sgrna_fa}

    ############################################################################
    # on-target alignment, map to human chr21                                  #
    ############################################################################
    log_progress "Performing on-target alignment to human chromosome 21"
    # Step2a. map to target genome (human chr21)
    local human_chr21_bam="${output_dir}/${prefix}.human_chr21.bam"
    local human_chr21_id_file="${output_dir}/${prefix}.human_chr21_ids.txt"
    map_to_genome ${sgrna_fa} ${human_chr21_bam} ${human_chr21_bowtie2_idx} ${n_cpu}
    log_progress "Identifying sgRNAs that map to human chromosome 21"
    find_off_targets ${human_chr21_id_file} ${human_chr21_bam} # sgRNAs mapped to human chr21

    ############################################################################
    # off-target alignment, map to mouse genome                                #
    ############################################################################
    log_progress "Performing off-target alignment to mouse genome"
    # Step3a. map to off-target genome (mouse genome)
    local mouse_bam="${output_dir}/${prefix}.mouse.bam"
    local mouse_id_file="${output_dir}/${prefix}.mouse_ids.txt"
    map_to_genome ${sgrna_fa} ${mouse_bam} ${mouse_bowtie2_idx} ${n_cpu}
    log_progress "Identifying sgRNAs that map to mouse genome"
    find_off_targets ${mouse_id_file} ${mouse_bam} # sgRNAs mapped to mouse genome

    ############################################################################
    # remove off-target ids                                                    #
    ############################################################################
    log_progress "Identifying on-target sgRNAs (map to human chr21 but not mouse genome)"
    # Step4. remove sgRNAs mapped to mouse genome, and keep sgRNAs mapped to human chr21
    local on_target_id_file="${output_dir}/${prefix}.on_target_ids.txt"
    awk 'NR==FNR {ids[$1]; next} ! ($1 in ids)' ${mouse_id_file} ${human_chr21_id_file} > ${on_target_id_file}

    # Step5. Extract on-target sgRNAs
    local on_target_file="${output_dir}/${prefix}.on_target.txt"
    log_progress "Extracting on-target sgRNAs from original file"
    # column-9: id, column-10: seq in sgrna_txt file
    awk 'NR==FNR {ids[$1]; next} ($9 in ids)' ${on_target_id_file} ${sgrna_txt} > ${on_target_file}

    # Summary
    local n_input=$(wc -l <${sgrna_txt})
    local n_aln=$(wc -l <${human_chr21_id_file})
    local n_off_target=$(wc -l <${mouse_id_file})
    local n_on_target=$(wc -l <${on_target_id_file})
    
    # report
    log_progress "Off-target removal complete"
    log_progress "  > sgRNA_stat"
    log_progress "  > input:      ${n_input}"
    log_progress "  > aligned:    ${n_aln}"
    log_progress "  > off-target: ${n_off_target}"
    log_progress "  > on-target:  ${n_on_target}"
    log_progress "On-target sgRNAs saved to: ${on_target_file}"
}

# Execute main function with all arguments
main "$@"
