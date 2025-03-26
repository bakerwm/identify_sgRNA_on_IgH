#!/bin/bash

# Date: 2025-03-25
# Author: Ming Wang
# Email: wangm08@hotmail.com
# Version: 1.0

# Description:
# Remove off-target sgRNAs
# by mapping sgRNA sequences to reference genome
#
# usage: ./map_sgrna_to_genome.sh <sgrna_fasta> <target_regions_bed> <output_prefix>

# see 

# Map sgRNA sequences to reference genome
# Input: sgRNA sequences in fasta format
# Output: alignment in bed format
# Args: <sgrna_fasta> <aln.bed>
map_to_genome() {
    [[ $# -lt 3 ]] && {
        echo "Usage: ./map_sgrna_to_genome.sh <sgrna.fa> <aln.bam> <bowtie2_idx> [n_cpu]"
        exit 1
    }
    local input_fa=$1
    local output_bam=$2
    local bowtie2_idx=$3
    local n_cpu=${4:-8} # default 8 threads
    # see environment variable N_CPU and BOWTIE2_IDX in workflow.sh

    # report all alignments
    local log_file="${output_bam%.bam}.bowtie2.log"
    if [ ! -f ${output_bam} ]; then
        # -a: output all alignments
        bowtie2 -a -f -p ${n_cpu} --end-to-end --no-unal \
            -x ${bowtie2_idx} -U ${input_fa} 2> ${log_file} | \
            samtools view -bhS - | \
            samtools sort -o ${output_bam} -
        samtools index ${output_bam}
    fi
}

# Find off-target sgRNAs
# if target_regions_bed is provided, only off-target sgRNAs mapped to target regions are reported
# otherwise, all off-target sgRNAs are reported
# Args: <off_target_file> <aln.bam> [target_regions_bed]
find_off_targets() {
    [[ $# -lt 2 ]] && {
        echo "Usage: ./find_off_targets <off_target_file> <aln.bam> [target_regions_bed]"
        exit 1
    }
    local off_target_file=$1
    local input_bam=$2
    local target_bed=${3:-""}
    
    local output_dir="$(dirname ${off_target_file})"
    mkdir -p ${output_dir}

    if [[ -f "${target_bed}" ]]; then
        # Extract off-target sgRNAs, mapped to outside target regions
        bedtools intersect -v -bed -a ${input_bam} -b ${target_bed} | \
            awk '{print $4}' | sort | uniq > ${off_target_file}
    else
        # Extract all alignments as off-target sgRNAs
        samtools view ${input_bam} | awk '{print $1}' | sort | uniq > ${off_target_file}
    fi
}

# Args: <output_dir> <sgrna.txt> <bowtie2_idx> [n_cpu] [target_regions_bed]
main() {
    [[ $# -lt 3 ]] && {
        echo "Usage: ./map_sgrna_to_genome.sh <output_dir> <sgrna.txt> <bowtie2_idx> [n_cpu] [target_regions_bed]"
        exit 1
    }
    local output_dir=$1
    local sgrna_txt=$2  # including id and seq in column-9 and column-10
    local bowtie2_idx=$3
    local n_cpu=${4:-8} # default 8 threads
    local target_bed=${5:-""} # optional

    local prefix="$(basename ${sgrna_txt} .txt)"
    mkdir -p ${output_dir}

    # Step 1: convert sgrna.txt to fasta format
    # column-9: id, column-10: seq in sgrna_txt file
    local sgrna_fa="${output_dir}/${prefix}.fa"
    awk '{print ">" $9 "\n" $10}' ${sgrna_txt} > ${sgrna_fa}

    # Step 2: map all sgRNAs to whole genome
    local aln_bam="${output_dir}/${prefix}.aln_all.bam"
    map_to_genome ${sgrna_fa} ${aln_bam} ${bowtie2_idx} ${n_cpu}
    local aln_id_file="${output_dir}/${prefix}.aln_all_ids.txt"
    samtools view ${aln_bam} | awk '{print $1}' | sort | uniq > ${aln_id_file}

    # Step 3: Extract aligned sgRNAs
    local off_target_id_file="${output_dir}/${prefix}.off_target_ids.txt"
    find_off_targets ${off_target_id_file} ${aln_bam} ${target_bed}

    # Step 4: Extract on-target sgRNAs ids
    local on_target_id_file="${output_dir}/${prefix}.on_target_ids.txt"
    awk 'NR==FNR {ids[$1]; next} ! ($1 in ids)' ${off_target_id_file} ${aln_id_file} > ${on_target_id_file}

    # Step 5: Extract on-target sgRNAs
    local on_target_file="${output_dir}/${prefix}.on_target.txt"
    awk 'NR==FNR {ids[$1]; next} ($9 in ids)' ${on_target_id_file} ${sgrna_txt} > ${on_target_file}
    
    # Step 6: Summary
    local n_input=$(wc -l <${sgrna_txt})
    local n_aln=$(wc -l <${aln_id_file})
    local n_off_target=$(wc -l <${off_target_id_file})
    local n_on_target=$(wc -l <${on_target_id_file})
    # report
    echo "  > sgRNA_stat"
    echo "  > input:      ${n_input}"
    echo "  > aligned:    ${n_aln}"
    echo "  > off-target: ${n_off_target}"
    echo "  > on-target:  ${n_on_target}"
}

main $@
