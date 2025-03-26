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
# Args: <off_target_file> <aln.bam> <target_regions_bed>
find_off_targets() {
    [[ $# -ne 3 ]] && {
        echo "Usage: ./map_sgrna_to_genome.sh <off_target_file> <aln.bam> <target_regions_bed>"
        exit 1
    }
    local off_target_file=$1
    local input_bam=$2
    local target_bed=${3:-""}
    
    local prefix="$(basename ${input_bam} .bam)"
    local output_dir="$(dirname ${off_target_file})"
    mkdir -p ${output_dir}

    local aln_bed="${output_dir}/${prefix}.bed"
    if [ ! -f ${aln_bed} ]; then
        bedtools bamtobed -i ${input_bam} > ${aln_bed}
    fi

    # filter alignments mapped only to target regions
    local ext_bed="${output_dir}/${prefix}.off_target.bed"
    if [ -n "${target_bed}" ]; then
        # Extract off-target sgRNAs
        bedtools intersect -v -a ${aln_bed} -b ${target_bed} > ${ext_bed}
        # Extract unique sgRNA ids for off-target sgRNAs
        awk '{print $4}' ${ext_bed} | sort | uniq > ${off_target_file}
    else
        # all alignments are off-target
        awk '{print $4}' ${aln_bed} | sort | uniq > ${off_target_file}
    fi

    # report
    # echo "Off-target sgRNA ids: ${off_target_file}"
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

    # Step1. convert sgrna.txt to fasta format
    # column-9: id, column-10: seq in sgrna_txt file
    local sgrna_fa="${output_dir}/${prefix}.fa"
    awk '{print ">" $9 "\n" $10}' ${sgrna_txt} > ${sgrna_fa}

    # Step2. map sgRNA sequences to reference genome
    local aln_bam="${output_dir}/${prefix}.aln_all.bam"
    map_to_genome ${sgrna_fa} ${aln_bam} ${bowtie2_idx} ${n_cpu}

    # Step3. Extract aligned sgRNAs
    local aligned_id_file="${output_dir}/${prefix}.aligned_ids.txt"
    samtools view ${aln_bam} | awk '{print $1}' | sort | uniq > ${aligned_id_file}

    # Step4. Extract off-target sgRNAs, ids
    local off_target_id_file="${output_dir}/${prefix}.off_target_ids.txt"
    find_off_targets ${off_target_id_file} ${aln_bam} ${target_bed}
    
    # Step5.remove off-target ids
    local on_target_id_file="${output_dir}/${prefix}.on_target_ids.txt"
    local on_target_file="${output_dir}/${prefix}.on_target.txt"
    # cat ${off_target_file} ${aligned_id_file} | sort | uniq -u > ${on_target_ids}
    awk 'NR==FNR {ids[$1]; next} ! ($1 in ids)' ${off_target_id_file} ${aligned_id_file} > ${on_target_id_file}
    # column-9: id, column-10: seq in sgrna_txt file
    awk 'NR==FNR {ids[$1]; next} ($9 in ids)' ${on_target_id_file} ${sgrna_txt} > ${on_target_file}

    # Summary
    local n_input=$(wc -l <${sgrna_txt})
    local n_aln=$(wc -l <${aligned_id_file})
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
