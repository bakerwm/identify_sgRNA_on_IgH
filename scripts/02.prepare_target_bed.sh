#!/bin/bash

# Date: 2025-03-25
# Author: Ming Wang
# Email: wangm08@hotmail.com
# Version: 1.0

# Description:
#   Prepare target BED file, slop {left|right} flanking distance
#
# Usage:
#   bash 02.prepare_target_bed.sh <target_bed>

set -e

main () {
    if [ $# -ne 5 ]; then
        echo "Usage: $0 <out_bed> <target_bed> <genome.size> <flank_left> <flank_right>"
        exit 1
    fi
    local out_bed=$1     # output bed file
    local target_bed=$2  # BED3, chr star end
    local genome_fasta=$3   # genome fasta file
    local flank_left=$4  # bp
    local flank_right=$5 # bp

    # remove non-numeric characters
    flank_left=$(echo ${flank_left} | sed 's/[^0-9]//g')
    flank_right=$(echo ${flank_right} | sed 's/[^0-9]//g')

    # prepare output files
    local prefix=$(basename ${out_bed} .bed)
    local output_dir=$(dirname ${out_bed})
    mkdir -p ${output_dir}

    # 1. copy raw_bed file
    local bed_raw=${output_dir}/${prefix}.raw.bed
    # remove ^chr, as Ensembl chromosomes are like: 1, 2, ...
    sed 's/^chr//' ${target_bed} | \
        sed -E 's/ +/\t/g' | \
        sort -k1,1 -k2,2n | \
        uniq > ${bed_raw}

    # 2. slop target bed file
    if [[ ! -f ${genome_fasta}.fai ]] ; then
        samtools faidx ${genome_fasta}
    fi
    bedtools slop -i ${bed_raw} -g ${genome_fasta}.fai -l ${flank_left} -r ${flank_right} > ${out_bed}

    # 3. get fasta file
    local out_fa=${output_dir}/${prefix}.fa
    if [[ ! -f ${out_fa} ]]; then
        bedtools getfasta -fi ${genome_fasta} -bed ${out_bed} -fo ${out_fa}
    fi
}

main $@
