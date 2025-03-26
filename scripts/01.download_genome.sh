#!/bin/bash

# Date: 2025-03-25
# Author: Ming Wang
# Email: wangm08@hotmail.com
# Version: 1.0

# Description:
#   Download human reference genome (GRCh38.102) and mouse reference genome (GRCm38.102)
#
#  1. download mouse genome, build bowtie2 index
#  2. download human genome, build bowtie2 index for chr21
# Usage:
#   bash 01.download_genome.sh <genome_directory> <genome_build>

set -e  # Exit on error

build_bowtie2_index () {
    [[ $# -ne 2 ]] && { echo "Usage: ${FUNCNAME[0]} <genome_fa> <genome_dir>"; exit 1; }
    local genome_fa=$1
    local bowtie2_idx=$2

    if bowtie2-inspect -n ${bowtie2_idx} > /dev/null 2>&1; then
        echo "Bowtie2 index ${bowtie2_idx} already exists"
    else
        bowtie2-build --threads ${N_CPU} -q ${genome_fa} ${bowtie2_idx} 1>/dev/null
    fi
}

# Download mouse genome, build bowtie2 index
# Download human genome, build bowtie2 index for chr21
# Args: <mouse_build> <human_build> <output_dir> <release>
main () {
    if [[ $# -ne 4 ]]; then
        echo "Usage: $0 <mouse_build> <human_build> <output_dir> <release>"
        exit 1
    fi
    local mouse_build=$1 # eg: GRCm38
    local human_build=$2 # eg: GRCh38
    local output_dir=$3  # eg: genome/
    local release=$4     # eg: 102
    
    # 1. Download mouse genome, build bowtie2 index
    local mouse_fa="${output_dir}/${mouse_build}.fa"
    local mouse_bowtie2_idx="${output_dir}/${mouse_build}"    
    local mouse_fa_url="https://ftp.ensembl.org/pub/release-${release}/fasta/mus_musculus/dna/Mus_musculus.${mouse_build}.dna.primary_assembly.fa.gz"
    local mouse_gtf_url="https://ftp.ensembl.org/pub/release-${release}/gtf/mus_musculus/Mus_musculus.${mouse_build}.${release}.gtf.gz"
    
    # 1.1 Download mouse genome fasta file
    if [ ! -f ${mouse_fa} ]; then
        local fa_file=${output_dir}/$(basename ${mouse_fa_url})
        curl -o ${fa_file} ${mouse_fa_url}
        # unzip fa file
        gunzip -c ${fa_file} > ${mouse_fa}
        # index genome with samtools
        samtools faidx ${mouse_fa}
        # create genome dictionary
        samtools dict ${mouse_fa} > ${mouse_fa}.dict
    fi

    # 1.2 Download mouse genome gtf file
    if [ ! -f ${mouse_gtf} ]; then
        local gtf_file=${output_dir}/$(basename ${mouse_gtf_url})
        curl -o ${gtf_file} ${mouse_gtf_url}
        gunzip -c ${gtf_file} > ${mouse_gtf}
    fi

    # 1.3 Build bowtie2 index for mouse genome
    if bowtie2-inspect -n ${MOUSE_BOWTIE2_IDX} > /dev/null 2>&1; then
        echo "  > Mouse Bowtie2 index: ${MOUSE_BOWTIE2_IDX}"
    else
        echo "  > Building mouse Bowtie2 index..."
        MOUSE_BOWTIE2_IDX=${mouse_bowtie2_idx} # update global variable
        build_bowtie2_index ${mouse_fa} ${MOUSE_BOWTIE2_IDX}
        echo "  > Mouse Bowtie2 index: ${MOUSE_BOWTIE2_IDX}"
    fi

    # 2. Download human genome, build bowtie2 index for chr21
    local human_fa="${output_dir}/${human_build}.fa"
    # local human_chr21_fa="${output_dir}/${human_build}_chr21.fa"
    ## see environment variable: HUMAN_CHR21_FA
    # local human_chr21_bowtie2_idx="${output_dir}/${human_build}_chr21"
    local human_fa_url="https://ftp.ensembl.org/pub/release-${release}/fasta/homo_sapiens/dna/Homo_sapiens.${genome_build}.dna.primary_assembly.fa.gz"
    local human_gtf_url="https://ftp.ensembl.org/pub/release-${release}/gtf/homo_sapiens/Homo_sapiens.${genome_build}.${release}.gtf.gz"
    
    # 2.1 Download human genome fasta file
    if [ ! -f ${human_fa} ]; then
        local fa_file=${output_dir}/$(basename ${human_fa_url})
        curl -o ${fa_file} ${human_fa_url}
        gunzip -c ${fa_file} > ${human_fa}
    fi

    # 2.2 Extract chr21 from human genome fasta file
    echo "  > Extracting chr21 from human genome fasta file..."
    if [ ! -f ${HUMAN_CHR21_FA} ]; then
        samtools faidx ${human_fa} 21 > ${HUMAN_CHR21_FA}
    fi

    # 2.2 Download human genome gtf file
    if [ ! -f ${human_gtf} ]; then
        local gtf_file=${output_dir}/$(basename ${human_gtf_url})
        curl -o ${gtf_file} ${human_gtf_url}
        gunzip -c ${gtf_file} > ${human_gtf}
    fi

    # 2.3 Build bowtie2 index for human genome
    if ! bowtie2-inspect -n ${HUMAN_CHR21_BOWTIE2_IDX} > /dev/null 2>&1; then
        echo "  > Building human chr21 Bowtie2 index..."
        build_bowtie2_index ${HUMAN_CHR21_FA} ${HUMAN_CHR21_BOWTIE2_IDX}
        echo "  > Human chr21 Bowtie2 index: ${HUMAN_CHR21_BOWTIE2_IDX}"
    fi
    echo "  > Human Bowtie2 index: ${HUMAN_CHR21_BOWTIE2_IDX}"
}

main $@