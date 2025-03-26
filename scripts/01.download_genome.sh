#!/bin/bash

# Date: 2025-03-25
# Author: Ming Wang
# Email: wangm08@hotmail.com
# Version: 1.0

# Description:
#   Download mouse reference genome (GRCm38.102) and gene annotation (GRCm38.102.gtf)
#   Build bowtie2 index for the genome if required
#
# Usage:
#   bash 01.download_genome.sh <genome_directory> <genome_build>

set -e  # Exit on error

build_bowtie2_index () {
    [[ $# -ne 2 ]] && { echo "Usage: $0 <genome_fa> <genome_dir>"; exit 1; }
    local genome_fa=$1
    local bowtie2_idx=$2

    if bowtie2-inspect -n ${bowtie2_idx} > /dev/null 2>&1; then
        echo "Bowtie2 index ${bowtie2_idx} already exists"
    else
        bowtie2-build --threads ${N_CPU} -q ${genome_fa} ${bowtie2_idx}
    fi
}

main () {
    if [ $# -lt 2 ]; then
        echo "Usage: $0 <genome_build> <output_dir> [release:102]"
        exit 1
    fi
    local genome_build=$1  # eg: GRCm38, GRCm39
    local genome_dir=$2    # eg: genome/
    local release=${3:-102}      # eg: 102
  
    if [ -z "${release}" ]; then
        case ${genome_build} in
            GRCm38)
                release=102
                ;;
            GRCm39)
                release=111
                ;;
            *)
                echo "Error: Unsupported genome build: ${genome_build}"
                exit 1
                ;;
        esac
    fi
    
    local genome_fa="${genome_dir}/${genome_build}.fa"
    local genome_gtf="${genome_dir}/${genome_build}.gtf"

    # urls of Ensembl release
    local release=102
    local url_fa="https://ftp.ensembl.org/pub/release-${release}/fasta/mus_musculus/dna/Mus_musculus.${genome_build}.dna.primary_assembly.fa.gz"
    local url_gtf="https://ftp.ensembl.org/pub/release-${release}/gtf/mus_musculus/Mus_musculus.${genome_build}.${release}.gtf.gz"

    # download genome fasta file
    mkdir -p ${genome_dir}

    if [ ! -f ${genome_fa} ]; then
        fa_file=${genome_dir}/$(basename ${url_fa})
        curl -o ${fa_file} ${url_fa}
        # unzip fa file
        gunzip -c ${fa_file} > ${genome_fa}
        # index genome with samtools
        samtools faidx ${genome_fa}
        # create genome dictionary
        samtools dict ${genome_fa} > ${genome_fa}.dict
    fi

    # download gene annotation file
    if [ ! -f ${genome_gtf} ]; then
        gtf_file=${genome_dir}/$(basename ${url_gtf})
        curl -o ${gtf_file} ${url_gtf}
        gunzip -c ${gtf_file} > ${genome_gtf}
    fi

    # build bowtie2 index
    local bowtie2_idx=${genome_dir}/${genome_build}
    if bowtie2-inspect -n ${BOWTIE2_IDX} > /dev/null 2>&1; then
        echo "  > Bowtie2 index: ${BOWTIE2_IDX}"
    else
        echo "  > Building bowtie2 index..."
        build_bowtie2_index ${genome_fa} ${bowtie2_idx}
        BOWTIE2_IDX=${bowtie2_idx} # update global variable
        echo "  > Bowtie2 index: ${BOWTIE2_IDX}"
    fi
}

main $@