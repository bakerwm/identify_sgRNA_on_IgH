#!/bin/bash

# Date: 2025-03-25
# Author: Ming Wang
# Email: wangm08@hotmail.com
# Version: 1.0

# Description:
#   Download human reference genome (GRCh38.102) and gene annotation (GRCh38.102.gtf)
#
# Usage:
#   bash 01.download_genome.sh <genome_directory> <genome_build>

set -e  # Exit on error

main () {
    if [ $# -lt 2 ]; then
        echo "Usage: $0 <genome_build> <output_dir> [release:102]"
        exit 1
    fi
    local genome_build=$1  # eg: GRCh38, GRCh39
    local genome_dir=$2    # eg: genome/
    local release=${3:-102}      # eg: 102
  
    if [ -z "${release}" ]; then
        case ${genome_build} in
            GRCh38)
                release=102
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
    local url_fa="https://ftp.ensembl.org/pub/release-${release}/fasta/homo_sapiens/dna/Homo_sapiens.${genome_build}.dna.primary_assembly.fa.gz"
    local url_gtf="https://ftp.ensembl.org/pub/release-${release}/gtf/homo_sapiens/Homo_sapiens.${genome_build}.${release}.gtf.gz"
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
    
}

main $@