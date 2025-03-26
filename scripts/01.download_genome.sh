#!/bin/bash

# Date: 2025-03-25
# Author: Ming Wang
# Email: wangm08@hotmail.com
# Version: 1.0

# Description:
#   Download genome (fasta)
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
        bowtie2-build --threads ${N_CPU} -q ${genome_fa} ${bowtie2_idx} 1> /dev/null
    fi
}

# Get Ensembl URL to download fasta or gtf file
# supported:
# human (GRCh39)
# mouse (GRCm38, GRCm39)
# fruitfly (BDGP6)
# Args: genome_build, release, [type: fasta|gtf]
get_ensembl_url () {
    [[ $# -ne 3 ]] && { 
        echo "Usage: $0 <genome_build> <release> [type: fasta|gtf]"
        exit 1
    }
    local genome_build=$1
    local release=$2
    local type=$3     # fasta or gtf

    # urls
    local url_root="https://ftp.ensembl.org/pub/release-${release}"

    case ${genome_build} in
        GRCh39)
            local species="homo_sapiens"
            local dna_level="primary_assembly"
            # local url_fa="${url_root}/fasta/homo_sapiens/dna/Homo_sapiens.${genome_build}.dna.primary_assembly.fa.gz"
            # local url_gtf="${url_root}/gtf/homo_sapiens/Homo_sapiens.${genome_build}.${release}.gtf.gz"
            ;;
        GRCm38|GRCm39)
            local species="mus_musculus"
            local dna_level="primary_assembly"
            # local url_fa="${url_root}/fasta/mus_musculus/dna/Mus_musculus.${genome_build}.dna.primary_assembly.fa.gz"
            # local url_gtf="${url_root}/gtf/mus_musculus/Mus_musculus.${genome_build}.${release}.gtf.gz"
            ;;
        BDGP6)
            genome_build="${genome_build}.28" # update to latest release
            local species="drosophila_melanogaster"
            local dna_level="toplevel"
            # local url_fa="${url_root}/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.${genome_build}.28.dna.toplevel.fa.gz"
            # local url_gtf="${url_root}/gtf/drosophila_melanogaster/Drosophila_melanogaster.${genome_build}.28.${release}.gtf.gz"
            ;;
        *)
            echo "Error: Unsupported genome build: ${genome_build}"
            exit 1
    esac

    # ${var^}:  capitalize the first letter
    # ${var^^}: capitalize all the letters
    # ${var,}:  lowercase the first letter
    # ${var,,}: lowercase all the letters
    local url_fa="${url_root}/fasta/${species}/dna/${species^}.${genome_build}.dna.${dna_level}.fa.gz"
    local url_gtf="${url_root}/gtf/${species}/${species^}.${genome_build}.${release}.gtf.gz"
    
    # echo "${url_fa} ${url_gtf}"
    case ${type,,} in
        fasta|fa)
            echo "${url_fa}"
            ;;
        gtf)
            echo "${url_gtf}"
            ;;
        *)
            echo "Error: Unsupported type: ${type}, expected: fasta|gtf"
            exit 1
    esac
}

main () {
    if [ $# -lt 2 ]; then
        echo "Usage: $0 <genome_build> <output_dir> [release:102]"
        exit 1
    fi
    local genome_build=$1   # eg: GRCh38
    local genome_dir=$2     # eg: genome/
    local release=${3:-102} # eg: 102
      
    # download genome fasta file
    mkdir -p ${genome_dir}
    local genome_fa="${genome_dir}/${genome_build}.fa"
    local url_fa=$(get_ensembl_url ${genome_build} ${release} fasta)

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

    # # download gene annotation file
    # local genome_gtf="${genome_dir}/${genome_build}.gtf"
    # local url_gtf=$(get_ensembl_url ${genome_build} ${release} gtf)
    # if [ ! -f ${genome_gtf} ]; then
    #     gtf_file=${genome_dir}/$(basename ${url_gtf})
    #     curl -o ${gtf_file} ${url_gtf}
    #     gunzip -c ${gtf_file} > ${genome_gtf}
    # fi

    # build bowtie2 index
    local bowtie2_idx=${genome_dir}/${genome_build}
    if ! bowtie2-inspect -n ${bowtie2_idx} > /dev/null 2>&1; then
        echo "  > Building bowtie2 index..."
        build_bowtie2_index ${genome_fa} ${bowtie2_idx}
        export BOWTIE2_IDX=${bowtie2_idx} # update global variable
        echo "  > Bowtie2 index: ${BOWTIE2_IDX}"
    fi
}

main $@
