#!/bin/bash

# Date: 2025-03-25
# Author: Ming Wang
# Email: wangm08@hotmail.com
# Version: 1.1

# Description:
#   Configuration file for the chromosome 21 sgRNA design pipeline

################################################################################
# Directory structure                                                           #
################################################################################
GENOME_DIR="genome"
OUTPUT_DIR="results"
SCRIPTS_DIR="scripts"

# Create necessary directories
mkdir -p $GENOME_DIR $OUTPUT_DIR $SCRIPTS_DIR

################################################################################
# System configuration                                                          #
################################################################################
export N_CPU=12      # Number of CPU cores to use for parallel processing

################################################################################
# Genome configuration                                                          #
################################################################################
export MOUSE_BUILD="GRCm38"       # Mouse genome build (used for off-target screening)
export HUMAN_BUILD="GRCh38"       # Human genome build
export RELEASE=102                # Ensembl release number

################################################################################
# Analysis parameters                                                           #
################################################################################
export MIN_COPY=10                # Minimum copy number of repeat sequences
export MIN_LENGTH=13              # Minimum length of repeat sequences

################################################################################
# Genome files and indices                                                      #
################################################################################
# Human chr21 fasta and bowtie2 index
export HUMAN_CHR21_FA="${GENOME_DIR}/${HUMAN_BUILD}_chr21.fa"
export HUMAN_CHR21_BOWTIE2_IDX="${GENOME_DIR}/${HUMAN_BUILD}_chr21" 

# Mouse genome bowtie2 index (for off-target screening)
export MOUSE_BOWTIE2_IDX="${GENOME_DIR}/${MOUSE_BUILD}" 
