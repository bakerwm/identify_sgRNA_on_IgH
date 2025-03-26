#!/bin/bash

# Date: 2025-03-25
# Author: Ming Wang
# Email: wangm08@hotmail.com
# Version: 1.0

# Description:
#   Configuration file for the sgRNA design pipeline

################################################################################
# BEGIN: DO NOT MODIFY THE FOLLOWING VARIABLES                                  #
# Configuration                                                                #
GENOME_DIR="genome"
OUTPUT_DIR="results"
SCRIPTS_DIR="scripts"
# Create necessary directories
mkdir -p $GENOME_DIR $OUTPUT_DIR $SCRIPTS_DIR
# END: custom variables here (optional)                                        #
################################################################################



################################################################################
# BEGIN: Required VARIABLES                                                    #
# DO NOT REMOVE THE FOLLOWING VARIABLES                                        #
# System configuration
export N_CPU=12

# Genome configuration
export MOUSE_BUILD="GRCm38"
export HUMAN_BUILD="GRCh38"
export RELEASE=102
export MIN_COPY=10      # Minimum copy number of repeat sequences
export MIN_LENGTH=13    # Minimum length of repeat sequences
export HUMAN_CHR21_FA="${GENOME_DIR}/${HUMAN_BUILD}_chr21.fa"
export HUMAN_CHR21_BOWTIE2_IDX="${GENOME_DIR}/${HUMAN_BUILD}_chr21" # HUMAN chr21
export MOUSE_BOWTIE2_IDX="${GENOME_DIR}/${MOUSE_BUILD}" # MOUSE genome
# END: Required VARIABLES                                                      #
################################################################################
