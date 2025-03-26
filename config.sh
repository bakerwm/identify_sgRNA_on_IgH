#!/bin/bash

# Date: 2025-03-25
# Author: Ming Wang
# Email: wangm08@hotmail.com
# Version: 1.0

# Description:
#   Configuration file for the sgRNA design pipeline

################################################################################
# BEGIN: custom variables here (optional)                                      #

# END: custom variables here (optional)                                        #
################################################################################



################################################################################
# BEGIN: Required VARIABLES                                                    #
# DO NOT REMOVE THE FOLLOWING VARIABLES                                        #
# System configuration
export N_CPU=12

# Genome configuration
export GENOME_BUILD="GRCh38"
export RELEASE=102
export BOWTIE2_IDX="genome/${GENOME_BUILD}" 

# Analysis parameters
export TARGET_BED="data/igh_regions.human.bed"  # Required, BED3 format
export FLANKING_LEFT=0         # 2 Mbp around IgH genes
export FLANKING_RIGHT=2000000              # 2 Mbp around IgH genes
export MIN_COPY=10                   # Minimum copy number of repeat sequences
export MIN_LENGTH=13                 # Minimum length of repeat sequences
# END: Required VARIABLES                                                      #
################################################################################
