#!/bin/bash

# Date: 2025-03-25
# Author: Ming Wang
# Email: wangm08@hotmail.com
# Version: 1.0

# Description:
# Function to find tandom repeat sequences in the given FASTA file, using TRF
# TRF, Tandom Repeats Finder, is a tool to find tandem repeats in DNA sequences.
#
# TRF command line:
# # TRF (Tandem Repeats Finder) parameters explanation:
# trf <input.fa> 2 7 7 80 10 50 500 -f -d -m -h
# 
# Detailed breakdown of parameters:
# 1. 2   = Matching score (match weight)
# 2. 7   = Mismatch penalty (mismatch weight)
# 3. 7   = Indel penalty (indel weight)
# 4. 80  = Minimum alignment score to report a match
# 5. 10  = Maximum period size to report
# 6. 50  = Minimum percent matches to report
# 7. 500 = Maximum length of a repeat to report
#
# Additional flags:
# -f     = Create additional output files
# -d     = Display alignment details
# -m     = Create HTML alignment summary
# -h     = Display help information
#
# Format of TRF output:
# see documentation: https://tandem.bu.edu/trf/definitions#table

# Example:
# Sequence: 12:111254829-116010093
# 584 644 20 3.3 20 70 25 71 11 31 39 18 1.86 GCCTGGGTCAGCATGGCCCT GCCTGGGTCAGCATGGCCCTGCCTGGGAAGCATGGCCTGGGTCAGCATGGCCCTGCCTGGG

# column-1: start (584)
# column-2: end (644)
# column-3: repeat size (20)
# column-4: copy number (3.3)
# column-5: consensus size (20)
# column-6: percent of matches (70)
# column-7: percent of indels (25)
# column-8: alignment score (71)
# column-9-12: percent composition for each of the four nucleotides, A, C, G, T (11 31 39 18)
# column-13: entropy measure based on percent composition (1.86)
# column-14: consensus sequence (GCCTGGGTCAGCATGGCCCT)
# column-15: repeat sequence (GCCTGGGTCAGCATGGCCCTGCCTGGGAAGCATGGCCTGGGTCAGCATGGCCCTGCCTGGG)

# Source:
# Table Explanation:
# The summary table includes the following information:

# Indices of the repeat relative to the start of the sequence.
# Period size of the repeat.
# Number of copies aligned with the consensus pattern.
# Size of consensus pattern (may differ slightly from the period size).
# Percent of matches between adjacent copies overall.
# Percent of indels between adjacent copies overall.
# Alignment score.
# Percent composition for each of the four nucleotides.
# Entropy measure based on percent composition.

set -e  # Exit on error

# check command
if ! command -v trf &> /dev/null; then
    echo "Error: trf not found, please switch to find_sgrna env"
    exit 1
fi

# Args: <input.fa> <output_dir>
run_trf() {
    [[ $# -ne 2 ]] && {
        echo "Usage: $0 <input_fasta> <output_dir>"
        return 1
    }

    local intput_fa=$1
    local output_dir=$2
    
    echo "  > Running Tandem Repeats Finder (TRF)..."

    trf_out="${output_dir}/trf_output.dat"
    if [ -f ${trf_out} ]; then
        echo "  > TRF Done."
    else 
        mkdir -p ${output_dir}

        # remove existing .dat file
        if ls *.dat > /dev/null 2>&1; then
            rm *.dat
        fi

        # remove existing .mask file
        if ls *.mask > /dev/null 2>&1; then
            rm *.mask
        fi

        # run TRF
        local fa_name=$(basename ${input_fa})
        if trf ${input_fa} 2 7 7 80 10 50 500 -f -d -m -h \
            1> ${output_dir}/trf.stdout \
            2> ${output_dir}/trf.stderr ; then
            echo "  > TRF Done."
        else 
            echo "  > TRF Failed?"
        fi
        # rename the .dat file to trf_output.dat
        mv ${fa_name}.2.7.7.80.10.50.500.dat ${trf_out}
        mv ${fa_name}.2.7.7.80.10.50.500.mask ${output_dir}
    fi
}
export -f run_trf

# Args: <trf_out> <output_dir>
main() {
    [[ $# -lt 2 ]] && {
        echo "Usage: $0 <input.fa> <output.bed> <min_length> <min_copy_number>"
        return 1
    }

    local input_fa=$1
    local output_bed=$2
    local min_length=${3:-13} # default: 13
    local min_copy_number=${4:-10} # default: 10

    local output_dir=$(dirname ${output_bed})
    mkdir -p ${output_dir}

    # 1. Run TRF
    run_trf ${input_fa} ${output_dir}
    local trf_out="${output_dir}/trf_output.dat"

    # Combine results from TRF and RepeatMasker only
    python $(dirname $0)/04.combine_repeats.py \
        --trf ${trf_out} \
        --min_length ${min_length} \
        --min_copy_number ${min_copy_number} \
        --output_bed ${output_bed}

    echo "Repeat finding complete. Results in ${output_dir}" 
}

main $@
