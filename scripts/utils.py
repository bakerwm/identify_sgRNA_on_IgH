#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Utility functions for the chr21 sgRNA design pipeline.

This module contains common utility functions used across the pipeline scripts.

Author: Ming Wang
Version: 1.0
Date: 2025-03-26
"""

import os
import sys
import logging
import subprocess
from typing import List, Dict, Tuple, Union, Optional, Any
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq


def setup_logger(log_file: str = None, log_level: int = logging.INFO) -> logging.Logger:
    """Set up and return a logger with file and console handlers.
    
    Args:
        log_file: Path to log file, if None, only console output is used
        log_level: Logging level (default: INFO)
        
    Returns:
        Logger object
    """
    # Create logger
    logger = logging.getLogger('chr21_sgRNA_pipeline')
    logger.setLevel(log_level)
    
    # Create formatters
    file_formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    console_formatter = logging.Formatter('%(levelname)s: %(message)s')
    
    # Create console handler
    console_handler = logging.StreamHandler()
    console_handler.setLevel(log_level)
    console_handler.setFormatter(console_formatter)
    logger.addHandler(console_handler)
    
    # Create file handler if log_file is provided
    if log_file:
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(log_level)
        file_handler.setFormatter(file_formatter)
        logger.addHandler(file_handler)
    
    return logger


def run_command(cmd: List[str], logger: logging.Logger = None) -> Tuple[int, str, str]:
    """Run a shell command and handle errors.
    
    Args:
        cmd: Command to run as a list of strings
        logger: Logger object for messages (optional)
        
    Returns:
        Tuple of (return_code, stdout, stderr)
    """
    if logger:
        logger.info(f"Running command: {' '.join(cmd)}")
    
    process = subprocess.Popen(
        cmd, 
        stdout=subprocess.PIPE, 
        stderr=subprocess.PIPE, 
        text=True
    )
    
    stdout, stderr = process.communicate()
    
    if process.returncode != 0 and logger:
        logger.error(f"Command failed with return code {process.returncode}")
        logger.error(f"Error message: {stderr}")
    
    return process.returncode, stdout, stderr


def load_fasta(fasta_file: str) -> Dict[str, str]:
    """Load a FASTA file and return a dictionary of sequences.
    
    Args:
        fasta_file: Path to FASTA file
        
    Returns:
        Dictionary mapping sequence IDs to sequences
    """
    if not os.path.exists(fasta_file):
        raise FileNotFoundError(f"FASTA file not found: {fasta_file}")
    
    sequences = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences[record.id] = str(record.seq)
    
    return sequences


def load_bed(bed_file: str, header: List[str] = None) -> pd.DataFrame:
    """Load a BED file into a pandas DataFrame.
    
    Args:
        bed_file: Path to BED file
        header: Column names (optional)
        
    Returns:
        DataFrame with BED entries
    """
    if not os.path.exists(bed_file):
        raise FileNotFoundError(f"BED file not found: {bed_file}")
    
    try:
        if header is None:
            df = pd.read_csv(bed_file, sep='\t', header=None)
        else:
            df = pd.read_csv(bed_file, sep='\t', names=header)
        return df
    except Exception as e:
        print(f"Error parsing BED file {bed_file}: {e}")
        return None


def write_bed(df: pd.DataFrame, output_file: str, header: bool = False) -> None:
    """Write a DataFrame to a BED file.
    
    Args:
        df: DataFrame to write
        output_file: Path to output file
        header: Whether to include header
    """
    df.to_csv(output_file, sep='\t', header=header, index=False)


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence.
    
    Args:
        seq: DNA sequence
        
    Returns:
        Reverse complement of seq
    """
    return str(Seq(seq).reverse_complement())


def find_sgRNAs(sequence: str, pam_seq: str = "GG", sgrna_length: int = 20) -> List[Dict]:
    """Find all possible sgRNAs in a sequence.
    
    Args:
        sequence: DNA sequence to search
        pam_seq: PAM sequence (default: GG for SpCas9)
        sgrna_length: Length of sgRNA without PAM (default: 20)
        
    Returns:
        List of dictionaries with sgRNA information
    """
    sgrnas = []
    
    # Search forward strand
    for i in range(len(sequence) - sgrna_length - len(pam_seq) + 1):
        window = sequence[i:i+sgrna_length+len(pam_seq)]
        if window[-len(pam_seq):] == pam_seq:
            sgrnas.append({
                'start': i,
                'end': i + sgrna_length + len(pam_seq),
                'seq': window,
                'strand': '+'
            })
    
    # Search reverse strand
    rev_comp = reverse_complement(sequence)
    for i in range(len(rev_comp) - sgrna_length - len(pam_seq) + 1):
        window = rev_comp[i:i+sgrna_length+len(pam_seq)]
        if window[-len(pam_seq):] == pam_seq:
            original_pos = len(sequence) - (i + sgrna_length + len(pam_seq))
            sgrnas.append({
                'start': original_pos,
                'end': original_pos + sgrna_length + len(pam_seq),
                'seq': window,
                'strand': '-'
            })
    
    return sgrnas


def calculate_gc_content(seq: str) -> float:
    """Calculate the GC content of a DNA sequence.
    
    Args:
        seq: DNA sequence
        
    Returns:
        GC content as a percentage
    """
    g_count = seq.upper().count('G')
    c_count = seq.upper().count('C')
    return (g_count + c_count) / len(seq) * 100


def ensure_directory(directory: str) -> None:
    """Ensure a directory exists, create if it doesn't.
    
    Args:
        directory: Directory path
    """
    if not os.path.exists(directory):
        os.makedirs(directory)


def is_valid_sgRNA(seq: str) -> bool:
    """Check if a sequence is a valid sgRNA (no homopolymers, balanced GC).
    
    Args:
        seq: DNA sequence of sgRNA
        
    Returns:
        True if valid, False otherwise
    """
    # Check for homopolymers (4+ of same base)
    for base in 'ATGC':
        if base * 4 in seq:
            return False
    
    # Check GC content is in reasonable range (30-70%)
    gc_content = calculate_gc_content(seq)
    if gc_content < 30 or gc_content > 70:
        return False
    
    return True


def main():
    """Main function for testing the module."""
    print("This is a utility module and not meant to be run directly.")
    print("Import the functions from this module in your scripts.")


if __name__ == "__main__":
    main() 