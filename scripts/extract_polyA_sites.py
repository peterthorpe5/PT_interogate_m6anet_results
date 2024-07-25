import pysam
import re
import pandas as pd
import argparse
import os
from scipy.stats import wasserstein_distance, mannwhitneyu
import numpy as np
import logging
from statsmodels.stats.multitest import multipletests

def get_args():
    parser = argparse.ArgumentParser(description="Extract poly(A) sites from nanopore direct RNAseq data",
                                     add_help=False)
    file_directory = os.path.realpath(__file__).split("extract_polyA_sites.py")[0]
    optional = parser.add_argument_group('optional arguments')

    optional.add_argument("--bam", dest='bam',
                          action="store",
                          nargs='+',
                          required=True,
                          type=str,
                          help="List of BAM files to be parsed e.g. --bam file1.bam file2.bam file3.bam")

    optional.add_argument("--output", dest='output',
                          action="store",
                          type=str,
                          required=True,
                          help="Output TSV file to store the poly(A) sites")

    optional.add_argument("--fasta", dest='fasta',
                          action="store",
                          type=str,
                          required=True,
                          help="FASTA file of the reference genome")

    optional.add_argument("--groups", dest='groups',
                          action="store",
                          nargs='+',
                          required=True,
                          type=str,
                          help="List of group names corresponding to BAM files e.g. --groups WT WT WT MUT MUT MUT")

    optional.add_argument("--fdr", dest='fdr',
                          action="store",
                          type=float,
                          default=0.05,
                          help="False discovery rate threshold for multiple testing correction")

    optional.add_argument("--log", dest='log',
                          action="store",
                          type=str,
                          default="script.log",
                          help="Log file to store the logging information")

    optional.add_argument("-h", "--help", action="help", default=argparse.SUPPRESS,
                          help="Show this help message and exit")

    return parser.parse_args()

def extract_polyA_sites(bam_file, fasta_file, group):
    logging.info(f"Processing file: {bam_file} as {group}")
    bam = pysam.AlignmentFile(bam_file, "rb")
    fasta = pysam.FastaFile(fasta_file)
    
    polyA_sites = []

    for read in bam.fetch():
        if not read.is_unmapped:
            seq = read.query_sequence
            match = re.search(r'(A{10,})$', seq)
            
            if match:
                read_name = read.query_name
                polyA_start = read.reference_start + match.start()
                polyA_length = len(match.group(0))
                transcript_id = read.reference_name
                chrom = read.reference_name
                coordinate = read.reference_start + match.start()
                
                if coordinate >= 20:
                    region_start = coordinate - 20
                    region_end = coordinate
                    pre_polyA_seq = fasta.fetch(chrom, region_start, region_end)
                else:
                    pre_polyA_seq = fasta.fetch(chrom, 0, coordinate)
                
                polyA_sites.append([read_name, transcript_id, coordinate, polyA_start, polyA_length, pre_polyA_seq])
    
    logging
