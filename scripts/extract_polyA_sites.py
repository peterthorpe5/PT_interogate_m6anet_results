import pysam
import re
import pandas as pd
import argparse
import os
from scipy.stats import wasserstein_distance
import numpy as np
import logging


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
    
    return polyA_sites


def main():
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    args = get_args()

    polyA_data = {'WT': [], 'MUT': []}

    for bam_file, group in zip(args.bam, args.groups):
        polyA_sites = extract_polyA_sites(bam_file, args.fasta, group)
        polyA_data[group].extend(polyA_sites)

    # Convert results to DataFrame
    polyA_df_wt = pd.DataFrame(polyA_data['WT'], columns=['Read_Name', 'TranscriptID', 'Genomic_Coordinate', 'PolyA_Start', 'PolyA_Length', 'Pre_PolyA_Sequence'])
    polyA_df_mut = pd.DataFrame(polyA_data['MUT'], columns=['Read_Name', 'TranscriptID', 'Genomic_Coordinate', 'PolyA_Start', 'PolyA_Length', 'Pre_PolyA_Sequence'])
    
    # Save to TSV
    polyA_df_wt.to_csv(f"WT_{args.output}", sep='\t', index=False)
    polyA_df_mut.to_csv(f"MUT_{args.output}", sep='\t', index=False)
    logging.info(f"Poly(A) sites have been extracted and saved to {args.output}")

    # Perform statistical comparison of poly(A) site locations
    logging.info("Starting statistical analysis of poly(A) site locations")
    wt_sites = polyA_df_wt['Genomic_Coordinate'].values
    mut_sites = polyA_df_mut['Genomic_Coordinate'].values

    w_distance = wasserstein_distance(wt_sites, mut_sites)
    logging.info(f"Wasserstein distance between WT and MUT poly(A) sites: {w_distance}")

    # Additional statistical tests
    from scipy.stats import mannwhitneyu
    u_statistic, p_value = mannwhitneyu(wt_sites, mut_sites, alternative='two-sided')
    logging.info(f"Mann-Whitney U test p-value: {p_value}")


if __name__ == "__main__":
    main()
