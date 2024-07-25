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
    
    return polyA_sites


def perform_statistical_analysis(polyA_data, fdr_threshold):
    results = []

    for transcript_id in set([site[1] for group in polyA_data.values() for site in group]):
        wt_sites = np.array([site[2] for site in polyA_data['WT'] if site[1] == transcript_id])
        mut_sites = np.array([site[2] for site in polyA_data['MUT'] if site[1] == transcript_id])

        if len(wt_sites) > 0 and len(mut_sites) > 0:
            w_distance = wasserstein_distance(wt_sites, mut_sites)
            u_statistic, p_value = mannwhitneyu(wt_sites, mut_sites, alternative='two-sided')
            results.append((transcript_id, w_distance, p_value))

    # Multiple testing correction
    p_values = [result[2] for result in results]
    reject, pvals_corrected, _, _ = multipletests(p_values, alpha=fdr_threshold, method='fdr_bh')

    significant_results = [result for result, reject_flag in zip(results, reject) if reject_flag]
    return significant_results


def main():
    args = get_args()
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s', filename=args.log, filemode='w')

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
    u_statistic, p_value = mannwhitneyu(wt_sites, mut_sites, alternative='two-sided')
    logging.info(f"Mann-Whitney U test p-value: {p_value}")

    # Summary statistics
    summary_stats = {
        'WT_Count': len(wt_sites),
        'MUT_Count': len(mut_sites),
        'WT_Mean': np.mean(wt_sites),
        'MUT_Mean': np.mean(mut_sites),
        'WT_Median': np.median(wt_sites),
        'MUT_Median': np.median(mut_sites),
        'WT_Std': np.std(wt_sites),
        'MUT_Std': np.std(mut_sites)
    }

    logging.info(f"Summary Statistics: {summary_stats}")

    # Perform statistical comparison of poly(A) site locations
    logging.info("Starting statistical analysis of poly(A) site locations")
    significant_results = perform_statistical_analysis(polyA_data, args.fdr)

    # Save significant results
    significant_df = pd.DataFrame(significant_results, columns=['TranscriptID', 'Wasserstein_Distance', 'Adjusted_p_value'])
    significant_df.to_csv(f"significant_{args.output}", sep='\t', index=False)
    logging.info(f"Significant transcripts with different poly(A) site locations have been saved to significant_{args.output}")

    logging.info(f"Total significant transcripts: {len(significant_results)}")


if __name__ == "__main__":
    main()

"""
# Example usage:
# python extract_polyA_sites.py --bam WT_rep1.bam WT_rep2.bam WT_rep3.bam MUT_rep1.bam MUT_rep2.bam MUT_rep3.bam --output polyA_sites.tsv --fasta reference_genome.fasta --groups WT WT WT MUT MUT MUT

# Output:
# The script will generate a TSV file containing the poly(A) sites for both WT and MUT groups. 
# It will also perform a statistical comparison of the poly(A) site locations, reporting the Wasserstein distance 
# and the p-value from the Mann-Whitney U test. Additionally, summary statistics for the poly(A) site locations 
# will be logged, including count, mean, median, and standard deviation for both WT and MUT groups.
# Transcripts with significantly different poly(A) site locations (based on the specified FDR threshold) 
# will be saved in a separate TSV file named significant_<output>.tsv.
"""
