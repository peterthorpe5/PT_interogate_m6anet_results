#!/usr/bin/env python3
#

import pandas as pd
from scipy.stats import chi2_contingency
import numpy as np


    # chi-squared test can be used to compare the observed distribution 
    # to an expected distribution. If we don't assume equal distribution,
    # the expected frequencies should be based on the proportion of each category 
    # observed in the data
    # Calculate the expected frequencies based on the overall distribution

# The chi-squared statistic is a measure of the difference between the observed and expected frequencies in your data.
# A chi-squared value of 0.0 indicates that there is no difference between the observed and expected frequencies.
# P-Value (p_value):
# The p-value indicates the probability that the observed difference (or a more extreme one) would occur if the null hypothesis were true.
# A p-value of 1.0 indicates that there is no evidence to reject the null hypothesis, 


import pandas as pd
from scipy.stats import chi2_contingency
import numpy as np

def benjamini_hochberg(p_values):
    """
    Apply the Benjamini-Hochberg correction for multiple hypothesis testing.

    Parameters:
    p_values (list): List of p-values to correct.

    Returns:
    list: List of adjusted p-values.
    """
    p_values = np.array(p_values)
    n = len(p_values)
    sorted_indices = np.argsort(p_values)
    sorted_p_values = p_values[sorted_indices]
    adjusted_p_values = np.zeros(n)
    
    cummin = sorted_p_values[-1]
    adjusted_p_values[sorted_indices[-1]] = cummin
    for i in range(n-2, -1, -1):
        cummin = min(cummin, sorted_p_values[i] * n / (i + 1))
        adjusted_p_values[sorted_indices[i]] = cummin
    
    return adjusted_p_values

def summarise_methylation_sites(results_df, output_file, logger):
    """
    Summarize the number of methylation sites per transcript and perform statistical comparison.

    Parameters:
    results_df (DataFrame): DataFrame containing the methylation site annotations.
    output_file (str): Path to the output file for the summary.
    logger (Logger): Logger for logging information.
    """
    # Summarize the data
    summary = results_df.groupby('transcript_id').apply(lambda df: pd.Series({
        'total_sites': len(df),
        'non_last_exon_sites': len(df[(df['exon_number'] != 'UTR') & (df['is_last_exon'] == False)]),
        'last_exon_sites': len(df[df['is_last_exon'] == True]),
        'utr_sites': len(df[df['exon_number'] == 'UTR']),
        'last_exon_and_utr_sites': len(df[(df['exon_number'] == 'UTR') | (df['is_last_exon'] == True)])
    })).reset_index()

    # Calculate ratios for each transcript
    summary['ratio_non_last_to_last_and_utr'] = summary.apply(
        lambda row: row['non_last_exon_sites'] / row['last_exon_and_utr_sites'] if row['last_exon_and_utr_sites'] != 0 else 0,
        axis=1
    )

    # Initialize columns for chi-squared test results
    summary['chi2'] = np.nan
    summary['p_value'] = np.nan

    # Perform chi-squared test for each transcript
    for index, row in summary.iterrows():
        non_last_exon_sites = row['non_last_exon_sites']
        last_exon_and_utr_sites = row['last_exon_and_utr_sites']
        total_sites = row['total_sites']
        if total_sites == 0:
            continue

        # Define observed frequencies
        observed = [non_last_exon_sites, last_exon_and_utr_sites]

        # Calculate the expected frequencies based on the overall distribution
        expected_non_last_exon = total_sites * (non_last_exon_sites / total_sites)
        expected_last_exon_and_utr = total_sites * (last_exon_and_utr_sites / total_sites)
        
        expected = [expected_non_last_exon, expected_last_exon_and_utr]

        # Check if expected frequencies have zero elements
        if any(e == 0 for e in expected):
            continue
        
        try:
            chi2, p, _, _ = chi2_contingency([observed, expected])
            summary.at[index, 'chi2'] = chi2
            summary.at[index, 'p_value'] = p
        except ValueError as e:
            logger.warning(f"Chi-squared test could not be performed for transcript {row['transcript_id']}: {e}")

    # Apply Benjamini-Hochberg correction to the p-values
    p_values = summary['p_value'].dropna().tolist()
    if p_values:
        adjusted_p_values = benjamini_hochberg(p_values)
        summary.loc[summary['p_value'].notna(), 'adjusted_p_value'] = adjusted_p_values

    # Write summary to a file
    summary.to_csv(output_file, index=False, sep="\t")
    out_note = f"Summary saved to {output_file}"
    logger.info(out_note)

    # Print the summary DataFrame for visual confirmation
    print(summary)

    # Generate overall summary statistics
    overall_summary = {
        'total_transcripts': len(summary),
        'total_sites': summary['total_sites'].sum(),
        'mean_sites_per_transcript': summary['total_sites'].mean(),
        'median_sites_per_transcript': summary['total_sites'].median(),
        'std_sites_per_transcript': summary['total_sites'].std(),
        'max_sites_per_transcript': summary['total_sites'].max(),
        'mean_ratio_non_last_to_last_and_utr': summary['ratio_non_last_to_last_and_utr'].mean(),
        'median_ratio_non_last_to_last_and_utr': summary['ratio_non_last_to_last_and_utr'].median(),
        'std_ratio_non_last_to_last_and_utr': summary['ratio_non_last_to_last_and_utr'].std(),
        'total_non_last_exon_sites': summary['non_last_exon_sites'].sum(),
        'mean_non_last_exon_sites': summary['non_last_exon_sites'].mean(),
        'median_non_last_exon_sites': summary['non_last_exon_sites'].median(),
        'std_non_last_exon_sites': summary['non_last_exon_sites'].std(),
        'max_non_last_exon_sites': summary['non_last_exon_sites'].max(),
        'total_last_exon_sites': summary['last_exon_sites'].sum(),
        'mean_last_exon_sites': summary['last_exon_sites'].mean(),
        'median_last_exon_sites': summary['last_exon_sites'].median(),
        'std_last_exon_sites': summary['last_exon_sites'].std(),
        'max_last_exon_sites': summary['last_exon_sites'].max(),
        'total_utr_sites': summary['utr_sites'].sum(),
        'mean_utr_sites': summary['utr_sites'].mean(),
        'median_utr_sites': summary['utr_sites'].median(),
        'std_utr_sites': summary['utr_sites'].std(),
        'max_utr_sites': summary['utr_sites'].max()
    }

    overall_summary_df = pd.DataFrame(overall_summary.items(), columns=['Statistic', 'Value'])

    output_summary_file = output_file + ".overall.summary"
    overall_summary_df.to_csv(output_summary_file, index=False, sep="\t")
    logger.info(f"Overall summary saved to {output_summary_file}")
    logger.info(overall_summary_df)

# Example usage
# results_df = pd.DataFrame(results)  # Assuming 'results' is the list of result dictionaries
# output_summary = 'methylation_summary.txt'
# summarize_methylation_sites(results_df, output_summary, logger)
