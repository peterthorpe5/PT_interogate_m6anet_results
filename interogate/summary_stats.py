import pandas as pd
from scipy.stats import chi2_contingency


# The chi-squared statistic is a measure of the difference between the observed and expected frequencies in your data.
# A chi-squared value of 0.0 indicates that there is no difference between the observed and expected frequencies.
# P-Value (p_value):
# The p-value indicates the probability that the observed difference (or a more extreme one) would occur if the null hypothesis were true.
# A p-value of 1.0 indicates that there is no evidence to reject the null hypothesis, 


def summarize_methylation_sites(results_df, output_file, logger):
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
        'utr_sites': len(df[df['exon_number'] == 'UTR'])
    })).reset_index()

    # Calculate ratios for each transcript
    summary['ratio_non_last_to_last'] = summary.apply(
        lambda row: row['non_last_exon_sites'] / row['last_exon_sites'] if row['last_exon_sites'] != 0 else float('inf'),
        axis=1
    )

    # Calculate overall counts
    total_sites = summary['total_sites'].sum()
    non_last_exon_sites = summary['non_last_exon_sites'].sum()
    last_exon_sites = summary['last_exon_sites'].sum()
    utr_sites = summary['utr_sites'].sum()


    # chi-squared test can be used to compare the observed distribution 
    # to an expected distribution. If we don't assume equal distribution,
    # the expected frequencies should be based on the proportion of each category 
    # observed in the data
    # Calculate the expected frequencies based on the overall distribution
    total_observed = non_last_exon_sites + last_exon_sites + utr_sites
    expected_non_last_exon = total_observed * (non_last_exon_sites / total_sites)
    expected_last_exon = total_observed * (last_exon_sites / total_sites)
    expected_utr = total_observed * (utr_sites / total_sites)

    observed = [non_last_exon_sites, last_exon_sites, utr_sites]
    expected = [expected_non_last_exon, expected_last_exon, expected_utr]

    chi2, p, _, _ = chi2_contingency([observed, expected])

    # Add statistical test results to the summary
    summary.loc['Total'] = summary.sum(numeric_only=True)
    summary.at['Total', 'transcript_id'] = 'Overall'
    summary['chi2'] = chi2
    summary['p_value'] = p

    # Write summary to a file
    summary.to_csv(output_file, index=False, sep="\t")
    out_note = f"Summary saved to {output_file}"
    logger.info(out_note)

    # Print the summary DataFrame for visual confirmation
    print(summary)