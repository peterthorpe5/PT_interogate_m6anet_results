import pandas as pd
from scipy.stats import chi2_contingency


def summarize_methylation_sites(results_df, output_file):
    """
    Summarize the number of methylation sites per transcript and perform statistical comparison.

    Parameters:
    results_df (DataFrame): DataFrame containing the methylation site annotations.
    output_file (str): Path to the output file for the summary.
    """
    # Summarize the data
    summary = results_df.groupby('transcript_id').apply(lambda df: pd.Series({
        'total_sites': len(df),
        'non_last_exon_sites': len(df[(df['exon_number'] != 'UTR') & (df['is_last_exon'] == False)]),
        'last_exon_sites': len(df[df['is_last_exon'] == True]),
        'utr_sites': len(df[df['exon_number'] == 'UTR'])
    })).reset_index()

    # Statistical comparison: Chi-squared test
    total_sites = summary['total_sites'].sum()
    non_last_exon_sites = summary['non_last_exon_sites'].sum()
    last_exon_sites = summary['last_exon_sites'].sum()
    utr_sites = summary['utr_sites'].sum()

    contingency_table = [
        [non_last_exon_sites, last_exon_sites, utr_sites],
        [total_sites - non_last_exon_sites, total_sites - last_exon_sites, total_sites - utr_sites]
    ]

    chi2, p, _, _ = chi2_contingency(contingency_table)

    # Add statistical test results to the summary
    summary.loc['Total'] = summary.sum(numeric_only=True)
    summary.at['Total', 'transcript_id'] = 'Overall'
    summary['chi2'] = chi2
    summary['p_value'] = p

    # Write summary to a file
    summary.to_csv(output_file, index=False, sep="\t")
    print(f"Summary saved to {output_file}")


