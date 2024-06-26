import matplotlib.pyplot as plt


def plot_methylation_distribution(results_df, output_file):
    """
    Plot the frequency distribution of methylation sites in non-last exons, last exons, and UTRs.

    Parameters:
    results_df (DataFrame): DataFrame containing the methylation site annotations.
    output_file (str): Path to the output PDF file for the plot.
    """
    # Count the occurrences of each category
    category_counts = {
        'non_last_exon': len(results_df[(results_df['exon_number'] != 'UTR') & (results_df['is_last_exon'] == False)]),
        'last_exon': len(results_df[results_df['is_last_exon'] == True]),
        'UTR': len(results_df[results_df['exon_number'] == 'UTR'])
    }

    # Create a bar plot
    categories = list(category_counts.keys())
    counts = list(category_counts.values())

    plt.figure(figsize=(8, 6))
    plt.bar(categories, counts, color=['blue', 'green', 'red'])
    plt.xlabel('Category')
    plt.ylabel('Frequency')
    plt.title('Frequency Distribution of Methylation Sites')
    plt.tight_layout()

    # Save the plot to a PDF file
    plt.savefig(output_file)
    plt.close()

# Example usage
# results_df = pd.DataFrame(results)  # Assuming 'results' is the list of result dictionaries
# output_file = 'methylation_distribution.pdf'
# plot_methylation_distribution(results_df, output_file)
