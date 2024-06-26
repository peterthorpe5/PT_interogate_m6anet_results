#!/usr/bin/env python3
#

import matplotlib.pyplot as plt
import seaborn as sns

def plot_methylation_distribution(results_df, output_file, transcript_lengths):
    """
    Plot the frequency distribution of methylation sites in non-last exons, last exons, and UTRs.

    Parameters:
    results_df (DataFrame): DataFrame containing the methylation site annotations.
    output_file (str): Path to the output PDF file for the plot.
    transcript_lengths (dict): Dictionary mapping transcript IDs to their lengths.
    """
    # Debugging: Print unique values in 'exon_number'
    # print("Unique values in 'exon_number':", results_df['exon_number'].unique())

    # Debugging: Print the DataFrame before filtering
    # print("DataFrame before filtering:")
    # print(results_df)

    # Check if the required columns exist
    if 'exon_number' not in results_df.columns or 'is_last_exon' not in results_df.columns:
        raise KeyError("'exon_number' or 'is_last_exon' column not found in the DataFrame.")

    # Count the occurrences of each category
    category_counts = {
        'non_last_exon': len(results_df[(results_df['exon_number'] != 'UTR') & (results_df['is_last_exon'] == False)]),
        'last_exon': len(results_df[results_df['is_last_exon'] == True]),
        'UTR': len(results_df[results_df['exon_number'] == 'UTR'])
    }

    # Debugging: Print category counts
    print("Category counts:", category_counts)

    # Normalize positions by transcript length
    results_df['normalised_position'] = results_df.apply(
        lambda row: row['position'] / transcript_lengths[row['transcript_id']], axis=1
    )

    # Create a bar plot
    categories = list(category_counts.keys())
    counts = list(category_counts.values())

    plt.figure(figsize=(12, 6))
    
    plt.subplot(1, 2, 1)  # Create subplot for bar plot
    plt.bar(categories, counts, color=['blue', 'green', 'red'])
    plt.xlabel('Category')
    plt.ylabel('Frequency')
    plt.title('Frequency Distribution of Methylation Sites')
    
    # Create a violin plot
    plt.subplot(1, 2, 2)  # Create subplot for violin plot
    results_df['category'] = results_df.apply(
        lambda row: 'last_exon' if row['is_last_exon'] else ('UTR' if row['exon_number'] == 'UTR' else 'non_last_exon'), axis=1
    )
    sns.violinplot(x='category', y='normalised_position', data=results_df, palette=['blue', 'green', 'red'])
    plt.xlabel('Category')
    plt.ylabel('Normalised Position')
    plt.title('Distribution of Methylation Sites (Normalised)')
    
    plt.tight_layout()

    # Save the plot to a PDF file
    plt.savefig(output_file)
    plt.close()

