import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import kruskal
import scikit_posthocs as sp
import os
import sys

# STATs analysis of m6anet output. 
# sys arg var [1] is a file list of genes to filter, if required. 

# File paths for data.site_proba.csv files
site_file_paths = {
    "vir_1": "./vir1_1_1/data.site_proba.csv",
    "vir_2": "./vir1_1_2/data.site_proba.csv",
    "vir_3": "./vir1_1_3/data.site_proba.csv",
    "vir_4": "./vir1_1_4/data.site_proba.csv",
    "VIRc_1": "./VIRc_1/data.site_proba.csv",
    "VIRc_2": "./VIRc_2/data.site_proba.csv",
    "VIRc_3": "./VIRc_3/data.site_proba.csv",
    "VIRc_4": "./VIRc_4/data.site_proba.csv"
}

# File paths for data.indiv_proba.csv files
indiv_file_paths = {
    "vir_1": "./vir1_1_1/data.indiv_proba.csv",
    "vir_2": "./vir1_1_2/data.indiv_proba.csv",
    "vir_3": "./vir1_1_3/data.indiv_proba.csv",
    "vir_4": "./vir1_1_4/data.indiv_proba.csv",
    "VIRc_1": "./VIRc_1/data.indiv_proba.csv",
    "VIRc_2": "./VIRc_2/data.indiv_proba.csv",
    "VIRc_3": "./VIRc_3/data.indiv_proba.csv",
    "VIRc_4": "./VIRc_4/data.indiv_proba.csv"
}

# Loading data
site_data_frames = {key: pd.read_csv(path) for key, path in site_file_paths.items()}
indiv_data_frames = {key: pd.read_csv(path) for key, path in indiv_file_paths.items()}

# Process a single site_proba data frame
def process_site_proba(data, condition, replicate):
    data['is_modified'] = data['probability_modified'] >= 0.9
    summary = data.groupby(['transcript_id', 'is_modified']).size().unstack(fill_value=0).reset_index()
    
    if True not in summary.columns:
        summary[True] = 0
    if False not in summary.columns:
        summary[False] = 0
    
    summary = summary.rename(columns={True: 'Modified', False: 'Non-Modified'})
    summary['total_sites'] = summary['Modified'] + summary['Non-Modified']
    summary['mod_ratio'] = summary['Modified'] / summary['total_sites']
    summary['condition'] = condition
    summary['replicate'] = replicate
    return summary



def process_indiv_proba(data, condition, replicate):
    """Process a single indiv_proba data frame"""
    data['is_modified'] = data['probability_modified'] >= 0.9
    summary = data.groupby(['transcript_id', 'is_modified']).size().unstack(fill_value=0).reset_index()
    
    if True not in summary.columns:
        summary[True] = 0
    if False not in summary.columns:
        summary[False] = 0
    
    summary = summary.rename(columns={True: 'Modified', False: 'Non-Modified'})
    summary['mod_ratio'] = summary['Modified'] / (summary['Modified'] + summary['Non-Modified'])
    summary['condition'] = condition
    summary['replicate'] = replicate
    return summary


def pairwise_dunn_test(data, group_col, value_col):
    """
    Perform pairwise Dunn's test with Benjamini-Hochberg correction.

    Parameters:
    data (pd.DataFrame): DataFrame containing the data to be tested.
    group_col (str): Name of the column containing the group labels.
    value_col (str): Name of the column containing the values to be tested.

    Returns:
    pd.DataFrame: DataFrame containing the pairwise comparison p-values.
    """
    comparisons = sp.posthoc_dunn(data, val_col=value_col, 
                                  group_col=group_col, p_adjust='fdr_bh')
    return comparisons


# Check if a gene list file is provided
if len(sys.argv) > 1:
    gene_list_file_path = sys.argv[1]
    gene_list_file_name = os.path.basename(gene_list_file_path).split('.')[0]
    filter_genes = True
else:
    gene_list_file_path = None
    gene_list_file_name = "all_genes"
    filter_genes = False

print("Looking at   -->  ", gene_list_file_path)


def read_gene_list(file_path):
    """
    Read a list of genes from a file, one gene per line.

    Parameters:
    file_path (str): Path to the file containing the list of genes.

    Returns:
    set: Set of genes read from the file.
    """
    gene_set = set()
    with open(file_path, 'r') as file:
        genes = file.read().splitlines()
    for gene in genes:
        if not gene.endswith(".1"):
            gene = gene + ".1"
        gene_set.add(gene)
    return gene_set

# Read the list of genes
if filter_genes:
    gene_list = read_gene_list(gene_list_file_path)

# Process all site_proba data frames
site_proba_data = pd.concat([
    process_site_proba(site_data_frames['vir_1'], 'vir1', 1),
    process_site_proba(site_data_frames['vir_2'], 'vir1', 2),
    process_site_proba(site_data_frames['vir_3'], 'vir1', 3),
    process_site_proba(site_data_frames['vir_4'], 'vir1', 4),
    process_site_proba(site_data_frames['VIRc_1'], 'VIRc', 1),
    process_site_proba(site_data_frames['VIRc_2'], 'VIRc', 2),
    process_site_proba(site_data_frames['VIRc_3'], 'VIRc', 3),
    process_site_proba(site_data_frames['VIRc_4'], 'VIRc', 4)
])

# Process all indiv_proba data frames
indiv_proba_data = pd.concat([
    process_indiv_proba(indiv_data_frames['vir_1'], 'vir1', 1),
    process_indiv_proba(indiv_data_frames['vir_2'], 'vir1', 2),
    process_indiv_proba(indiv_data_frames['vir_3'], 'vir1', 3),
    process_indiv_proba(indiv_data_frames['vir_4'], 'vir1', 4),
    process_indiv_proba(indiv_data_frames['VIRc_1'], 'VIRc', 1),
    process_indiv_proba(indiv_data_frames['VIRc_2'], 'VIRc', 2),
    process_indiv_proba(indiv_data_frames['VIRc_3'], 'VIRc', 3),
    process_indiv_proba(indiv_data_frames['VIRc_4'], 'VIRc', 4)
])

# Filter the data to include only the genes in the provided gene list
if filter_genes:
    site_proba_data = site_proba_data[site_proba_data['transcript_id'].isin(gene_list)]
    indiv_proba_data = indiv_proba_data[indiv_proba_data['transcript_id'].isin(gene_list)]



# Find common genes across all conditions
print(" I am filtering this for genes that are common in both datsets")
common_genes = set(site_proba_data['transcript_id'])
for condition in site_proba_data['condition'].unique():
    condition_genes = set(site_proba_data[site_proba_data['condition'] == condition]['transcript_id'])
    common_genes = common_genes.intersection(condition_genes)

print("there are %d common genes\n" % (len(condition_genes)))
gene_list_file_name = "common_%d_genes" % (len(condition_genes))

# Filter the data to include only common genes
site_proba_data = site_proba_data[site_proba_data['transcript_id'].isin(common_genes)]
indiv_proba_data = indiv_proba_data[indiv_proba_data['transcript_id'].isin(common_genes)]

#

# Ensure the order of columns is consistent
site_proba_data = site_proba_data[['transcript_id', 'Non-Modified', 'Modified', 'total_sites', 'mod_ratio', 'condition', 'replicate']]
indiv_proba_data = indiv_proba_data[['transcript_id', 'Modified', 'Non-Modified', 'mod_ratio', 'condition', 'replicate']]

# Write the first "number_of_lines" lines of each DataFrame to a file with tab-separated values
number_of_lines = 100000
site_proba_data.head(number_of_lines).to_csv(f'{gene_list_file_name}_site_proba_data_first_{number_of_lines}_lines.csv', index=False, sep='\t')
indiv_proba_data.head(number_of_lines).to_csv(f'{gene_list_file_name}_indiv_proba_data_first_{number_of_lines}_lines.csv', index=False, sep='\t')

#########################################################################
# Generate the box plot for site proba data with violin plot for density
plt.figure(figsize=(12, 8))
sns.violinplot(x='condition', y='mod_ratio', data=site_proba_data, inner=None, palette="vlag")
sns.boxplot(x='condition', y='mod_ratio', data=site_proba_data, whis=[0, 100], width=.2, palette="vlag")
sns.stripplot(x='condition', y='mod_ratio', data=site_proba_data, size=1, color=".4", linewidth=0)
plt.title(f'Ratio of Modified to Non-Modified Sites per Transcript per Experiment (Site Proba) - {gene_list_file_name}')
plt.ylabel('Ratio of Modified to Non-Modified Sites')
plt.xlabel('Condition')
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig(f'{gene_list_file_name}_ratio_of_modified_sites_boxplot_site_proba.png')
plt.show()


#########################################################################
# Generate the box plot for indiv proba data with violin plot for density
plt.figure(figsize=(12, 8))
sns.violinplot(x='condition', y='mod_ratio', data=indiv_proba_data, inner=None, palette="vlag")
sns.boxplot(x='condition', y='mod_ratio', data=indiv_proba_data, whis=[0, 100], width=.2, palette="vlag")
sns.stripplot(x='condition', y='mod_ratio', data=indiv_proba_data, size=1, color=".4", linewidth=0)
plt.title(f'Ratio of Modified to Non-Modified Reads per Transcript per Experiment (Indiv Proba) - {gene_list_file_name}')
plt.ylabel('Ratio of Modified to Non-Modified Reads')
plt.xlabel('Condition')
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig(f'{gene_list_file_name}_ratio_of_modified_reads_boxplot_indiv_proba.png')
plt.show()


#########################################################################
# Calculate mean modification ratio per site per transcript
mean_mod_ratio_data = site_proba_data.groupby(['transcript_id', 'condition']).agg(
    mean_mod_ratio=pd.NamedAgg(column='mod_ratio', aggfunc='mean')
).reset_index()


#########################################################################
# Generate the box plot for mean modification ratio per site per transcript
plt.figure(figsize=(12, 8))
sns.violinplot(x='condition', y='mean_mod_ratio', data=mean_mod_ratio_data, inner=None, palette="vlag")
sns.boxplot(x='condition', y='mean_mod_ratio', data=mean_mod_ratio_data, whis=[0, 100], width=.2, palette="vlag")
sns.stripplot(x='condition', y='mean_mod_ratio', data=mean_mod_ratio_data, size=1, color=".4", linewidth=0)
plt.title(f'Mean Modification Ratio per Site per Transcript per Experiment - {gene_list_file_name}')
plt.ylabel('Mean Modification Ratio per Site')
plt.xlabel('Condition')
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig(f'{gene_list_file_name}_mean_modification_ratio_per_site_boxplot.png')
plt.show()

print("""Summary of Comparisons:
Modification Ratios (Site Proba): Comparison of modification ratios of modified to non-modified sites across conditions.
Modification Ratios (Indiv Proba): Comparison of modification ratios of modified to non-modified reads across conditions.
Mean Modification Ratios (Per Site Per Transcript): Comparison of mean modification ratios per site per transcript across conditions.
Number of Modified Reads: Comparison of the number of modified reads per transcript across conditions\n\n""")

# Summary statistics for site proba
summary_stats_site = site_proba_data.groupby('condition')['mod_ratio'].describe()
print(summary_stats_site)

# Summary statistics for indiv proba
summary_stats_indiv = indiv_proba_data.groupby('condition')['mod_ratio'].describe()
print(summary_stats_indiv)

# Summary statistics for mean modification ratio per site per transcript
summary_stats_mean_mod_ratio = mean_mod_ratio_data.groupby('condition')['mean_mod_ratio'].describe()
print(summary_stats_mean_mod_ratio)

# Kruskal-Wallis test for site proba
kruskal_test_site = kruskal(*[group["mod_ratio"].values for name, group in site_proba_data.groupby("condition")])
print("\nComparison of modified sites")
print("This test compares the modification ratios (mod_ratio) of modified to non-modified sites across different conditions")
print(f"Kruskal-Wallis test for site proba: {kruskal_test_site}")

# Kruskal-Wallis test for indiv proba
kruskal_test_indiv = kruskal(*[group["mod_ratio"].values for name, group in indiv_proba_data.groupby("condition")])
print("\nComparison of modified reads")
print("This test compares the modification ratios (mod_ratio) of modified to non-modified reads across different conditions.")
print(f"Kruskal-Wallis test for indiv proba: {kruskal_test_indiv}")

# Kruskal-Wallis test for mean modification ratio per site per transcript
kruskal_test_mean_mod_ratio = kruskal(*[group["mean_mod_ratio"].values for name, group in mean_mod_ratio_data.groupby("condition")])
print("\nMean modification rate")
print("This test compares the mean modification ratios per site per transcript across different conditions.")
print(f"Kruskal-Wallis test for mean modification ratio per site per transcript: {kruskal_test_mean_mod_ratio}")

# Kruskal-Wallis test for the number of modified to non-modified reads per transcript
kruskal_test_reads = kruskal(*[group["Modified"].values for name, group in indiv_proba_data.groupby("condition")])
print("\nNumber of modified reads")
print("This test compares the number of modified reads per transcript across different conditions.")
print(f"Kruskal-Wallis test for the number of modified reads per transcript: {kruskal_test_reads}")

# Post hoc pairwise comparisons using Dunn's test with Benjamini-Hochberg correction
pairwise_comparisons_site = pairwise_dunn_test(site_proba_data, 'condition', 'mod_ratio')
pairwise_comparisons_indiv = pairwise_dunn_test(indiv_proba_data, 'condition', 'mod_ratio')
pairwise_comparisons_mean_mod_ratio = pairwise_dunn_test(mean_mod_ratio_data, 'condition', 'mean_mod_ratio')
pairwise_comparisons_reads = pairwise_dunn_test(indiv_proba_data, 'condition', 'Modified')

print("\n\t### SECTION: Dunn test with BH post hoc correction ###\n\n")

print("\tPairwise comparisons for site proba:")
print("Modification ratios of modified to non-modified SITES across conditions.")
print(pairwise_comparisons_site)

print("\nPairwise comparisons for indiv proba:")
print("Modification ratios of modified to non-modified READS across conditions.")
print(pairwise_comparisons_indiv)

print("\nPairwise comparisons for mean modification ratio per site per transcript:")
print("Mean modification ratios per site per transcript across conditions.")
print(pairwise_comparisons_mean_mod_ratio)

print("\nPairwise comparisons for the number of modified reads per transcript:")
print("Number of modified reads per transcript across conditions.")
print(pairwise_comparisons_reads)

# Save summary statistics and comparison results to TSV files
summary_stats_site.to_csv(f'{gene_list_file_name}_summary_statistics_site.tsv', sep='\t')
summary_stats_indiv.to_csv(f'{gene_list_file_name}_summary_statistics_indiv.tsv', sep='\t')
summary_stats_mean_mod_ratio.to_csv(f'{gene_list_file_name}_summary_statistics_mean_mod_ratio.tsv', sep='\t')
pairwise_comparisons_site.to_csv(f'{gene_list_file_name}_pairwise_comparisons_site.tsv', sep='\t')
pairwise_comparisons_indiv.to_csv(f'{gene_list_file_name}_pairwise_comparisons_indiv.tsv', sep='\t')
pairwise_comparisons_mean_mod_ratio.to_csv(f'{gene_list_file_name}_pairwise_comparisons_mean_mod_ratio.tsv', sep='\t')
pairwise_comparisons_reads.to_csv(f'{gene_list_file_name}_pairwise_comparisons_reads.tsv', sep='\t')
