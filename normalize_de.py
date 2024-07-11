import pandas as pd
import argparse
import pathlib
import numpy as np
from scipy.stats import ttest_ind
import math


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Normalizes counts, calculates fold change and p-values')
    parser.add_argument('-i', '--all_counts_file', type=pathlib.Path, help='counts .tsv file, must have design column')
    parser.add_argument('-n', '--normalization_method', type=str, help='method to normalize counts. Options: m = median, t = total reads')
    parser.add_argument('-p', '--sample_prefixes', type=str, help='prefix of sample (col) names in counts .tsv file (must be the same among replicates)')
    parser.add_argument('-c', '--control_condition_prefix', type=str, help='prefix of control condition (col) names in counts .tsv file')
    parser.add_argument('-o', '--output', type=str, help='name of output .tsv file')


    args = parser.parse_args()

    # Set up the counts DataFrame
    counts_df = pd.read_csv(args.all_counts_file, sep='\t')
    counts_df.set_index('design', inplace=True) 
    counts_df = counts_df.astype(int)

    # Normalize counts
    if args.normalization_method == 'm':
        """ Populate dictionary with geometric means for each design"""
        design_dict = {}
        for design in counts_df.index.to_list():
            if not all(counts_df.loc[design] == 0):
                series = counts_df.loc[design][counts_df.loc[design] != 0]
                product = 1
                for val in series:
                    product = product * val
                design_dict[design] = math.sqrt(product)
            elif all(counts_df.loc[design] == 0):
                design_dict[design] = 1
        """ Done """
        for sample in counts_df.columns:
            for design in counts_df.index.to_list():
                counts_df[design,'temp_' + sample] = counts_df.loc[design,sample] / design_dict[design]
            median_ratio = counts_df['temp_' + sample][counts_df['temp_' + sample] != 0].median()
            counts_df[design, sample + '_normalized'] = counts_df.loc[design,sample] / median_ratio

    elif args.normalization_method == 't':
        for sample in counts_df.columns:
            total_reads = counts_df[sample].sum()
            counts_df[sample + '_normalized'] = counts_df[sample] / total_reads


    # Calculate fold change and p-values
    ctr_colnames = [col for col in counts_df.columns if col.startswith(args.control_condition_prefix) and col.endswith('_normalized')]
    for gene in counts_df.index.to_list():
        for c_col in ctr_colnames:
            average_ctr = counts_df.loc[gene, c_col].mean()
            for sample in args.sample_prefixes.split(','):
                sample_colnames = [col for col in counts_df.columns if col.startswith(sample) and col.endswith('_normalized')]
                average_sample = counts_df.loc[gene, sample_colnames].mean()
                counts_df[gene, sample + '_log2FC'] = np.log2(average_sample / average_ctr)
                t_stat, p_val = ttest_ind(counts_df.loc[gene, sample_colnames], counts_df.loc[gene, ctr_colnames], equal_var=False)
                counts_df[gene, sample + '_pval'] = p_val
