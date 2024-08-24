import pandas as pd
import numpy as np
from scipy.stats import ttest_ind


def compare_abundance(counts_df: pd.DataFrame, sample_prefixes: list[str], control_prefix: str) -> None:
    """Calculate log2FC values between control and experimental conditions (combining replicates) and perform 2-tailed heteroscedastic t-tests.
    
    Args:
        counts_df: an indexed pandas dataframe with normalized counts for each guide (row) in each sample (column)
        sample_prefixes: list of prefixes for sample column names in counts_df
        control_prefix: prefix for control condition column names in counts_df
    """
    
    # Identify control condition columns
    ctr_colnames = [col for col in counts_df.columns 
                    if col.startswith(control_prefix) and col.endswith('_normalized')] # TODO: Incorporate option for using non-normalized counts (right now the column name has to end with _normalized)
    
    for design in counts_df.index.to_list():
        
        # For each design, calculate average counts from control condition replicates
        average_ctr = float(counts_df.loc[design, ctr_colnames].mean())

        # Compare every sample to each control condition
        for sample in sample_prefixes:
            
            # Identify experimental condition columns
            sample_colnames = [col for col in counts_df.columns 
                               if col.startswith(sample) and col.endswith('_normalized') 
                               and not (col.endswith('_log2FC') or col.endswith('_pval'))]
            
            # For this sample and this guide, calculate average counts from experimental condition replicates 
            average_sample = float(counts_df.loc[design, sample_colnames].mean())

            # For this guide, compare average experimental condition counts to average control condition counts
            # Calculate Log2FC: 
            counts_df.loc[design, sample + '_log2FC'] = np.log2(average_sample / average_ctr)
            # Perform two-sided heteroscedastic t-test:
            t_stat, p_val = ttest_ind(counts_df.loc[design, sample_colnames], 
                                      counts_df.loc[design, ctr_colnames], 
                                      equal_var=False)
            counts_df.loc[design, sample + '_pval'] = p_val


