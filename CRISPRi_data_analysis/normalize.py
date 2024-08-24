import pandas as pd
import pathlib
import math

def set_up_counts_df(all_counts_file: pathlib.Path, name_guide_column: str = 'design', sep: str = '\t') -> pd.DataFrame:
    """Set up a pandas dataframe for working with and normalizing counts data.
    Args:
        all_counts_file: path to counts .tsv file with column for guide design names and column for raw aligned read counts
        name_guide_column: name of the column containing guide design names. Default is 'design.'
        sep: separator for the counts .tsv file. Default is '\t.'
    """
    counts_df = pd.read_csv(all_counts_file, sep=sep) # initialize dataframe
    counts_df.set_index(name_guide_column, inplace=True) # set index to the guide design column
    counts_df = counts_df.astype(float) # correct data type to float for all values
    return counts_df


def median_ratio_normalize(counts_df: pd.DataFrame) -> None:
    """ Perform media ratio normalization method on counts dataframe. Adds normalized count columns to existing df.
    This method accounts for the factors of sequencing depth and RNA composition, and is the method incorporated into DESeq2.
    Good for gene count comparisons between samples and for DE analysis; NOT for within sample comparisons.
    Source: https://github.com/hbctraining/DGE_workshop/blob/master/lessons/02_DGE_count_normalization.md
    
    Args:
        counts_df: an indexed counts pandas dataframe with all raw counts for each guide (row) in each sample (column)
    """
    # Duplicate the counts df for reference
    counts_dup = counts_df.copy()
    
    # Initialize dictionary to store geometric means for each guide
    design_dict = {}

    # For every row (every guide) in the counts df, calculate the geometric mean of the non-zero values
    for design in counts_df.index.to_list():

        if not all(counts_df.loc[design] == 0):
            series = counts_df.loc[design][counts_df.loc[design] != 0]
            product = 1
            # 1st step in geom mean calculation: multiply counts across all samples
            for val in series:
                product = product * val
            # 2nd step: take the square root of the product
            design_dict[design] = math.sqrt(product) 

        # assign intermediate value of 1 to guides with all zero counts to avoid division error
        elif all(counts_df.loc[design] == 0):
            design_dict[design] = 1

        # 3rd step: calculate ratios, i.e. divide each count by the geometric mean across samples for that guide
        counts_dup.loc[design] = counts_dup.loc[design] / design_dict[design] 
    
    for sample in counts_df.columns:
        # 4th step: calculate size factors, i.e. find the median of the non-zero ratios for each sample
        median_ratio = counts_dup[sample][counts_dup[sample] != 0].median()
        # 5th step: divide original counts values by the size factors (i.e. the median ratios)
        counts_df[sample + '_normalized'] = counts_df[sample] / median_ratio


def total_reads_normalize(counts_df: pd.DataFrame) -> pd.DataFrame:
    """ Scale counts by total reads in each sample. Adds normalized count columns to existing df.
    This method accounts for sequencing depth only. Not recommended for DE analysis between samples.
    Source: https://github.com/hbctraining/DGE_workshop/blob/master/lessons/02_DGE_count_normalization.md
    
    Args:
        all_counts_df: an indexed counts pandas dataframe with raw counts for each guide (row) in each sample (column)
    """
    for sample in counts_df.columns:
        total_reads = counts_df[sample].sum()
        counts_df[sample + '_normalized'] = counts_df[sample] / total_reads

