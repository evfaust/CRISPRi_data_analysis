import pandas as pd
import glob

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Combines count files from counter.py')
    parser.add_argument('-c', '--count_files_folder', type=pathlib.Path, help='count files folder')


    args = parser.parse_args()

    # Get all .tsv files in the folder
    file_pattern = str(args.file / '*.tsv')
    file_list = glob.glob(file_pattern)

    # Initialize an empty DataFrame
    merged_df = pd.DataFrame()

    # Merge all the .tsv files
    for file in file_list:
        df = pd.read_csv(file, sep='\t')
        merged_df = pd.merge(merged_df, df, on='design', how='outer')

    merged_df.fillna(0, inplace=True)  # Replace NAs with 0s
    

    # Save the merged DataFrame to a new file
    merged_df.to_csv(args.count_files_folder/'merged_counts.tsv', sep='\t', index=False)