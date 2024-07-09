import pandas as pd
import os
import argparse
import pathlib

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Combines count files from counter.py')
    parser.add_argument('-c', '--count_files_folder', type=pathlib.Path, help='count files folder')


    args = parser.parse_args()

    # Initialize an empty DataFrame
    merged_df = pd.DataFrame(columns=['design'])
 
    # Get all files in the count_files_folder
    files = os.listdir(str(args.count_files_folder))

    # Merge all the .tsv files
    for file in files:
        if file.endswith('.tsv'):
            file_path = os.path.join(str(args.count_files_folder), file)
            df = pd.read_csv(file_path, sep='\t')
            df.rename(columns={'counts': file.replace('_counts.tsv','')}, inplace=True)
            merged_df = pd.merge(merged_df, df, on='design', how='outer')

    merged_df.fillna(0, inplace=True)  # Replace NAs with 0s

    # Save the merged DataFrame to a new file
    output_file = os.path.join(str(args.count_files_folder), 'all_counts.tsv')
    merged_df.to_csv(output_file, sep='\t', index=False)

