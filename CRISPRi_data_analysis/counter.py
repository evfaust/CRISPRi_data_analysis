# Adapted from William "Bill" Alexander's code

from Bio import SeqIO
import pathlib
import csv
import pandas as pd
import os

def count_abundance(file: pathlib.Path, design: pathlib.Path, out_file: str, wash_control: str, id_filter: float = 100.0) -> None:
	"""Count the number of reads aligned to each guide in the design library. Also counts the wash control for diagnostics."
	Args:
		file: .tsv alignment file from vsearch
		design: guide design library file in FASTA format
		output_folder: folder to save output files
		name_out_file: name of output count .tsv file
		wash_control: name of wash control
		id_filter: minimum alignment percent identity to include in count. Default is 100%

	"""
	#count = 0
	wash_count = 0

	designSet = set()
	for record in SeqIO.parse(design, 'fasta'):
		designSet.add(str(record.id))

	finalDataDict = {}
	for design in designSet:
		finalDataDict[design] = 0

	with open(file, 'r') as inCsv:
		read = csv.reader(inCsv, delimiter='\t')
		for row in read:
			matchID = float(row[2])
			cassette = row[1]

			if matchID >= id_filter:
				finalDataDict[cassette] += 1
			if cassette == wash_control:
				wash_count += 1
			#count += 1
			#if count % 100000 == 0:
				#print('[readCounter] %i lines processed' % count)

	with open(out_file, 'w+') as outFile:
		print('design\tcounts', file=outFile)
		for key, value in finalDataDict.items():
			#prop = '%.5E' % Decimal(value/args.readcounts)
			print('%s\t%i' % (key, value), file=outFile)
	if wash_count > 0:
		print(f'Warning: wash control present - count: {wash_count} in {file}')


def combine_count_files(count_files_folder: pathlib.Path, out_all_counts_file: str) -> None:
	"""Combine all count files in a folder into a single .tsv file.
	Args:
		count_files_folder: folder containing count files. Count files must end in '_counts.tsv'
		out_all_counts_file: name of output file (example: 'all_counts.tsv')
	"""
	
	# Initialize an empty DataFrame
	merged_df = pd.DataFrame(columns=['design'])
 
	# Get all files in the count_files_folder
	files = os.listdir(count_files_folder)

	# Merge all the .tsv files
	for file in files:
		if file.endswith('_counts.tsv'):
			file_path = os.path.join(count_files_folder, file)
			df = pd.read_csv(file_path, sep='\t')
			df.rename(columns={'counts': file.replace('_counts.tsv','')}, inplace=True)
			merged_df = pd.merge(merged_df, df, on='design', how='outer')

	# Replace NAs with 0s
	merged_df.fillna(0, inplace=True)

	# Save the merged DataFrame to a new file
	output_file = os.path.join(str(count_files_folder), out_all_counts_file)
	merged_df.to_csv(output_file, sep='\t', index=False)

#combine_count_files("/Users/evelynfaust/Desktop/ORNL/CRISPRi_data_analysis/NEW_TEST_ALIGN_OUT", "all_counts.tsv")

