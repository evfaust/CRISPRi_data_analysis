import os
import argparse
import pathlib
import math
import h5py
import multiprocessing as mp
import time
import csv
from Bio import SeqIO
from decimal import Decimal

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Counts reads aligned in preprocess.sh')
	parser.add_argument('-f', '--file', type=pathlib.Path, help='.tsv file from vsearch')
	parser.add_argument('-d', '--design', type=pathlib.Path, help='guide design library file in FASTA format')
	parser.add_argument('-o', '--output', type=str, help='name of output .tsv file')
	parser.add_argument('-i', '--id_filter', type=float, help='minimum alignment percent identity to include in count')
	parser.add_argument('-w', '--wash_control', type=str, help='name of wash control')
	#parser.add_argument('-c', '--readcounts', type=int, help='total reads against which to normalize')
# 'pWGA128_escapeControl'

	args = parser.parse_args()
	count = 0
	diagList = []
	designSet = set()
	for record in SeqIO.parse(args.design, 'fasta'):
		designSet.add(str(record.id))

	finalDataDict = {}
	for design in designSet:
		finalDataDict[design] = 0

	with open(args.file, 'r') as inCsv:
		read = csv.reader(inCsv, delimiter='\t')
		for row in read:
			matchID = float(row[2])
			cassette = row[1]

			if matchID >= args.id_filter:
				finalDataDict[cassette] += 1
			if cassette == args.wash_control:
				diagList.append(row)
			count += 1
			if count % 100000 == 0:
				print('[readCounter] %i lines processed' % count)

	with open('./' + args.output + '.tsv', 'w+') as outFile:
		print('design\tcounts', file=outFile)
		for key, value in finalDataDict.items():
			#prop = '%.5E' % Decimal(value/args.readcounts)
			print('%s\t%i' % (key, value), file=outFile)
	if len(diagList) > 0:
		with open('./diagnostics.tsv', 'w+') as diagOut:
			for i in diagList:
				print(i, file=diagOut)


exit()


