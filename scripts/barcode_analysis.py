import csv
import gzip
import os
from Bio import SeqIO
from collections import Counter


def frequency(file):
	"""
	Gets counts of each unique read in fastq file.

	Inputs: 
		- file: fastq.gz file

	Output: returns dictionary in the following format: {read: count}
	"""

	reads = []
	with gzip.open(file, 'rt') as filename:
		for rec in SeqIO.parse(filename, 'fastq'):
			reads.append(rec.seq)
	return dict(Counter(reads))


def check_unique_barcodes(data, barcodesfile, region_num, newfile):
	"""
	Gets individual barcode counts for barcodes in specific region region_num.

	Inputs:
		- data: dictionary representing read counts (output of frequency())
		- barcodesfile: path to fasta file containing barcode metadata
		- region_num: barcode region to check barcode counts (0, 1, or 2)
		- newfile: name of output csv file

	Output: csv file with individual barcode counts
	"""

	unique_barcodes = {}
	barcodes = [[] for _ in range(3)]
	for rec in SeqIO.parse(barcodesfile, 'fasta'):
		region = int(rec.name.partition('region')[2][0])
		barcodes[region-1].append(str(rec.seq))

	for barcode in barcodes[region_num]:
		for key, count in data.items():
			sequence = str(key)
			if barcode in sequence:
				if barcode in unique_barcodes.keys():
					unique_barcodes[barcode] += count
				else:
					unique_barcodes[barcode] = count
	
	total_counts = sum(unique_barcodes.values())

	with open(str(newfile), 'w') as csv_file:  
		writer = csv.writer(csv_file)
		writer.writerow(['Barcode', 'Count', 'Frequency'])
		for key, value in unique_barcodes.items():
			writer.writerow([key, value, (value/total_counts)])
	

def check_combinations(data, barcodesfile, newfile):
	"""
	Gets counts of each unique chimeric barcode sequence (BC1 + BC2 + BC3) present in data.
	Multiple barcodes have to be present in a single read.

	Inputs:
		- data: dictionary representing read counts (output of frequency())
		- barcodesfile: path to fasta file containing barcode metadata
		- newfile: name of output csv file

	Output: csv file with each row being a unique chimeric barcode sequence corresponding to its count
	"""

	barcodes = [[] for _ in range(3)]
	for rec in SeqIO.parse(barcodesfile, 'fasta'):
		region = int(rec.name.partition('region')[2][0])
		barcodes[region-1].append(str(rec.seq))

	barcode_combos = {}
	for key, count in data.items():
		combination = []
		sequence = str(key)
		one_barcode = False
		for region in barcodes:
			for barcode in region:
				if barcode in sequence:
					combination.append(barcode)
					one_barcode = True
					break
				elif barcode == region[-1]:
					combination.append('')
					break
				else:
					pass
		seq_str = '_'.join(combination)
		if one_barcode and '' not in seq_str.split('_'):
			if seq_str in barcode_combos.keys():
				barcode_combos[seq_str] += count
			else:
				barcode_combos[seq_str] = count

	
	total_counts = sum(barcode_combos.values())
	
	with open(str(newfile), 'w') as csv_file:  
		writer = csv.writer(csv_file)
		writer.writerow(['BC1', 'BC2', 'BC3', 'Count', 'Frequency'])
		for seq, count in barcode_combos.items():
			barcodes = seq.split('_')
			writer.writerow([barcodes[0], barcodes[1], barcodes[2], count, (count/total_counts)])



if __name__ == "__main__":
	barcodes_file = 'metadata/barcodes.fasta'

	filename = os.path.basename(snakemake.input[0])
	reads_file = snakemake.input[0]

	final_reads = frequency(reads_file) # calculates frequencies of reads

	for i in range(len(snakemake.output)-1): # finds individual counts of barcodes in all 3 regions
		check_unique_barcodes(final_reads, barcodes_file, i, snakemake.output[i])

	check_combinations(final_reads, barcodes_file, snakemake.output[-1]) # finds counts of chimeric barcode sequencies