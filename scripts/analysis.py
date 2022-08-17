import gzip
import os
from Bio import SeqIO
from collections import Counter


def frequency(file):
	reads = []
	print(f"Loading sequences: {os.path.basename(file)}")
	with gzip.open(file, 'rt') as filename:
		for rec in SeqIO.parse(filename, 'fastq'):
			reads.append(rec.seq)
	return dict(Counter(reads))


def check_unique_barcodes(data, barcodesfile, region_num, newfile):
	print(f"Checking unique barcodes: BC{region_num+1}")
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
	
	f = open(str(newfile), 'w')
	f.write('Barcode,Count\n')
	for key, value in unique_barcodes.items():
		f.write(key + ',' + str(value))
		f.write('\n')
	f.close()
	

def check_combinations(data, barcodesfile, newfile):
	print('Checking unique sequences')
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
		if one_barcode:
			if seq_str in barcode_combos.keys():
				barcode_combos[seq_str] += count
			else:
				barcode_combos[seq_str] = count

	f = open(str(newfile), 'w')
	f.write('BC1,BC2,BC3,Counts')
	f.write('\n')
	for seq, count in barcode_combos.items():
		barcodes = seq.split('_')
		f.write('{},{},{},{}'.format(barcodes[0], barcodes[1], barcodes[2], str(count)))
		f.write('\n')
	f.close()



if __name__ == "__main__":
	barcodes_file = 'metadata/barcodes.fasta'

	filename = os.path.basename(snakemake.input[0])
	reads_file = snakemake.input[0]

	final_reads = frequency(reads_file)

	for i in range(len(snakemake.output)-1):
		check_unique_barcodes(final_reads, barcodes_file, i, snakemake.output[i])

	check_combinations(final_reads, barcodes_file, snakemake.output[-1])