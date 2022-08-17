import csv
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
import gzip
import matplotlib.pyplot as plt
from collections import Counter
from tqdm import tqdm
import sequences as sq


def filter_quality(file):
	good_quality = []
	with gzip.open(file, 'rt') as filename:
		for rec in tqdm(SeqIO.parse(filename, 'fastq')):
#			if min(rec.letter_annotations['phred_quality']) >= 20:
			good_quality.append(rec)
	return good_quality


def find_matches(forward, reverse):
# returns list of Sequences objects (not SeqRecord)
	raw_seq_forward = [sequence.seq for sequence in forward]
	raw_seq_reverse = [sequence.seq.reverse_complement() for sequence in reverse]

	filtered_forward = []
	filtered_reverse = []
	intersect = list(set(raw_seq_forward) & set(raw_seq_reverse))

	for seq in tqdm(intersect):
		for _ in range(raw_seq_forward.count(seq)):
				filtered_forward.append(seq)
		for _ in range(raw_seq_reverse.count(seq)):
				filtered_reverse.append(seq.reverse_complement())

	return filtered_forward, filtered_reverse


def combine_forward_reverse(forward, reverse):
	# input is dictionary of histogram type ({sequence1: frequency, sequence2: frequency, etc})
	final_dict = forward.copy()
	for seq_reverse, count in reverse.items():
		final_dict[seq_reverse.reverse_complement()] += count
	return final_dict




def check_unique_barcodes(data, barcodesfile, region_num, newfile):
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
	
	f = open(str(newfile+'.csv'), 'w')
	for key, value in unique_barcodes.items():
		f.write(key + ',' + str(value))
		f.write('\n')
	f.close()
	

def check_combinations(data, barcodesfile, newfile):
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

	f = open(str(newfile+'.csv'), 'w')
	f.write('BC1,BC2,BC3,Counts')
	f.write('\n')
	for seq, count in barcode_combos.items():
		barcodes = seq.split('_')
		f.write('{},{},{},{}'.format(barcodes[0], barcodes[1], barcodes[2], str(count)))
		f.write('\n')
	f.close()



# ----Main code begins-------
barcodes_file = 'barcodes.fasta'
path = '/Users/architc/Documents/Rice/Research/Silberg Lab/08152022/'

forward_reads = path + 'test1.fastq.gz'
reverse_reads = path + 'test2.fastq.gz'

print('Filtering reads by quality')
quality_forward = filter_quality(forward_reads)
quality_reverse = filter_quality(reverse_reads)

print(len(quality_forward))
print(len(quality_reverse))

print('Finding forward and reverse read matches')
filtered_forward, filtered_reverse = find_matches(quality_forward, quality_reverse)

occurrences_forward = dict(Counter(filtered_forward))
occurrences_reverse = dict(Counter(filtered_reverse))

final_reads = combine_forward_reverse(occurrences_forward, occurrences_reverse)

#for i in range(3):
#	check_unique_barcodes(final_reads, barcodes_file, i, str(i))

check_combinations(final_reads, barcodes_file, 'test')

