from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
import itertools
import matplotlib.pyplot as plt
from collections import Counter
from tqdm import tqdm
import sequences as sq


def filter_quality(filename):
	good_quality = []

	for rec in SeqIO.parse(filename, 'fastq'):
		if min(rec.letter_annotations['phred_quality']) >= 20:
			good_quality.append(rec)

	return good_quality


def find_matches(forward, reverse):
# returns list of Sequences objects (not SeqRecord)
	raw_seq_forward = [sequence.seq for sequence in forward]
	raw_seq_reverse = [sequence.seq.reverse_complement() for sequence in reverse]

	found_forward = []
	found_reverse = []

	filtered_forward = []
	filtered_reverse = []

	for seq in raw_seq_forward:
		if seq in raw_seq_reverse and seq not in found_forward:
			found_forward.append(seq)

	for seq in raw_seq_reverse:
		if seq in raw_seq_forward and seq not in found_reverse:
			found_reverse.append(seq)

	for seqf in forward:
		if seqf.seq in found_forward:
			filtered_forward.append(seqf.seq)

	for seqr in reverse:
		if seqr.seq.reverse_complement() in found_reverse:
			filtered_reverse.append(seqr.seq)

	return filtered_forward, filtered_reverse


def combine_forward_reverse(forward, reverse):
	# input is dictionary of histogram type ({sequence1: frequency, sequence2: frequency, etc})
	final_dict = forward.copy()

	for seq_reverse, count in reverse.items():
		final_dict[seq_reverse.reverse_complement()] += count

	return final_dict


def combine_dictionaries(stock1, stock2):
# for combining final sequences from both glycerol stocks
	final_dict = stock1.copy()

	for seq2, count2 in stock2.items():
		if seq2 in stock1.keys():
			final_dict[seq2] += count2
		else:
			final_dict[seq2] = count2

	return final_dict




def histogram(seq_dict):
	# makes a histogram from sequence dictionary
	sequences = []
	counts = []

	for seq, count in seq_dict.items():
		sequences.append(str(seq))
		counts.append(count)

	plt.bar(sequences, counts)
	plt.show()



def recombined_barcodes(filename, num_regions):
# returns list of SeqRecord objects
	regions = [[] for _ in range(num_regions)]

	for rec in SeqIO.parse(filename, 'fasta'):
		region = int(rec.name.partition('region')[2][0])
		regions[region-1].append(rec)

	barcodes_combined = itertools.product(regions[2], regions[1], regions[0])
	final_barcodes = []

	for tple in barcodes_combined:
		sequence = ''
		oligo_id = ''

		for barcode in tple:
			sequence += str(barcode.seq)
			oligo_id += str(barcode.name).partition('barcode')[2][0]

			if barcode.seq != tple[-1].seq:
				oligo_id += ','

		seq = Seq(sequence)
		seq_rec = SeqRecord(seq, id=oligo_id)
		final_barcodes.append(seq_rec)

	return final_barcodes





def get_barcodes(csvfile, oligosfile, barcodesfile, alignment, raw_alignment):
	
	csv_file = open(csvfile, 'r')
	data = list(csv.reader(csv_file))
	sequences = []

	for index in tqdm(range(1, len(data))):
		protein = sq.get_protein(data[index], alignment)
		temp_sequence = sq.get_sequence(protein, alignment)
		sequence = sq.complete_ends(temp_sequence, 0, raw_alignment, alignment)
		sequence = sq.reverse_translate(sequence)
		sequences.append(sequence)

	csv_file.close()
	sequences = list(set(sequences))
	fragment_points = sq.get_fragment_points('conserved_regions.txt', alignment, raw_alignment)

	# get oligos and sort by region
	regions = [[] for _ in range(3)]
	for rec in SeqIO.parse(oligosfile, 'fasta'):
		region = int(rec.name.partition('region')[2][0])
		regions[region-1].append(rec)


	barcodes = [[] for _ in range(3)]
	for rec in SeqIO.parse(barcodesfile, 'fasta'):
		region = int(rec.name.partition('region')[2][0])
		barcodes[region-1].append(rec)


	temp_barcodes = []
	for chimera in tqdm(sequences):
		frag1 = chimera[:fragment_points[0]]
		frag2 = chimera[fragment_points[0]:fragment_points[1]]
		frag3 = chimera[fragment_points[1]:]

		combined_barcode = []

		for oligo_seq in regions[0]:
			if frag1 in str(oligo_seq.seq):
				oligo_name = str(oligo_seq.name)
				barcode_name = 'barcode' + oligo_name.partition('oligo')[2]

				for barcode in barcodes[0]:
					if str(barcode.name) == barcode_name:
						combined_barcode.append(barcode)
						break
				break

		for oligo_seq in regions[1]:
			if frag2 in str(oligo_seq.seq):
				oligo_name = str(oligo_seq.name)
				barcode_name = 'barcode' + oligo_name.partition('oligo')[2]

				for barcode in barcodes[1]:
					if str(barcode.name) == barcode_name:
						combined_barcode.append(barcode)
						break
				break

		for oligo_seq in regions[2]:
			if frag3 in str(oligo_seq.seq):
				oligo_name = str(oligo_seq.name)
				barcode_name = 'barcode' + oligo_name.partition('oligo')[2]

				for barcode in barcodes[2]:
					if str(barcode.name) == barcode_name:
						combined_barcode.append(barcode)
						break
				break

		temp_barcodes.append(combined_barcode)

	final_barcodes = []
	for bc_combos in temp_barcodes:
		bc1 = str(bc_combos[0].seq)
		bc2 = str(bc_combos[1].seq)
		bc3 = str(bc_combos[2].seq)
		final_bc = bc3+bc2+bc1 # this is wrong, should include filler in between
		final_barcodes.append(final_bc)

	return final_barcodes



def check_coverage(data, barcodes):
	# checks which complete barcodes were found in the NGS reads
	sequences = [str(seq) for seq in list(data.keys())]
	found_sequences = []
	found_barcodes = []

	for barcode in barcodes:
		for sequence in sequences:
			if barcode in sequence:
				found_barcodes.append(barcode)
				found_sequences.append(sequence)

	print(len(found_barcodes)) # how many barcodes were found
	print(len(barcodes))




# ----Main code begins-------


quality_forward_g180 = filter_quality('g180_S1_L001_R1_001.fastq')
quality_reverse_g180 = filter_quality('g180_S1_L001_R2_001.fastq')

quality_forward_g181 = filter_quality('g181_S2_L001_R1_001.fastq')
quality_reverse_g181 = filter_quality('g181_S2_L001_R2_001.fastq')

filtered_forward_g180, filtered_reverse_g180 = find_matches(quality_forward_g180, quality_reverse_g180)
filtered_forward_g181, filtered_reverse_g181 = find_matches(quality_forward_g181, quality_reverse_g181)

occurrences_forward_g180 = dict(Counter(filtered_forward_g180))
occurrences_reverse_g180 = dict(Counter(filtered_reverse_g180))

occurrences_forward_g181 = dict(Counter(filtered_forward_g181))
occurrences_reverse_g181 = dict(Counter(filtered_reverse_g181))

reads_g180 = combine_forward_reverse(occurrences_forward_g180, occurrences_reverse_g180)
reads_g181 = combine_forward_reverse(occurrences_forward_g181, occurrences_reverse_g181)

combined_reads = combine_dictionaries(reads_g180, reads_g181)
# combined_reads is dictionary containing unique sequences from both g180 and g181 after filtering
# format of dictionary is following: {sequence1: count, 'sequence2: count, etc.}, 
# where sequence1, sequence2 are SeqRecord objects (can be converted to string with str(sequence1))




oligos_file = 'oligos2.fasta'
barcodes_file = 'barcodes2.fasta'
chimera_file = 'chimeras.csv'
pdb_file = '4ZO0.pdb'
alignment_file = 'aav2_aav4_aav5_obd.clustal_num'

raw_alignment = AlignIO.read(alignment_file, "clustal")
modified_alignment, raw_alignment2 = sq.modify_alignment(pdb_file, 'A', raw_alignment)

barcodes = get_barcodes(chimera_file, oligos_file, barcodes_file, modified_alignment, raw_alignment) # gets all complete barcodes from chimera library
found_barcodes = check_coverage(combined_reads, barcodes) # gets barcodes from above were found in NGS data



# code below just writes unique sequences to csv file

#sequences = []
#for seq in combined_reads.keys():
#	sequence = str(seq)
#	sequences.append(sequence)


#a = open('unique_sequences.csv', 'w')
#for sequence in sequences:
#	a.write(sequence)
#	a.write('\n')
#a.close()




