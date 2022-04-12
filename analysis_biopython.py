from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
import itertools
import matplotlib.pyplot as plt
from collections import Counter
from tqdm import tqdm
import csv
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



def filter_barcodes(sequence_dict, barcodes):
# get the barcodes from all possilibities that match from the data
# sequence_dict is histogram like dictionary
	sequences = list(sequence_dict.keys())

	def longestSubstringFinder(string1, string2):
	    '''Find the longest matching word'''
	    answer = ""
	    len1, len2 = len(string1), len(string2)
	    for i in range(len1):
	        for j in range(len2):
	            lcs_temp=0
	            match=''
	            while ((i+lcs_temp < len1) and (j+lcs_temp<len2) and string1[i+lcs_temp] == string2[j+lcs_temp]):
	                match += string2[j+lcs_temp]
	                lcs_temp+=1         
	            if (len(match) > len(answer)):
	                answer = match              
	    return answer

	def listCheck(main):
	    '''control the input for finding substring in a list of words'''
	    string1 = main[0]
	    result = []
	    for i in range(1, len(main)):
	        string2 = main[i]
	        res1 = longestSubstringFinder(string1, string2)
	        res2 = longestSubstringFinder(string2, string1)
	        result.append(res1)
	        result.append(res2)
	    result.sort()
	    return result

	found_barcodes = []
	for barcode in tqdm(barcodes):
		barcode_seq = str(barcode.seq) # SeqRecord object so need .seq
		for seq in sequences:
			raw_seq = str(seq) # no need for .seq since already Seq object

			main = [barcode_seq, raw_seq]
			first_answer = listCheck(main)
			final_answer  = []

			for item1 in first_answer:
			    string1 = item1
			    double_check = True
			    for item2 in main:
			        string2 = item2
			        if longestSubstringFinder(string1, string2) != string1:
			            double_check = False
			    if double_check:
			        final_answer.append(string1)

			overlap = len(list(set(final_answer))[0])
			if overlap >= 50:
				found_barcodes.append(barcode)

	return found_barcodes



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
		final_bc = bc3+bc2+bc1
		final_barcodes.append(final_bc)

	return final_barcodes



def check_coverage(data, barcodes):
	sequences = [str(seq) for seq in list(data.keys())]
	found_sequences = []
	found_barcodes = []

	for barcode in barcodes:
		for sequence in sequences:
			if barcode in sequence:
				found_barcodes.append(barcode)
				found_sequences.append(sequence)

	print(len(found_barcodes))
	print(len(barcodes))




'''
possibly more efficient option cause current option for finding barcodes takes too long to run:

go through merged.csv chimeras, reverse translate, then find which 
oligos make up each chimera (going to have to remove extra stuff like fillers, bsmbi)
then, find corresponding barcode for each oligo (might run into a problem here bc of three copies)
since you know which barcode makes up each oligo then, create BC3+2+1 for each chimera sequence
then, run filter_barcodes with BC3+2+1 sequences and NGS sequences to see which of the
BC3+2+1 sequences are covered in the NGS data

hopefully this idea should be faster bc won't have to check over a million combinations
DON'T DELETE filter_barcodes
'''

# ASK: how much of the NGS barcode region has to overlap with th
# does the entire barcode region need to be there?
# if entire region does not need to be there,
# then can loop through each barcode region and check if barcode in NGS_barcode_sequence
# if barcodes from all three regions there, then record which ones and find corresponding oligos


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

combined_dict_g180 = combine_forward_reverse(occurrences_forward_g180, occurrences_reverse_g180)
combined_dict_g181 = combine_forward_reverse(occurrences_forward_g181, occurrences_reverse_g181)

combined_dict = combine_dictionaries(combined_dict_g180, combined_dict_g181)
#histogram(combined_dict_g180)

'''
oligos_file = '/Users/architc/Documents/Rice/Research/Silberg Lab/Scripts/Final_Scripts/OBD Oligos/oligos2.fasta'
barcodes_file = '/Users/architc/Documents/Rice/Research/Silberg Lab/Scripts/Final_Scripts/OBD Oligos/barcodes2.fasta'
csv_file = '/Users/architc/Documents/Rice/Research/Silberg Lab/Scripts/Calibration_Libraries/OBD Library/merged.csv'
pdb_file = '/Users/architc/Documents/Rice/Research/Silberg Lab/Scripts/PDB_Files/4ZO0.pdb'
alignment_file = '/Users/architc/Documents/Rice/Research/Silberg Lab/Scripts/Alignment_Files/aav2_aav4_aav5_obd.clustal_num'

raw_alignment = AlignIO.read(alignment_file, "clustal")
modified_alignment, raw_alignment2 = sq.modify_alignment(pdb_file, 'A', raw_alignment)

barcodes = get_barcodes(csv_file, oligos_file, barcodes_file, modified_alignment, raw_alignment)
#check_coverage(combined_dict, barcodes)
'''

sequences = []
for seq in combined_dict.keys():
	sequence = str(seq)
	sequences.append(sequence)


a = open('unique_sequences.csv', 'w')
for sequence in sequences:
	a.write(sequence)
	a.write('\n')
a.close()




