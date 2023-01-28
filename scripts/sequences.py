import Bio.PDB
from Bio import AlignIO
import argparse
import random
import csv
import math
import requests


amino_acids = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     		'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
     		'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
     		'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

codons = {
   	'A': ['GCC', 'GCT', 'GCA', 'GCG'],
   	'C': ['TGC', 'TGT'],
	'D': ['GAC', 'GAT'],
   	'E': ['GAG', 'GAA'],
   	'F': ['TTC', 'TTT'],
   	'G': ['GGC', 'GGG', 'GGA', 'GGT'],
   	'H': ['CAC', 'CAT'],
   	'I': ['ATC', 'ATT', 'ATA'],
   	'K': ['AAG', 'AAA'],
   	'L': ['CTG', 'CTC', 'CTT', 'TTG', 'TTA', 'CTA'],
   	'M': ['ATG'],
   	'N': ['AAC', 'AAT'],
   	'P': ['CCC', 'CCT', 'CCA', 'CCG'],
   	'Q': ['CAG', 'CAA'],
   	'R': ['CGG', 'AGG', 'AGA', 'CGC', 'CGA', 'CGT'], 
   	'S': ['AGC', 'TCC', 'TCT', 'TCA', 'AGT', 'TCG'],
   	'T': ['ACC', 'ACA', 'ACT', 'ACG'],
   	'V': ['GTG', 'GTC', 'GTT', 'GTA'],
   	'W': ['TGG'],
   	'Y': ['TAC', 'TAT'],
   	'_': ['TGA', 'TAA', 'TAG'],
   	'-': ['---']
   	}


adapters = [('GTGAAACCGTCTCCACAT', 'GGAGACGCGTCATTAAAG'),
			('TAACTACCAGTACGTCTCC', 'GGAGACGACACCAATACAA'),
			('AAACACTTCAGCGTCTCCAT', 'GATAGGAGACGGATTCGGAA'),
			('GGGAATTTATCCGTCTCCA', 'AACCGGAGACGCACCCTGAAGT')]



def modify_alignment(pdbfile, main_chain, alignment):
	parser = Bio.PDB.PDBParser(QUIET=True)
	structure = parser.get_structure(pdbfile.split(".")[0] + "_chain", pdbfile)

	principle_chain = structure[0][main_chain]

	sequence = ''
	for residue in principle_chain:
		if residue.get_id()[0] == ' ' and residue.get_resname() in amino_acids.keys():
			sequence += amino_acids[residue.get_resname()]


	temp_alignment = []

	for seq in alignment:
		temp = ''

		for residue in seq:
			temp += str(residue)

		temp_alignment.append(temp)

	missed_residues_front = temp_alignment[0].find(sequence[:3])

	new_alignment = []

	for seq in temp_alignment:
		new_alignment.append(seq[missed_residues_front:])

	missed_residues_back = new_alignment[0].find(sequence[-3:])
	missed_addon = missed_residues_back + 3

	new_alignment_2 = []

	for idx in range(len(new_alignment)):
		new_alignment_2.append(new_alignment[idx][:missed_addon])
		
	return new_alignment_2, temp_alignment




def read_crossover_points():
	f = open("metadata/crossover_points.txt", "r")
	t = f.read()
	f.close()
	
	temp_points = t.split('\n')[0]
	
	points = [int(i) for i in temp_points.split(',')]
	return points




def get_protein(row, alignment):

	protein_data = [int(float(value)) for value in row]
	crossovers = read_crossover_points()

	middle_idx = int((len(protein_data)-2)/2)
	regions = [protein_data[idx] for idx in range(middle_idx, len(protein_data)-2)]
	alignment_len = len(alignment[0])

	key_tuples = [(0, crossovers[0])]
	protein = {}
	
	for crossover in range(1, len(crossovers)):
		tple = (crossovers[crossover - 1] + 1, crossovers[crossover])
		key_tuples.append(tple)

	key_tuples.append((crossovers[-1] + 1, alignment_len - 1))

	for idx in range(len(regions)):
		protein[key_tuples[idx]] = regions[idx]

	return protein



def get_parent(protein, residue):
	'''
	Returns parent for residue position of a recombined protein
	{(crossover_position1, crossover_position2): parent in that range}
	'''

	for key in protein.keys():
		if residue >= key[0] and residue <= key[1]:
			return protein[key]



def A(alignment, parent, k):
	return alignment[parent][k]



def get_sequence(protein, alignment):

	sequence = ''
	parents = list(range(len(alignment)))


	for i in range(len(alignment[0])):

		s_i = get_parent(protein, i)
		residue = A(alignment, s_i, i)

		if residue != "-":
			sequence += residue
		else:

			for parent in parents:
				if parent != s_i:
					residue = A(alignment, parent, i)

					if residue != '-':
						sequence += residue
						break

				if parent == len(parents) - 1:
					residue = A(alignment, s_i, i)
					sequence += residue


	return sequence




def complete_ends(temp_sequence, parent, temp_raw_alignment, modified_alignment):

	raw_alignment = []
	for seq in temp_raw_alignment:
		str_seq = ''

		for residue in seq:
			str_seq += str(residue)

		raw_alignment.append(str_seq)

	if temp_sequence in modified_alignment:
		sequence = ''
		for residue in range(len(temp_sequence)):
			if temp_sequence[residue] == '-':
				for seq in modified_alignment:
					if seq != temp_sequence:
						if seq[residue] != '-':
							sequence += seq[residue]
							break
			else:
				sequence += temp_sequence[residue]
	else:
		sequence = temp_sequence



	front_index = raw_alignment[0].find(modified_alignment[0])
	end_index = len(modified_alignment[0]) + front_index

	front_sequence = raw_alignment[parent][:front_index]
	end_sequence = raw_alignment[parent][end_index:]

	temp_sequence = ''.join([front_sequence, sequence, end_sequence])
	new_sequence = temp_sequence.replace('-', '')
	num_residues = len(new_sequence) - len(sequence)

	return new_sequence



def reverse_translate(sequence):

	bsmbi_sites = ['CGTCTC', 'GAGACG']
	dna_sequence_lst = []

	for residue in sequence:
		codon_options = codons[residue]

		for codon in codon_options:
			dna_sequence_lst.append(codon)
			dna_sequence = ''.join(dna_sequence_lst)

			if any(x in dna_sequence for x in bsmbi_sites):
				if codon == codon_options[-1]:
					break
				else:
					dna_sequence_lst.pop(-1)
			else:
				break

	return dna_sequence



def fill_gaps(alignment):

	residues = []
	new_alignment = [list(sequence) for sequence in alignment]

	for i in range(len(alignment[0])):
		for parent in range(len(alignment)):
			residue = A(alignment, parent, i)
			
			if residue == '-':
				residues.append(i)
				break

	for residue in residues:
		for parent in range(len(alignment)):
			new_residue = A(alignment, parent, residue)

			if new_residue != '-':
				for temp_parent in range(len(alignment)):
					temp_residue = A(alignment, temp_parent, residue)
					if temp_residue == '-':
						new_alignment[temp_parent][residue] = new_residue
				break

	final_alignment = [''.join(sequence_lst) for sequence_lst in new_alignment]
	return final_alignment
				





def find_conserved_regions(alignment, temp_raw_alignment, primer_length, fname):
	dna_sequences_temp = [complete_ends(seq, 0, temp_raw_alignment, alignment) for seq in alignment]
	dna_sequences = [reverse_translate(sequence) for sequence in dna_sequences_temp]
	conserved_regions = []

	for idx in range(len(dna_sequences[0]) - primer_length):
		sequences = []
		
		for sequence in dna_sequences:
			primer = sequence[idx:idx+primer_length]
			sequences.append(primer)

		if len(set(sequences)) == 1:
			region = (idx, idx+primer_length)
			conserved_regions.append(region)

	out = open(fname, 'w')
	for region in range(len(conserved_regions)):
		new_tuple = str(conserved_regions[region]).replace('(', '')
		new_tuple = new_tuple.replace(')', '')
		new_tuple = new_tuple.replace(' ', '')
		out.write(new_tuple)

		if region != len(conserved_regions)-1:
			out.write('\n')

	out.close()




def get_fragment_points(regions_file, alignment, temp_raw_alignment):

	f = open(regions_file, 'r')
	regions = (f.read()).split('\n')
	f.close()

	modified_alignment = [complete_ends(seq, 0, temp_raw_alignment, alignment) for seq in alignment]

	conserved_regions = [(int(region.split(',')[0]), int(region.split(',')[1])) for region in regions]
	unusable_lens = []
	fragment_points = []
	points_found = False
	fragment_len = 0
	num_fragments = 0
	buffer = 50 # need to add functionality so that the program test different buffers (up to 50) before it increases num_fragments
	max_frag_len = 300
	
	while not points_found:
		for divisor in range(2, (len(modified_alignment[0])*3)+1):
			fragment_len = int(math.ceil((len(modified_alignment[0])*3)/divisor))
			if fragment_len <= max_frag_len and fragment_len not in unusable_lens:
				num_fragments = divisor
				break

		min_distance = fragment_len

		for region in conserved_regions:
			lower = region[0]
			upper = region[1]

			# smarter way of doing this:
			# 1. go through each region in conserved_regions
			# 2. measure its ("upper") distance to the calculated length
			# 3. use index with min distance to fragment

			target = (fragment_len*len(fragment_points))+fragment_len
			distance = abs(upper - target)

			if distance < min_distance:
				min_distance = distance
				previous_upper = upper
			else:
				if min_distance > 60: # goes to next num_fragments if fragment lengths are too far apart
					break
				fragment_points.append(previous_upper)
				min_distance = fragment_len
				if len(fragment_points) + 1 == num_fragments:
					points_found = True
					break


#			if upper in range((fragment_len*len(fragment_points))+fragment_len-buffer, (fragment_len*len(fragment_points))+fragment_len+buffer):
#				fragment_points.append(upper)
#				if len(fragment_points) + 1 == num_fragments:
#					points_found = True
#					break
#			^^^^^^^***NOTE***: this is the old way of fragmenting (uses a buffer thing)

		if not points_found:
			unusable_lens.append(fragment_len)
			fragment_points = []

	first_len = fragment_points[0]
	last_len = (len(modified_alignment[0])*3)-fragment_points[-1]
	frag_lens = []

	if len(fragment_points) > 1:
		for idx in range(1, len(fragment_points)):
			frag_len = fragment_points[idx] - fragment_points[idx-1] + 1 # add one to account for index start at 0
			frag_lens.append(frag_len)

	frag_lens.insert(0, first_len)
	frag_lens.append(last_len)

	print('Fragment Lengths: {}'.format(', '.join([str(a) for a in frag_lens])))
	return fragment_points


def add_filler(idx):

	bases = ['A', 'T', 'C', 'G']
	bsmbi_sites = ['CGTCTC', 'GAGACG']
	front_part = ['GGATTAGACGACAGGTGG', 'TGGGAGTCTATCACCCCT', 'TGCACAAGCAATTGACAA']
	end_part = ['ATTGCTAGCCCTTGAACG', 'AGGGCCCATATCTGGAAA', 'CTTCTACCCATCTTCCGA']



	sequence = ''

	while len(sequence) < 0:
		base = random.choice(bases)
		temp_sequence = front_part[idx]+sequence+base+end_part[idx]

		if not any(x in temp_sequence for x in bsmbi_sites) and 'ATG' not in temp_sequence and 'TAC' not in temp_sequence:
			sequence += base

	sequence = front_part[idx]+sequence+end_part[idx]
	print(len(sequence))
	return sequence





def generate_oligos(csvfile, alignment, raw_alignment):

	find_conserved_regions(alignment, raw_alignment, 16, 'conserved_regions.txt')
	fragment_points = get_fragment_points('conserved_regions.txt', alignment, raw_alignment)
	fragment_points.insert(0, 0)
	bsmbi_sites = ['CGTCTC', 'GAGACG']
	barcodes = []
	used_barcodes = []
	final_barcodes = []
	sequences = []
	oligos = []
	raw_oligos = []
	random.seed(20)

	csv_file = open(csvfile, 'r')
	data = list(csv.reader(csv_file))

	print('Reverse translating chimera amino acid sequences')
	for index in tqdm(range(1, len(data))):
		protein = get_protein(data[index], alignment)
		temp_sequence = get_sequence(protein, alignment)
		sequence = complete_ends(temp_sequence, 0, raw_alignment, alignment)
		sequence = reverse_translate(sequence)
		sequences.append(sequence)

		if index == 1:
			fragment_points.append(len(sequence))

	csv_file.close()
	sequences = list(set(sequences))
	sequences_copies = []

	for sequence in sequences: # used to make three copies of every oligo
		sequences_copies.append(sequence)
		sequences_copies.append(sequence)
		sequences_copies.append(sequence)

	sequences = sequences_copies

	# add on barcodes
	data = requests.get('https://raw.githubusercontent.com/churchlab/AAV_fitness_landscape/master/data/meta/AAV2scan_chip_lookup_table.txt')
	data = data.text
	data = data.split('\n')[1:800]

	print('Checking barcodes')
	for line in tqdm(data):
		barcode = line.split(',')[2]
		if not any(barcode in x for x in sequences) and len(str(barcode)) == 20:
			barcodes.append(str(barcode))

	for idx in range(1, len(fragment_points)):
		oligo_region = []
		raw_oligo_region = []
		sequence_region = []
		barcode_region = []
		filler = add_filler(idx-1)
		print('Generating oligos: region {}'.format(idx))
		for sequence in tqdm(sequences):
			adapter_pair = adapters[idx-1]

			if idx != len(fragment_points)-1:
				oligo = sequence[fragment_points[idx-1]:fragment_points[idx]] # ***IMPORTANT***: do sequence[fragment_points[idx-1]:fragment_points[idx]+4] if sticky ends between oligos are needed
			else:
				oligo = sequence[fragment_points[idx-1]:fragment_points[idx]]

			found_barcode = ''
			
			for barcode in barcodes:
				test_oligo = adapter_pair[0] + oligo + filler + barcode + adapter_pair[1]
				if barcode not in (oligo+filler) and barcode not in used_barcodes:
					occurences = 0
					for site in bsmbi_sites:
						occurences += test_oligo.count(site)
					if occurences <= 2:
						found_barcode = barcode
						break

			if sequence_region.count(oligo) < 3:
				sequence_region.append(oligo)

				raw_oligo_region.append(oligo)

				oligo = adapter_pair[0] + oligo + filler + barcode + adapter_pair[1]
				oligo_region.append(oligo)

				used_barcodes.append(found_barcode)
				barcode_region.append(found_barcode)

		raw_oligos.append(raw_oligo_region)
		oligos.append(oligo_region)
		final_barcodes.append(barcode_region)


	print('Writing oligo sequences to fasta file')
	f = open('oligos2.fasta', 'w')
	d = open('barcodes2.fasta', 'w')
#	a = open('oligos2.csv', 'w')
	b = open('oligos_raw.fasta', 'w')
	label = 1
	copy = 1
	for idx in range(len(oligos)):
		start = fragment_points[idx] + 1
		end = fragment_points[idx+1]

		for oligo_num in range(len(oligos[idx])):
			f.write('>oligo{}_region{}_{}-{}_copy{}'.format(label, (idx+1), start, end, copy))
			f.write('\n')
			f.write(oligos[idx][oligo_num])
			f.write('\n')

			b.write('>oligo{}_region{}_{}-{}_copy{}'.format(label, (idx+1), start, end, copy))
			b.write('\n')
			b.write(raw_oligos[idx][oligo_num])
			b.write('\n')

#			a.write(oligos[idx][oligo_num])
#			a.write('\n')

			d.write('>barcode{}_region{}_{}-{}_copy{}'.format(label, (idx+1), start, end, copy))
			d.write('\n')
			d.write(final_barcodes[idx][oligo_num])
			d.write('\n')

			copy += 1

			if copy > 3:
				label += 1
				copy = 1


	f.close()
	d.close()
	b.close()
#	a.close()
	print('Complete')



if __name__ == "__main__":

	parser = argparse.ArgumentParser()
	parser.add_argument('-a', '--align', required=True, help='Alignment file name')
	parser.add_argument('-p', '--pdb', action='store', required=True, help='PDB file name')
	parser.add_argument('--chain', action='store', default='A', help='Chain to make chimeras')
	parser.add_argument('-c', '--csv', action='store', required=True, help='CSV File Library')
	args = parser.parse_args()

	pdb_file = "/Users/architc/Documents/Rice/Research/Silberg Lab/Scripts/PDB_Files/{}".format(str(args.pdb))
#	pdb_file = "../../PDB_Files/{}".format(str(args.pdb))
	alignment_file = "/Users/architc/Documents/Rice/Research/Silberg Lab/Scripts/Alignment_Files/{}".format(str(args.align))
#	alignment_file = "../../Alignment_Files/{}".format(str(args.align))
	csv_file = "../metadata/{}".format(str(args.csv))
#	csv_file = "../../Calibration_Libraries/OBD Library/{}".format(str(args.csv))

	raw_alignment = AlignIO.read(alignment_file, "clustal")
	modified_alignment, raw_alignment2 = modify_alignment(pdb_file, args.chain, raw_alignment)
#	alignment = fill_gaps(temp_alignment)
	
	generate_oligos(csv_file, modified_alignment, raw_alignment)
#	csv_file = open(csv_file, 'r')
#	data = list(csv.reader(csv_file))
#	sequences = []

#	for index in tqdm(range(1, len(data))):
#		protein = get_protein(data[index], modified_alignment)
#		temp_sequence = get_sequence(protein, modified_alignment)
#		sequence = complete_ends(temp_sequence, 0, raw_alignment, modified_alignment)
	#	sequence = reverse_translate(sequence)
#		sequences.append(sequence)

#	sequences = list(set(sequences))
#	a = open('sequences_test.csv', 'w')
#	csv_file.close()
#	for sequence in sequences:
#		a.write(sequence)
#		a.write('\n')
#	a.close()






