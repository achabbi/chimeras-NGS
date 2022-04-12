import itertools

sequences = []
quality = []

f = open('g180_S1_L001_R1_001.fastq', 'r')
lines = f.readlines()
f.close()
data = [line.splitlines()[0] for line in lines]

for idx in range(len(data)):
	if data[idx] == '+':
		sequences.append(data[idx-1])
		quality.append(data[idx+1])


f = open('/Users/architc/Documents/Rice/Research/Silberg Lab/Scripts/Final_Scripts/OBD Oligos/barcodes2.fasta')
lines = f.readlines()
f.close()
data = [line.splitlines()[0] for line in lines]

identifiers = []
barcodes = []

for idx in range(len(data)):
	if '>' in data[idx]:
		identifiers.append(data[idx])
		barcodes.append(data[idx+1])

regions = []

current_region = identifiers[0].partition('region')[2][0]
region_list = []

for idx in range(len(identifiers)):
	if current_region != identifiers[idx].partition('region')[2][0] or idx == len(identifiers)-1:
		regions.append(region_list)
		current_region = identifiers[idx].partition('region')[2][0]
		region_list = []

	region_list.append(barcodes[idx])

oligo_barcodes_tuples = list(itertools.product(regions[0], regions[1], regions[2]))
oligo_barcodes = [tple[0]+tple[1]+tple[2] for tple in oligo_barcodes_tuples]
# note these are all oligo barcodes (bc1+bc2+bc3), need ones that are covered by library only
# need to read original csv file, get the sequence, see which oligo combinations are used