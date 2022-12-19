import re
import csv
from Bio import SeqIO

enrichment_file = ''
oligos_file = ''
mutants = {}

bsmbi_sites = ['CGTCTC', 'GAGACG']

adapters = [('GTGAAACCGTCTCCACAT', 'GGAGACGCGTCATTAAAG'),
			('TAACTACCAGTACGTCTCC', 'GGAGACGACACCAATACAA'),
			('AAACACTTCAGCGTCTCCAT', 'GATAGGAGACGGATTCGGAA'),
			('GGGAATTTATCCGTCTCCA', 'AACCGGAGACGCACCCTGAAGT')]


def remove_adapters(seq, id):
    # adapter_pair[0] + oligo + filler + barcode + adapter_pair[1]
    
    pass


with open(enrichment_file) as csvfile:
    reader = csv.reader(csvfile)
    for row in reader:
        try:
            mutants[(int(row[0]), int(row[1]), int(row[2]))] = float(row[4]) # assigns id combination to its score
        except:
            pass
            
        if row[0] == 'WT':
            wt_naive_counts = float(row[4])
        elif row[0] == 'K340H':
            k340h_naive_counts = float(row[4])

oligos = {}
for rec in SeqIO.parse(oligos_file, 'fasta'):
    id = int(re.compile(r"oligo(\d+)_").search(rec.name).group(1))
    if id not in oligos.keys():
        seq = str(rec.seq)


        oligos[id] = seq


# 1. load all sequences from oligos.fasta into dictionary -> {id: sequence}
#       - need to remove extra DNA sequences from oligo (barcode, adapters, cut sites, etc.)
# 2. convert mutant id combo into a DNA sequence and store all sequences
#       - use dictionary from 1. to get oligo sequence
#       - join oligo sequences together to make chimera sequence
#       - match chimera sequence with its enrichment score in new dictionary
# 3. go through chimeras.csv and convert each row to a chimera sequence with functions from sequences.py
#       - store sequences in list (may take up too much memory)
#       - other option: convert each row to chimera sequence, then immediately check new dictionary from 2. to check if present
# 4. create new csv with chimera sequence (or crossovers/inherited regions from chimeras.csv), SCHEMA score, distance m, and enrichment score