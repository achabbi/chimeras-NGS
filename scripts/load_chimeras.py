import re
import csv
import itertools
import sequences as sq
from Bio import SeqIO
from Bio import AlignIO

adapters = [('GTGAAACCGTCTCCACAT', 'GGAGACGCGTCATTAAAG'),
			('TAACTACCAGTACGTCTCC', 'GGAGACGACACCAATACAA'),
			('AAACACTTCAGCGTCTCCAT', 'GATAGGAGACGGATTCGGAA')]

filler_front = ['GGATTAGACGACAGGTGG', 'TGGGAGTCTATCACCCCT', 'TGCACAAGCAATTGACAA']
filler_end = ['ATTGCTAGCCCTTGAACG', 'AGGGCCCATATCTGGAAA', 'CTTCTACCCATCTTCCGA']


def get_oligo_sequence(ids):
    oligos = {}
    for rec in SeqIO.parse(oligos_file, 'fasta'):
        id = int(re.compile(r"oligo(\d+)_").search(rec.name).group(1))
        if id not in oligos.keys():
            seq = str(rec.seq)
            for idx in range(len(adapters)):
                if adapters[idx][0] in seq and filler_front[idx] in seq:
                    start = seq.index(adapters[idx][0]) + len(adapters[idx][0])
                    end = seq.index(filler_front[idx])
                    seq = seq[start:end]
                    break
            oligos[id] = seq
    
    sequence = ''
    for id in ids:
        sequence += oligos[id]
    return sequence


def get_barcode_combinations(ids):
    barcodes = {}
    for rec in SeqIO.parse(barcodes_file, 'fasta'):
        id = int(re.compile(r"barcode(\d+)_").search(rec.name).group(1))
        seq = str(rec.seq)
        if id not in barcodes.keys():
            barcodes[id] = [seq]
        else:
            barcodes[id].append(seq)

    combination = []
    for id in ids:
        combination.append(barcodes[id])
    return combination
    


def get_barcode_distance(barcode_combinations):
    wt = 'ATATAGTGAACCCCGCCCCATTGGCACCAGATACCTGACTCGTAATCTGTAATAATTGCTTGTTAATCAATAAACCGTTTAATTCGTTTCAGTTGAACTTTGGTCTCTGCGAAGGGCGAATTCGTTTAAACCTGC'
    combinations = list(itertools.product(*barcode_combinations))
    hamming_distances = []

    for combination in combinations:
        complete_barcode = ''.join(combination)
        distances = []
        window = len(complete_barcode)

        for i in range(len(wt)-window):
            subseq = wt[i:i+window]
            distance = sum(complete_barcode[i] != subseq[i] for i in range(len(subseq)))
            distances.append(distance)
        
        min_distance = min(distances)
        hamming_distances.append(min_distance)

    return min(hamming_distances)





if __name__ == "__main__":
    pdb_file = 'metadata/4ZO0.pdb'
    alignment_file = 'metadata/aav2_aav4_aav5_obd.clustal_num'
    oligos_file = 'metadata/oligos.fasta'
    barcodes_file = 'metadata/barcodes.fasta'
    chimeras_file = 'metadata/chimeras.csv'
    mutants_csv = snakemake.input[0]
    file_path = snakemake.output[0]

    raw_alignment = AlignIO.read(alignment_file, "clustal")
    alignment, _ = sq.modify_alignment(pdb_file, 'A', raw_alignment)

    mutants_file = open(mutants_csv, 'r')
    mutants_data = list(csv.reader(mutants_file))
    chimeras_data = csv.reader(open(chimeras_file, 'r'))

    with open(file_path, 'w') as csv_file:  
        writer = csv.writer(csv_file)
        writer.writerow(['OBD1 ID', 'OBD2 ID', 'OBD3 ID', 'region_1', 'region_2', 'region_3', 'region_4', 'region_5', 'region_6', 'region_7', 'region_7', 'region_8', 'region_8',
                            'region_9', 'region_10', 'region_11', 'contacts_broken', 'chimera_distance', 'barcode_distance', 'frequency'])

        for row1 in mutants_data[1:-2]:
            ids = [int(id) for id in row1[:3]]
            frequency = row1[4]
            chimera_sequence = get_oligo_sequence(ids)
            barcode_combinations = get_barcode_combinations(ids)

            for row2 in chimeras_data:
                if "crossover_1" in row2:
                    continue
                protein = sq.get_protein(row2, alignment)
                temp_sequence = sq.get_sequence(protein, alignment)
                sequence = sq.complete_ends(temp_sequence, 0, raw_alignment, alignment)
                sequence = sq.reverse_translate(sequence)

                if chimera_sequence == sequence:
                    barcode_distance = get_barcode_distance(barcode_combinations)
                    row = []
                    row.extend(ids)
                    row.extend(row2[10:])
                    row.append(barcode_distance)
                    row.append(frequency)
                    writer.writerow(row)
                    break