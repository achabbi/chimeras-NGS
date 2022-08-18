import re
import csv
from Bio import SeqIO


if __name__ == "__main__":

    oligos_file = 'metadata/oligos.fasta'
    barcodes_file = 'metadata/barcodes.fasta'
    barcodes_csv = snakemake.input[0]
    file_path = snakemake.output[0]

    csv_file = open(barcodes_csv, 'r')
    data = list(csv.reader(csv_file))
    data.pop(0)

    barcode_ref = {'':0}
    for rec in SeqIO.parse(barcodes_file, 'fasta'):
        id = int(re.compile(r"barcode(\d+)_").search(rec.name).group(1))
        seq = str(rec.seq)
        barcode_ref[seq] = id

    mutants_dict = {}
    for entry in data:
        id1 = str(barcode_ref[entry[0]])
        id2 = str(barcode_ref[entry[1]])
        id3 = str(barcode_ref[entry[2]])
        count = int(entry[3])

        mutant = "_".join([id1, id2, id3])
        if mutant in mutants_dict.keys():
            mutants_dict[mutant] += count
        else:
            mutants_dict[mutant] = count
    
    with open(file_path, 'w') as csv_file:  
        writer = csv.writer(csv_file)
        writer.writerow(['OBD1 ID', 'OBD2 ID', 'OBD3 ID', 'Count'])
        for key, value in mutants_dict.items():
            ids = key.split('_')
            ids = ['' if id == '0' else id for id in ids]
            writer.writerow([ids[0], ids[1], ids[2], value])