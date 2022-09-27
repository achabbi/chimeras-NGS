import re
import csv
from Bio import SeqIO


if __name__ == "__main__":
    """
    Gets counts of each individual oligo instead of barcode.
    Necessary as there are 3 barcodes per oligo.

    Input: csv file with individual barcode sequences and counts through Snakemake
    Output: csv file with individual oligo ids and counts through Snakemake
    """

    oligos_file = 'metadata/oligos.fasta'
    barcodes_file = 'metadata/barcodes.fasta'
    barcodes_csv = snakemake.input[0] # Snakemake input barcode counts csv file
    file_path = snakemake.output[0] # Snakemake output oligo counts csv file

    csv_file = open(barcodes_csv, 'r')
    data = list(csv.reader(csv_file))
    data.pop(0)

    barcode_dict = {}
    for entry in data:
        barcode_dict[entry[0]] = int(entry[1])

    oligo_counts = {}
    total_counts = 0
    for rec in SeqIO.parse(barcodes_file, 'fasta'):
        id = int(re.compile(r"barcode(\d+)_").search(rec.name).group(1))
        seq = str(rec.seq)
        total_counts += 1

        if seq in barcode_dict.keys():
            count = barcode_dict[seq]

            if id in oligo_counts.keys():
                oligo_counts[id] += count
            else:
                oligo_counts[id] = count
    
    with open(file_path, 'w') as csv_file:  
        writer = csv.writer(csv_file)
        writer.writerow(['ID', 'Count', 'Frequency'])
        for key, value in oligo_counts.items():
            writer.writerow([key, value, (value/total_counts)])
    
