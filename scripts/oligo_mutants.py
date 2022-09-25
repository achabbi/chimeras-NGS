import re
import csv
import gzip
from Bio import SeqIO


def barcodes_present(seq, barcodes):
    for barcode in barcodes:
        if barcode in seq:
            return True
    return False


if __name__ == "__main__":
    """
    Gets counts of each unique mutant formed by oligos instead of barcodes.
    Necessary as there are 3 barcodes per oligo.

    Input: csv file with chimeric barcode sequences and counts through Snakemake
    Output: csv file with chimeric oligo sequences and counts through Snakemake,
            with each row being a chimeric oligo sequence formed from ids corresponding 
            to its count. Last two rows are wildtype and K340H counts.
    """

    oligos_file = 'metadata/oligos.fasta'
    barcodes_file = 'metadata/barcodes.fasta'
    barcodes_csv = snakemake.input[0] # Snakemake input barcodes csv file
    reads_file = snakemake.input[1] # Snakemake input fastq.gz reads file
    file_path = snakemake.output[0] # Snakemake output mutants csv file

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


    # ------CALCULATES WILDTYPE and K340H COUNTS------
    wt_count = 0
    k340h_count = 0
    k340h_barcode = 'TACGTAGCTAGACGACTGCT'

    barcodes = []
    for rec in SeqIO.parse(barcodes_file, 'fasta'):
        barcodes.append(str(rec.seq))

    with gzip.open(reads_file, 'rt') as filename:
        for rec in SeqIO.parse(filename, 'fastq'):
            seq = str(rec.seq)
            if k340h_barcode in seq:
                k340h_count += 1
            elif not barcodes_present(seq, barcodes):
                wt_count += 1
            else:
                pass


    
    total_counts = sum(mutants_dict.values()) + k340h_count + wt_count

    with open(file_path, 'w') as csv_file:  
        writer = csv.writer(csv_file)
        writer.writerow(['OBD1 ID', 'OBD2 ID', 'OBD3 ID', 'Count', 'Frequency'])
        for key, value in mutants_dict.items():
            ids = key.split('_')
            ids = ['' if id == '0' else id for id in ids]
            writer.writerow([ids[0], ids[1], ids[2], value, (value/total_counts)])

        writer.writerow(['WT', '', '', wt_count, (wt_count/total_counts)]) # adds wildtype counts
        writer.writerow(['K340H', '', '', k340h_count, (k340h_count/total_counts)]) # adds k340h counts