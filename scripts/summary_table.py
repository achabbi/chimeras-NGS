import re
import csv

summary_file = snakemake.output[0]

num_reps = max([int(re.compile(r"replicate(\d+)").search(file).group(1)) for file in snakemake.input])

enrichment_files = [file for file in snakemake.input[:num_reps]]
mutants_naive_files = [file for file in snakemake.input[num_reps:2*num_reps]]
mutants_selected_files = [file for file in snakemake.input[2*num_reps:3*num_reps]]

# OBD1 ID, OBD2 ID, OBD3 ID, Naive Frequency Replicate 1, Naive Frequency Replicate 2, Naive Frequency Replicate 3, 
# Selected Frequency Replicate 1, Selected Frequency Replicate 2, Selected Frequency Replicate 3, Enrichment Score Replicate 1,
# Enrichment Score Replicate 2, Enrichment Score Replicate 3

enrichment_chimeras = [{}, {}, {}]
mutants_naive = [{}, {}, {}]
mutants_selected = [{}, {}, {}]

for i in range(len(enrichment_files)):
    with open(enrichment_files[i]) as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
                try:
                    ids = (int(row[0]), int(row[1]), int(row[2]))
                    enrichment_chimeras[i][ids] = int(row[3])
                except:
                    pass

    with open(mutants_naive_files[i]) as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
                try:
                    ids = (int(row[0]), int(row[1]), int(row[2]))
                    mutants_naive[i][ids] = int(row[4])
                except:
                    pass
    
    with open(mutants_selected_files[i]) as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
                try:
                    ids = (int(row[0]), int(row[1]), int(row[2]))
                    mutants_selected[i][ids] = int(row[4])
                except:
                    pass

summary = []
ids_checked = []

for rep in enrichment_chimeras:
    for ids in rep.keys():
        if ids in ids_checked:
            continue

        data_naive = {}
        data_selected = {}
        enrichment = {}

        for i in range(len(enrichment_chimeras)): # checks mutants in each replicate
            try:
                data_naive[i] = mutants_naive[i][ids]
            except:
                data_naive[i] = ''
        
            try:
                data_selected[i] = mutants_selected[i][ids]
            except:
                data_selected[i] = ''

            try:
                enrichment[i] = enrichment_chimeras[i][ids]
            except:
                enrichment[i] = ''
        
        ids_checked.append(ids)

        row = [id for id in ids]
        row.extend([data_naive[i] for i in range(len(data_naive))])
        row.extend([data_selected[i] for i in range(len(data_selected))])
        row.extend([enrichment[i] for i in range(len(enrichment))])
        summary.append(row)


with open(summary_file, 'w') as csv_file:  
        writer = csv.writer(csv_file)
        writer.writerow(['OBD1 ID', 'OBD2 ID', 'OBD3 ID', 
                        'Naive Frequency Replicate 1', 'Naive Frequency Replicate 2, Naive Frequency Replicate 3',
                        'Selected Frequency Replicate 1', 'Selected Frequency Replicate 2, Selected Frequency Replicate 3',
                        'Enrichment Score Replicate 1', 'Enrichment Score Replicate 2', 'Enrichment Score Replicate 3'])
        for row in summary:
            writer.writerow(row)

        


