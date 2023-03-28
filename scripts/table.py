import csv

num_reps = 3
#naive_files = [file for file in snakemake.input[0:num_reps]]
#selected_files = [file for file in snakemake.input[num_reps:2*num_reps]]
#table_file = snakemake.output[0]
naive_samples = ['N1', 'N2', 'N3']
selected_samples = ['1i', '2i', '3i']

path = 'data/results/oligo_mutant_counts/'
naive_files = [''.join([path, f"{name}_oligo_mutants.csv"]) for name in naive_samples]
selected_files = [''.join([path, f"{name}_oligo_mutants.csv"]) for name in selected_samples]
table_file = 'data/results/table.csv'

naive = [{} for _ in range(num_reps)]
selected = [{} for _ in range(num_reps)]

for i in range(len(naive_files)):
    with(open(naive_files[i])) as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            if 'OBD' in row[0]:
                continue
            else:
                ids = (row[0], row[1], row[2])
                freq = row[4]
                naive[i][ids] = freq

for i in range(len(selected_files)):
    with(open(selected_files[i])) as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            if 'OBD' in row[0]:
                continue
            else:
                ids = (row[0], row[1], row[2])
                freq = row[4]
                selected[i][ids] = freq

data = {}
for i in range(len(naive)):
    rep = naive[i]
    for id, freq in rep.items():
        if id not in data.keys():
            data[id] = [0 for _ in range(num_reps*2)]
        data[id][i] = freq

for i in range(len(selected)):
    rep = selected[i]
    for id, freq in rep.items():
        if id not in data.keys():
            data[id] = [0 for _ in range(num_reps*2)]
        data[id][i+num_reps] = freq

rows = []
for key, value in data.items():
    row = [id for id in key]
    row.extend(value)
    rows.append(row)

with open(table_file, 'w') as csv_file:  
        writer = csv.writer(csv_file)
        writer.writerow(['OBD1', 'OBD2', 'OBD3', 
                        'Frequency Naive Replicate 1', 'Frequency Naive Replicate 2', 'Frequency Naive Replicate 3',
                        'Frequency Selected Replicate 1', 'Frequency Selected Replicate 2', 'Frequency Selected Replicate 3'])
        for row in rows:
            writer.writerow(row)