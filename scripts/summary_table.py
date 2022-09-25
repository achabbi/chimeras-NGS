import csv

enrichment_file_rep1 = ''
enrichment_file_rep2 = ''
enrichment_file_rep3 = ''
enrichment_files = [enrichment_file_rep1, enrichment_file_rep2, enrichment_file_rep3]

mutants_naive_file_rep1 = ''
mutants_naive_file_rep2 = ''
mutants_naive_file_rep3 = ''
mutants_naive_files = [mutants_naive_file_rep1, mutants_naive_file_rep2, mutants_naive_file_rep3]

mutants_selected_file_rep1 = ''
mutants_selected_file_rep2 = ''
mutants_selected_file_rep3 = ''
mutants_selected_files = [mutants_selected_file_rep1, mutants_selected_file_rep2, mutants_selected_file_rep3]

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

# data is loaded
summary = []
ids_checked = []

for rep in enrichment_chimeras:
    for ids, score in rep.items():

        if ids in ids_checked:
            continue

        data_naive = {}
        data_selected = {}

        for i in range(len(enrichment_chimeras)):
            try:
                data_naive[i] = mutants_naive[i][ids]
            except:
                data_naive[i] = ''
        
            try:
                data_selected[i] = mutants_selected[i][ids]
            except:
                data_selected[i] = ''
        
        enrichment_score = score
        ids_checked.append(ids)



    

        


