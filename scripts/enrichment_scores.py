import csv
import math

def L(counts_naive, counts_selected):
    """
    Calculates unscaled enrichment score for variant.

    Inputs:
        - counts_naive: variant counts in naive sample
        - counts_selected: variant counts in selected sample
        - wt_naive_counts: wildtype counts in naive sample
        - wt_selected_counts: wildtype counts in selected sample
    
    Output: returns unscaled enrichment score
    """

#    t1 = math.log((counts_selected+0.5)/(wt_selected_counts+0.5), 2)
#    t2 = math.log((counts_naive+0.5)/(wt_naive_counts+0.5), 2)
#    return t1 - t2
    return math.log(counts_selected/counts_naive, 2)
    
    
if __name__ == '__main__':
    naive_mutants_file = snakemake.input[0]
    selected_mutants_file = snakemake.input[1]
    newfile = snakemake.output[0]

    naive = {}
    selected = {}
    enrichments = {}

    with open(naive_mutants_file) as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            try:
                naive[(int(row[0]), int(row[1]), int(row[2]))] = float(row[4]) # assigns id combination to its count
            except:
                pass
            
            if row[0] == 'WT':
                wt_naive_counts = float(row[4])
            elif row[0] == 'K340H':
                k340h_naive_counts = float(row[4])

        
    with open(selected_mutants_file) as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            try:
                ids = (int(row[0]), int(row[1]), int(row[2]))
                if ids in naive.keys():
                    selected[ids] = float(row[4])
            except:
                pass
        
            if row[0] == 'WT':
                wt_selected_counts = float(row[4])
            elif row[0] == 'K340H':
                k340h_selected_counts = float(row[4])
    
#    L_wt = L(wt_naive_counts, wt_selected_counts, wt_naive_counts, wt_selected_counts)
#    L_k340h = L(k340h_naive_counts, k340h_selected_counts, wt_naive_counts, wt_selected_counts)

    L_wt = L(wt_naive_counts, wt_selected_counts)
    L_k340h = L(k340h_naive_counts, k340h_selected_counts)

    for ids, count in selected.items():
        counts_selected = count
        counts_naive = naive[ids]

#        L_v = L(counts_naive, counts_selected, wt_naive_counts, wt_selected_counts)
#        L_vscaled = (L_v - L_k340h)/(L_wt - L_k340h)

        L_v = L(counts_naive, counts_selected)
        enrichments[ids] = L_v


    with open(str(newfile), 'w') as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(['OBD1 ID', 'OBD2 ID', 'OBD3 ID', 'Enrichment'])
        for ids, score in enrichments.items():
            writer.writerow([ids[0], ids[1], ids[2], score])

    


        


