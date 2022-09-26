SAMPLES = ['OBDfinallibrary_replicate1_naive',
            'OBDfinallibrary_replicate1_selected',
            'OBDfinallibrary_replicate2_naive',
            'OBDfinallibrary_replicate2_selected',
            'OBDfinallibrary_replicate3_naive',
            'OBDfinallibrary_replicate3_selected']


def names_enrichment(names):
    return ['_'.join(name.split('_')[:2]) for name in names]

def names_summary(names):
    return [name.split('_')[0] for name in names]

#def replicates_present(names):
#    for name in names:
#        if 'replicate' not in name:
#            return "expand('data/results/enrichment_scores/{sample}_enrichment.csv', sample=names_enrichment(SAMPLES))"
#    return str(expand('data/results/summary_tables/{sample}_summary.csv', sample=names_summary(SAMPLES)))


rule all:
    input:
        expand('data/results/oligo_counts/{sample}_BC1_oligo_counts.csv', sample=SAMPLES),
        expand('data/results/oligo_counts/{sample}_BC2_oligo_counts.csv', sample=SAMPLES),
        expand('data/results/oligo_counts/{sample}_BC3_oligo_counts.csv', sample=SAMPLES),
        expand('data/results/oligo_mutant_counts/{sample}_oligo_mutants.csv', sample=SAMPLES),
#        expand('data/results/enrichment_scores/{sample}_enrichment.csv', sample=names_enrichment(SAMPLES))
#        replicates_present(SAMPLES)


rule filter_merge_sequences:
    input:
        R1 = 'data/raw_fastq_files/{sample}_R1.fastq.gz',
        R2 = 'data/raw_fastq_files/{sample}_R2.fastq.gz'
    params:
        qs='20'
    log:
        "data/logs/{sample}.html"
    output:
        R1 = "data/fastq_files_unpaired/{sample}_R1_unpaired.fastq.gz",
        R2 = "data/fastq_files_unpaired/{sample}_R2_unpaired.fastq.gz",
        M = "data/fastq_files_filtered_merged/{sample}_filtered_merged.fastq.gz"
    shell:
        "fastp --in1 {input.R1} --in2 {input.R2} --out1 {output.R1} --out2 {output.R2} "
        " --qualified_quality_phred {params.qs} -A --merge --merged_out {output.M} --html {log}"


rule get_barcode_counts:
    input:
        'data/fastq_files_filtered_merged/{sample}_filtered_merged.fastq.gz'
    output:
        'data/results/barcode_counts/{sample}_BC1_barcode_counts.csv',
        'data/results/barcode_counts/{sample}_BC2_barcode_counts.csv',
        'data/results/barcode_counts/{sample}_BC3_barcode_counts.csv',
        'data/results/barcode_mutant_counts/{sample}_barcode_mutants.csv'
    script:
        'scripts/barcode_analysis.py'


rule get_oligo_counts:
    input:
        'data/results/barcode_counts/{sample}_barcode_counts.csv'
    output:
        'data/results/oligo_counts/{sample}_oligo_counts.csv'
    script:
        'scripts/oligo_counts.py'


rule get_oligo_mutant_counts:
    input:
        'data/results/barcode_mutant_counts/{sample}_barcode_mutants.csv',
        'data/fastq_files_filtered_merged/{sample}_filtered_merged.fastq.gz'
    output:
        'data/results/oligo_mutant_counts/{sample}_oligo_mutants.csv'
    script:
        'scripts/oligo_mutants.py'


rule calculate_enrichments:
    input:
        'data/results/oligo_mutant_counts/{sample}_naive_oligo_mutants.csv',
        'data/results/oligo_mutant_counts/{sample}_selected_oligo_mutants.csv'
    output:
        'data/results/enrichment_scores/{sample}_enrichment.csv'
    script:
        'scripts/enrichment_scores.py'


rule summary:
    input:
        'data/results/enrichment_scores/{sample}_replicate1_enrichment.csv',
        'data/results/enrichment_scores/{sample}_replicate2_enrichment.csv',
        'data/results/enrichment_scores/{sample}_replicate3_enrichment.csv',

        'data/results/oligo_mutant_counts/{sample}_replicate1_naive_oligo_mutants.csv',
        'data/results/oligo_mutant_counts/{sample}_replicate2_naive_oligo_mutants.csv',
        'data/results/oligo_mutant_counts/{sample}_replicate3_naive_oligo_mutants.csv',

        'data/results/oligo_mutant_counts/{sample}_replicate1_selected_oligo_mutants.csv',
        'data/results/oligo_mutant_counts/{sample}_replicate2_selected_oligo_mutants.csv',
        'data/results/oligo_mutant_counts/{sample}_replicate3_selected_oligo_mutants.csv'
    output:
        'data/results/summary_tables/{sample}_summary.csv'
    script:
        'scripts/summary_table.py'