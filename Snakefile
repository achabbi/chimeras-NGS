SAMPLES = ['OBD1', 'OBD12', 'OBD123A', 'OBD123B', 'OBD123repcapTA', 'OBD123repcapLW', 'OBD123finallibrary']

rule all:
    input:
        expand('data/results/barcode_counts/{sample}_BC1_counts.csv', sample=SAMPLES),
        expand('data/results/barcode_counts/{sample}_BC2_counts.csv', sample=SAMPLES),
        expand('data/results/barcode_counts/{sample}_BC3_counts.csv', sample=SAMPLES),
        expand('data/results/mutant_counts/{sample}_mutants_oligos.csv', sample=SAMPLES),
        expand('data/results/oligo_counts/{sample}_BC1_oligo_counts.csv', sample=SAMPLES),
        expand('data/results/oligo_counts/{sample}_BC2_oligo_counts.csv', sample=SAMPLES),
        expand('data/results/oligo_counts/{sample}_BC3_oligo_counts.csv', sample=SAMPLES)


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


rule get_barcode_sequence_counts:
    input:
        'data/fastq_files_filtered_merged/{sample}_filtered_merged.fastq.gz'
    output:
        'data/results/barcode_counts/{sample}_BC1_counts.csv',
        'data/results/barcode_counts/{sample}_BC2_counts.csv',
        'data/results/barcode_counts/{sample}_BC3_counts.csv',
        'data/results/sequence_counts/{sample}_mutants.csv'
    script:
        'scripts/analysis.py'


rule get_oligo_counts:
    input:
        'data/results/barcode_counts/{sample}_counts.csv'
    output:
        'data/results/oligo_counts/{sample}_oligo_counts.csv'
    script:
        'scripts/oligo_counts.py'


rule get_mutant_counts:
    input:
        'data/results/sequence_counts/{sample}_mutants.csv'
    output:
        'data/results/mutant_counts/{sample}_mutants_oligos.csv'
    script:
        'scripts/mutants.py'