# chimeras-NGS
[Snakemake](https://snakemake.github.io/) pipeline to streamline analysis of NGS deep sequencing data on AAV protein chimeras.

## Setup
1. Run `git clone https://github.com/archit-c/chimeras-NGS.git` to clone repository.
2. Install [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html#regular-installation).
3. Run the following commands to create conda environment and install required packages:
```
$ conda create --name chimeras python=3.9
$ conda activate chimeras
$ conda install -c bioconda snakemake
$ conda install -c bioconda fastp
$ conda install -c conda-forge biopython
```

## Usage
### Copying files and naming conventions
1. Copy `fastq.gz` files to be processed into `data/raw_fastq_files`.
2. Rename `fastq.gz` files to match naming convention `[name]_R1.fastq.gz` & `[name]_R2.fastq.gz`, where `[name]` is `[identifier]_naive` or `[identifier]_selected`.
    - Example file pair: `OBD123A_naive_R1.fastq.gz` and `OBD123A_naive_R2.fastq.gz`.
3. Open `Snakefile`  and update `SAMPLES` list with the names of all files to be processed excluding `_[R1/R2].fastq.gz`.
    - Example: `SAMPLES = ['OBD12_naive', 'OBD123A_naive']` will run the pipeline on the following pairs of files from `data/raw_fastq_files`:
        - `OBD12_naive_R1.fastq.gz` and `OBD12_naive_R2.fastq.gz`
        - `OBD123A_naive_R1.fastq.gz` and `OBD123A_naive_R2.fastq.gz`
    - If samples are in replicates (pipeline is currently configured for 3 replicates per sample), then the following naming convention should be used: `SAMPLES = ['OBD12_replicate1_naive', 'OBD12_replicate2_naive', 'OBD12_replicate3_naive', 'OBD123A_replicate1_naive', 'OBD123A_replicate2_naive', 'OBD123A_replicate3_naive']`. The files in `data/raw_fastq_files` should be named accordingly:
        - `OBD12_replicate1_naive_R1.fastq.gz` and `OBD12_replicate1_naive_R2.fastq.gz`
        - `OBD12_replicate2_naive_R1.fastq.gz` and `OBD12_replicate2_naive_R2.fastq.gz`
        - `OBD12_replicate3_naive_R1.fastq.gz` and `OBD12_replicate3_naive_R2.fastq.gz`
        - `OBD123A_replicate1_naive_R1.fastq.gz` and `OBD123A_replicate1_naive_R2.fastq.gz`
        - `OBD123A_replicate2_naive_R1.fastq.gz` and `OBD123A_replicate2_naive_R2.fastq.gz`
        - `OBD123A_replicate3_naive_R1.fastq.gz` and `OBD123A_replicate3_naive_R2.fastq.gz`

### Running the pipeline
To run the pipeline from the home directory, run the following commands:
```
$ cd chimeras-NGS/
$ conda activate chimeras
$ snakemake --cores
```

## Results
All results are `csv` files contained in the `data/results` directory. Each directory in `results` contain files with different analysis results.
- `data/results/barcode_counts`: counts of individual barcodes from each region.
- `data/results/oligo_counts`: counts of individual oligos from each region.
- `data/results/barcode_mutant_counts`: counts of chimeric barcode mutants.
- `data/results/oligo_mutant_counts`: counts of chimeric oligo mutants.
- `data/results/enrichment_scores`: enrichment scores of chimeric oligo mutants.
    - **NOTE**: both naive and selected files for each sample must be present in `data/raw_fastq_files` for enrichment results.
- `data/results/summary_tables`: merged file with enrichment and frequency data for each oligo mutant in each replicate.
