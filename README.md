# RNAseq Snakemake Pipeline

## How to use
1. Install conda or mamba
2. Install snakemake in it's own conda env
3. Create `data/fastq` directory and place raw RNAseq reads in the following format: Forward: `sample1_R1.fastq.gz` ; Reverse: `sample1_R2.fastq.gz`
4. Download the appropriate `reference.fna` and `reference.gtf` and place them in `bin` (some are already provided)
5. Open `config.yml` and adjust the values accordingly
6. Activate the snakemake conda env and run `snakemake -c{thread_number} --use-conda --conda-frontend {conda_client}
7. All results will be output to `data/`