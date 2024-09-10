configfile: 'config.yml'

rule all:
    input:
        'data/featureCounts/featureCounts_table.txt'
        
rule trim_galore:
    input:
        'data/fastq/{sample}_R1.fastq.gz',
        'data/fastq/{sample}_R2.fastq.gz'
    output:
        directory('data/fastq_trimmed/{sample}')
    conda:
        'env/trim-galore.yml'
    threads: 4
    shell:
        'trim_galore --paired {input} -q 20 --phred33 --fastqc --length 50 -o {output} -j {threads}'

rule bbmap:
    input:
        'data/fastq_trimmed/{sample}'
    output:
        out1='data/fastQ_trimmed_norRNA/{sample}_1.fastq.gz',
        out2='data/fastQ_trimmed_norRNA/{sample}_2.fastq.gz',
        outm1='data/fastQ_trimmed_rRNA/{sample}_1.fastq.gz',
        outm2='data/fastQ_trimmed_rRNA/{sample}_2.fastq.gz'
    conda:
        'env/bbmap.yml'
    params:
        ribokmers='bin/ribokmers.fa',
        k='31'
    threads: 8
    shell:
        'bbduk.sh in={input}/{wildcards.sample}_R1_val_1.fq.gz in2={input}/{wildcards.sample}_R2_val_2.fq.gz out={output.out1} out2={output.out2} outm={output.outm1} outm2={output.outm2} k={params.k} ref={params.ribokmers}'

rule seqkit:
    input:
        norRNA1='data/fastQ_trimmed_norRNA/{sample}_1.fastq.gz',
        norRNA2='data/fastQ_trimmed_norRNA/{sample}_2.fastq.gz',
        rRNA1='data/fastQ_trimmed_rRNA/{sample}_1.fastq.gz',
        rRNA2='data/fastQ_trimmed_rRNA/{sample}_2.fastq.gz'
    output:
        norRNA='data/fastQ_trimmed_norRNA/{sample}_seqkit_norRNA.txt',
        rRNA='data/fastQ_trimmed_rRNA/{sample}_seqkit_rRNA.txt'
    conda:
        'env/seqkit.yml'
    threads: 4
    shell:
        'seqkit stats fastQ_trimmed_norRNA/*.fastq.gz > {output.norRNA} ;'
        'seqkit stats fastQ_trimmed_rRNA/*.fastq.gz > {output.rRNA}'

rule bowtie:
    input:
        norRNA1='data/fastQ_trimmed_norRNA/{sample}_1.fastq.gz',
        norRNA2='data/fastQ_trimmed_norRNA/{sample}_2.fastq.gz'
    output:
        'data/Bowtie2_SAM/{sample}.sam'
    conda:
        'env/bowtie.yml'
    threads: 8
    params: 
        ref=config['reference_genome'],
        strain=config['strain']
    shell:
        'bowtie2-build --threads {threads} {params.ref} {params.strain} &&'
        'bowtie2 -x {params.strain} -1 {input.norRNA1} -2 {input.norRNA2} -S {output} --no-mixed --threads {threads} &&'
        'mv *.bt2 data/Bowtie2_SAM/'

rule bam:
    input:
        sam='data/Bowtie2_SAM/{sample}.sam'
    output:
        bam='data/BAM/{sample}.bam'
    conda:
        'env/samtools.yml'
    threads: 8 
    shell:
        'samtools view -bS {input.sam} -@ {threads} > {output.bam}'

rule sort:
    input:
        'data/BAM/{sample}.bam'
    output:
        'data/BAM_sorted/{sample}_sorted.bam'
    conda:
        'env/samtools.yml'
    threads: 4 
    shell:
        'samtools sort {input} -@ {threads} -o {output}'

rule feature_count:
    input:
        expand('data/BAM_sorted/{sample}_sorted.bam',sample=config['samples'])
    output:
        'data/featureCounts/featureCounts_table.txt'
    conda:
        'env/subread.yml'
    threads:
        config['threads']
    params:
        gtf=config['gtf']
    shell:
        'featureCounts -a {params.gtf} -p -T {threads} -t CDS -g gene_id -o {output} data/BAM_sorted/*.bam'