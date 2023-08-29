configfile='config.yml'

rule trim_galore:
    input:
        'data/{sample}_1.fastq',
        'data/{sample}_2.fastq'
    output:
        directory('results/fastq_trimmed/{sample}')
    conda:
        'env/trim_galore.yml'
    log:
        'log/trim_galore/{sample}.log'
    threads: 
        config['cores']
    shell:
        ('trim_galore --paired {input} -q 20 --phred33 --fastqc --length 50 -o {output} -j {threads}) > {log}'

rule bbmap:
    input:
        'results/fastq_trimmed/{sample}'
    output:
        out1='results/fastQ_trimmed_norRNA/{sample}_1.fastq.gz',
        out2='results/fastQ_trimmed_norRNA/{sample}_2.fastq.gz',
        outm1='results/fastQ_trimmed_rRNA/{sample}_1.fastq.gz',
        outm2='results/fastQ_trimmed_rRNA/{sample}_2.fastq.gz'
    conda:
        'env/bbmap.yml'
    param:
        ribokmers='bin/ribokmers.fa',
        k='31'
    log:
        'log/bbmap/{sample}.log'
    shell:
        ('bbduk.sh in={input}_1_val_1.fq in2={input}_2_val_2.fq' 
        'out={output.out1} out2={output.out2}'
        'outm={output.outm1} outm2={output.outm2}'
        'k={param.k} ref={param.ribokmers}) > {log}

rule seqkit:
    input:
        norRNA1='results/fastQ_trimmed_norRNA/{sample}_1.fastq.gz',
        norRNA2='results/fastQ_trimmed_norRNA/{sample}_2.fastq.gz',
        rRNA1='results/fastQ_trimmed_rRNA/{sample}_1.fastq.gz',
        rRNA2='results/fastQ_trimmed_rRNA/{sample}_2.fastq.gz'
    output:
        norRNA='results/fastQ_trimmed_norRNA/seqkit_norRNA.txt',
        rRNA='results/fastQ_trimmed_rRNA/seqkit_rRNA.txt'
    conda:
        'env/seqkit.yml'
    log:
        'log/seqkit/{sample}.log
    shell:
        '(seqkit stats fastQ_trimmed_norRNA/*.fastq.gz > {output.norRNA} ;'
        'seqkit stats fastQ_trimmed_rRNA/*.fastq.gz > {output.rRNA}) > {log}'

rule bowtie:
    input:
        norRNA1='results/fastQ_trimmed_norRNA/{sample}_1.fastq.gz',
        norRNA2='results/fastQ_trimmed_norRNA/{sample}_2.fastq.gz'
    output:
        'results/Bowtie2_SAM/{sample}.sam'
    conda:
        'env/bowtie.yml'
    log:
        'log/bowtie/{sample}.log'
    threads: 
        config['threads']
    param: 
        ref=config['reference_genome'],
        strain=config['strain']
    shell:
    '(bowtie2-build --threads {threads} {param.ref} {param.strain} &&'
    'bowtie2 -x {param.strain} -1 {input.norRNA1} -2 {input.norRNA2}'
    ' -S {output} --no-mixed --threads {threads}) > {log}'

rule bam:
    input:
        sam='results/Bowtie2_SAM/{sample}.sam'
    output:
        bam='results/BAM/{sample}.bam'
    conda:
        'env/samtools.yml'
    log:
        'log/bam/{sample}.log'
    threads:
        config['threads']
    shell:
        'samtools view -bS {input.sam} -@ {threads} > {output.bam}'

rule sort:
    input:
    output:
    conda:
        'env/samtools.env'
    log:
        'log/sort/{sample}.log'
    shell: