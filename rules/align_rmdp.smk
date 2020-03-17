rule trimming:
    input:
        "samples/raw/{sample}_R1.fq.gz",
	"samples/raw/{sample}_R2.fq.gz"
    output:
        "samples/trimmed/{sample}_R1_val_1.fq.gz",
	"samples/trimmed/{sample}_R2_val_2.fq.gz",
        "samples/fastqc/{sample}_R1_val_1_fastqc.zip",
	"samples/fastqc/{sample}_R2_val_2_fastqc.zip",
        "samples/fastqc/{sample}_R1_val_1_fastqc.html",
	"samples/fastqc/{sample}_R2_val_2_fastqc.html",
        "samples/trimmed/{sample}_R1.fq.gz_trimming_report.txt",
	"samples/trimmed/{sample}_R2.fq.gz_trimming_report.txt"
    params:
        adapter = config["adapter"]
    conda:
        "../envs/trimG.yaml"
    message:
        """--- Trimming."""
    shell:
        """trim_galore --gzip --paired -o samples/trimmed/ --fastqc_args "--outdir samples/fastqc/" {input[0]} {input[1]}"""

rule fastqscreen:
    input:
        "samples/trimmed/{sample}_R1_val_1.fq.gz",
	"samples/trimmed/{sample}_R2_val_2.fq.gz"
    output:
        "samples/fastqscreen/{sample}/{sample}_R1_val_1_screen.html",
	"samples/fastqscreen/{sample}/{sample}_R2_val_2_screen.html",
        "samples/fastqscreen/{sample}/{sample}_R1_val_1_screen.png",
	"samples/fastqscreen/{sample}/{sample}_R2_val_2_screen.png",
        "samples/fastqscreen/{sample}/{sample}_R1_val_1_screen.txt",
	"samples/fastqscreen/{sample}/{sample}_R2_val_2_screen.txt"
    params:
        conf = config["conf"]
    conda:
        "../envs/fastqscreen.yaml"
    shell:
        """fastq_screen --aligner bowtie2 --conf {params.conf} --outdir samples/fastqscreen/{wildcards.sample} {input[0]} {input[1]}"""

rule Hisat2:
    input:
        "samples/trimmed/{sample}_R1_val_1.fq.gz",
	"samples/trimmed/{sample}_R2_val_2.fq.gz"
    output:
        "samples/hisat2/{sample}_output.bam"
    threads: 12
    params:
        gtf=config["gtf_file"],
	pathToGenomeIndex=config["hisat2_index"]
    conda:
        "../envs/Mapping.yaml"
    shell:
        """
	hisat2 -q -x {params.pathToGenomeIndex} -1 {input[0]} -2 {input[1]} -p {threads} --dta --sp 1000,1000 --no-mixed --no-discordant -S samples/hisat2/{wildcards.sample}_output.sam 
	samtools view -S -b samples/hisat2/{wildcards.sample}_output.sam > {output}         
        rm samples/hisat2/{wildcards.sample}_output.sam
        """

rule feature_count:
    input:
        expand("samples/hisat2/{sample}_output.bam", sample = SAMPLES)
    output:
        "data/counts/raw_counts_.tsv",
        "data/counts/raw_counts_.tsv.summary"
    params:
        anno=config["gtf_file"],	
    conda:
        "../envs/featureCounts.yaml"
    shell:
        """featureCounts -p -t exon -g gene_id -a {params.anno} -o {output[0]} {input}"""
	
rule filter_counts:
    input:
        countsFile="data/counts/raw_counts_.tsv"
    output:
        "data/counts/raw_counts_.filt.tsv",
	"data/counts/sample_metadata.tsv",
	"data/gene_metadata.tsv"	
    params:
        anno=config["filter_anno"],
        biotypes=config["biotypes"],
        mito=config['mito']
    script:
        "../scripts/RNAseq_filterCounts.R"
