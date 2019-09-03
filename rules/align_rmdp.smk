rule trimming:
    input:
        "samples/raw/{sample}_R1.fastq.gz",
	"samples/raw/{sample}_R2.fastq.gz"
    output:
        "samples/trimmed/{sample}_R1_val_1.fq.gz",
	"samples/trimmed/{sample}_R2_val_2.fq.gz",
        "samples/fastqc/{sample}_R1_val_1_fastqc.zip",
	"samples/fastqc/{sample}_R2_val_2_fastqc.zip",
        "samples/fastqc/{sample}_R1_val_1_fastqc.html",
	"samples/fastqc/{sample}_R2_val_2_fastqc.html",
        "samples/trimmed/{sample}_R1.fastq.gz_trimming_report.txt",
	"samples/trimmed/{sample}_R2.fastq.gz_trimming_report.txt"
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
        "samples/hisat2/{sample}_output.sam"
    threads: 12
    params:
        gtf=config["gtf_file"]
    run:
        HiSat2=config["hisat2_tool"],
        pathToGenomeIndex = config["hisat2_index"]

        shell("""
                {HiSat2} -q -x {pathToGenomeIndex} \
                -1 {input[0]} -2 {input[1]} -p {threads} \
                --dta --sp 1000,1000 --no-mixed \
                --no-discordant -S {output}           
                """)


rule feature_count:
    input:
        expand("samples/hisat2/{sample}_output.sam", sample = SAMPLES)
    output:
        "data/{project_id}_counts.txt".format(project_id=config["project_id"]),
        "data/{project_id}_counts.txt.summary".format(project_id=config["project_id"])
    conda:
        "../envs/featureCounts.yaml"
    shell:
        """featureCounts -t exon -g gene_id -a /home/groups/CEDAR/anno/gtf/hg19_ens87.chr.gtf -o {output[0]} {input}"""


rule filter_counts:
    input:
        countsFile="data/{project_id}_counts.txt".format(project_id=config["project_id"])
    output:
        "data/{project_id}_counts.filt.txt".format(project_id=config["project_id"])
    params:
        anno=config["filter_anno"],
        biotypes=config["biotypes"],
        mito=config['mito']
    script:
        "../scripts/RNAseq_filterCounts.R"