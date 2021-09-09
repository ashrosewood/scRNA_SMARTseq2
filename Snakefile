__author__ = "Aaron Doe"
__email__ = "doe@ohsu.edu"
__license__ = "MIT"

"""Computation Hub omic data processing pipeline for SE Bulk RNA-sequencing"""

import datetime
import sys
import os
import pandas as pd
import json

timestamp = ('{:%Y-%m-%d_%H:%M:%S}'.format(datetime.datetime.now()))

configfile:"omic_config.yaml"
project_id = config["project_id"]

SAMPLES, = glob_wildcards("samples/raw/{sample}_R1.fq.gz")

ext = ['r','pdf','xls']
fastqscreen_ext = ['html','png','txt']
fastqc_ext = ['html','zip']

with open('cluster.json') as json_file:
     json_dict = json.load(json_file)

rule_dirs = list(json_dict.keys())

for rule in rule_dirs:
    if not os.path.exists(os.path.join(os.getcwd(),'logs',rule)):
       log_out = os.path.join(os.getcwd(), 'logs', rule)
       os.makedirs(log_out)
       print(log_out)


result_dirs = ['diffexp','tables']
for rule in result_dirs:
    if not os.path.exists(os.path.join(os.getcwd(),'results',rule)):
       log_out = os.path.join(os.getcwd(), 'results', rule)
       os.makedirs(log_out)
       print(log_out)


def message(mes):
    sys.stderr.write("|--- " + mes + "\n")


for sample in SAMPLES:
    message("Sample " + sample + " will be processed")

rule all:
    input:
#       expand("samples/trimmed/{sample}_R1.fq.gz_trimming_report.txt", sample = SAMPLES),
#       expand("samples/trimmed/{sample}_R2.fq.gz_trimming_report.txt", sample = SAMPLES),
       expand("samples/fastqc/{sample}_R1_val_1_fastqc.{fastqc_ext}", sample = SAMPLES, fastqc_ext=fastqc_ext),
       expand("samples/fastqc/{sample}_R2_val_2_fastqc.{fastqc_ext}", sample = SAMPLES, fastqc_ext=fastqc_ext),
       expand("samples/fastqscreen/{sample}/{sample}_R1_val_1_screen.{fastqscreen_ext}", sample=SAMPLES, fastqscreen_ext=fastqscreen_ext),
       expand("samples/fastqscreen/{sample}/{sample}_R2_val_2_screen.{fastqscreen_ext}", sample=SAMPLES, fastqscreen_ext=fastqscreen_ext),
       expand("samples/hisat2/{sample}_output.bam", sample = SAMPLES),
#       expand("samples/hisat2/{sample}.saturation.pdf", sample = SAMPLES),
#       "donefile.txt",
       "samples/hisat2/all_merged.saturation.r",
       "data/counts/raw_counts_.tsv",
       "data/counts/raw_counts_.tsv.summary",       
       "data/counts/raw_counts_.filt.tsv",
       "data/counts/sample_metadata.tsv",
       "data/gene_metadata.tsv",
       "plots/seurat/pre_QCviolin.png",
       "plots/seurat/post_QCviolin.png",
       "plots/seurat/post_elbowPlot.png",       
       "plots/seurat/post_dimHeatPCA.png",
       "plots/seurat/post_SO_UMAP.png",
       "data/seurat/post_var_genes.tsv",
       "plots/seurat/post_var_genes_scatter.png",
       "data/seurat/post_DEgenes.tsv",
       "plots/seurat/post_DEheatmap.png",
       "data/seurat/SeuratObject.rds"

include: "rules/align_rmdp.smk"
include: "rules/analysis.smk"
