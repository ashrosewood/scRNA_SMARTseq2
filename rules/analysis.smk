rule get_genes:
     input: rules.filter_counts.output
     output:
        "data/gene_hg19.cellRanger_metadata.tsv"
     shell:
        "Rscript scripts/getEnsemblGenes.R"

rule create_Seurat:
     input:
        "data/gene_hg19.cellRanger_metadata.tsv",
        "data/counts/raw_counts_.filt.tsv"
     output:
        "data/SeuratObject.rds"
     shell:
        "Rscript scripts/create_Seurat.R"