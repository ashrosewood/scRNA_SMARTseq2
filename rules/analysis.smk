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
        "data/SeuratObject.rds",
	"plots/SO_UMAP.png",
	"plots/var_genes_scatter.png",
	"plots/DE_heatmap.png",
	"tables/DE_genes.tsv",
	"tables/var_genes.tsv"
     shell:
        "Rscript scripts/create_Seurat.R"