rule create_Seurat:
     input:
        "data/gene_metadata.tsv",
        "data/counts/raw_counts_.filt.tsv"
     output:
        "plots/seurat/pre_QCviolin.png",
        "plots/seurat/rna_libsize_barplot.pdf",
        "plots/seurat/rna_features_barplot.pdf",
        "plots/seurat/post_QCviolin.png",
        "plots/seurat/post_elbowPlot.png",       
        "plots/seurat/post_dimHeatPCA.png",
        "plots/seurat/post_SO_UMAP.png",
        "data/seurat/post_var_genes.tsv",
        "plots/seurat/post_var_genes_scatter.png",
        "data/seurat/post_DEgenes.tsv",
        "plots/seurat/post_DEheatmap.png",
        "data/seurat/SeuratObject.rds"
     params:
        coverage_threshold = config["coverage_threshold"],
        features_threshold = config["features_threshold"],
        top50_threshold = config["top50_threshold"],
        MT_threshold = config["MT_threshold"],
        Feature_lowerQuantile = config["Feature_lowerQuantile"],
        Feature_upperQuantile = config["Feature_upperQuantile"],
        percentMT_upperQuantile = config["percentMT_upperQuantile"]
     conda:
        "../envs/NMT_Create_Seurat.yaml"        	
     shell:
        "Rscript scripts/create_Seurat.R"