rule create_Seurat:
     input:
        "data/gene_metadata.tsv",
        "data/counts/raw_counts_.filt.tsv"
     output:
        "plots/seurat/pre_QCviolin.png",
        "plots/seurat/post_QCviolin.png",
        "plots/seurat/post_elbowPlot.png",
	"plots/seurat/post_elbowPlot_CCreduced.png",
        "plots/seurat/post_dimHeatPCA.png",
        "plots/seurat/post_SO_UMAP.png",
	"plots/seurat/post_SO_UMAP_CCreduced.png",
	"plots/seurat/CCPhase_UMAP.png",
	"plots/seurat/CCPhase_UMAP_CCreduced.png",
        "data/seurat/post_var_genes.tsv",
        "plots/seurat/post_var_genes_scatter.png",
	"plots/seurat/post_var_genes_scatter_CCreduced.png",
        "data/seurat/post_DEgenes.tsv",
	"data/seurat/post_DEgenes_CCreduced.tsv",
        "plots/seurat/post_DEheatmap.png",
	"plots/seurat/post_DEheatmap_CCreduced.png",
        "data/seurat/SeuratObject.rds",
        "data/seurat/SeuratObject_CCred.rds"
     params:
        coverage_threshold = config["coverage_threshold"],
        features_threshold = config["features_threshold"],
        top50_threshold = config["top50_threshold"],
        Count_upperQuantile = config["Count_upperQuantile"],
        Feature_lowerQuantile = config["Feature_lowerQuantile"],
        Feature_upperQuantile = config["Feature_upperQuantile"],
        percentMT_upperQuantile = config["percentMT_upperQuantile"],
	integrateTF = config["integrateTF"],
	clusterRes = config["clusterRes"]
     conda:
        "../envs/NMT_Create_Seurat.yaml"        	
     shell:
        "Rscript scripts/create_Seurat.R --covThresh={params.coverage_threshold} --featThresh={params.features_threshold} --top50={params.top50_threshold} --featureLQ={params.Feature_lowerQuantile} --featureUQ={params.Feature_upperQuantile} --CountUQ={params.Count_upperQuantile} --MTUQ={params.percentMT_upperQuantile} --integrate={params.integrateTF} --res={params.clusterRes}"
