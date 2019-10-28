conda create -n smartseq2 r-seurat bioconductor-singlecellexperiment r-data.table bioconductor-rsubread r-purrr bioconductor-rtracklayer r-ggplot2 r-cowplot bioconductor-scater r-png

conda env export > envs/scAnalysis.yaml
