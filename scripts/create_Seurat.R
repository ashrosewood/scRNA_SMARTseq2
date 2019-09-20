library(Seurat)
library(data.table)
library(purrr)
library(scater)
library(ggplot2)
library(cowplot)
library(png)

barplot_theme <- function() {
  p <- theme(
    plot.title = element_text(size=20, hjust=0.5),
    # axis.title.x = element_text(colour="black", size=25, vjust=1.5),
    axis.title.x = element_text(colour="black", size=18),
    axis.title.y = element_text(colour="black", size=18),
    # axis.text.x = element_text(colour="black",size=rel(1.6)),
    axis.text.y = element_text(colour="black",size=rel(1.5)),
    axis.line = element_line(colour="black", size=rel(0.7)),
   # axis.ticks.x = element_line(colour="black", size=rel(0.7)),
    axis.ticks.y = element_line(colour="black", size=rel(0.7)),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    #legend.position="none",
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank()
  )
}

fread_df <- partial(fread, data.table = FALSE)

## Options ##
opts <- list()

# Stringent thresholds
#opts$coverage_threshold <- 2e5    # Minimum library size (coverage)
#opts$features_threshold <- 1000   # Minimum number of expressed features
#opts$top50_threshold <- 0.75      # Maximum fraction of reads accounting for the top 50 features

# Lenient thresholds
opts$coverage_threshold <- 1e3    # Minimum library size (coverage)
opts$features_threshold <- 500   # Minimum number of expressed features
opts$top50_threshold <- 0.75      # Maximum fraction of reads accounting for the top 50 features

# Seurat Lenient Thresholds
opts$Feature_lowerQuantile <- .01
opts$Feature_upperQuantile <- .99
opts$Count_upperQuantile <- .99
opts$percentMT_upperQuantile <- .95


opts$MT_threshold <- 0.25         # Maximum fraction of reads mapping to mithocondrial genes

## I/O ##
io <- list()
io$in.gene_metadata <- "data/gene_hg19.cellRanger_metadata.tsv"
#io$in.sample_metadata <- "Box/My_NMT-seq_work/Stephen/code/MCF7_matched-selected/rna/counts_hg19/sample_metadata.tsv"
io$in.raw_counts <- "data/counts/raw_counts_.filt.tsv"
io$out.file <- "data/SeuratObject.rds"
#io$out.sample_metadata <- "Box/My_NMT-seq_work/Stephen/code/MCF7_matched-selected/sample_metadata.tsv"



## read in counts data and associated metadata ##

counts <- fread_df(io$in.raw_counts) %>% 
  tibble::column_to_rownames("ens_id") %>% 
  as.matrix()

feature_metadata <- fread_df(io$in.gene_metadata)
rownames(feature_metadata) <- feature_metadata$ens_id

genes <- rownames(feature_metadata[rownames(feature_metadata) %in% rownames(counts),])
feature_metadata <- feature_metadata[genes,]
counts <- counts[rownames(feature_metadata),]

## Create Seurat Object ##

SO <- CreateSeuratObject(counts = counts, min.cells = 5, min.features = opts$features_threshold)

## Calculate quality metrics ##

SO[["percent.mt"]] <- PercentageFeatureSet(object = SO, pattern = "^MT-")

png("plots/QC/pre_QCviolin.png")
VlnPlot(object=SO, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

## Filter by library size ##

SO <- subset(x = SO, subset = nFeature_RNA >= quantile(SO$nFeature_RNA, opts$Feature_lowerQuantile) & nFeature_RNA <= quantile(SO$nFeature_RNA, opts$Feature_upperQuantile) & nCount_RNA <= quantile(SO$nCount_RNA, opts$Count_upperQuantile) & nCount_RNA >= opts$coverage_threshold & percent.mt <= quantile(SO$percent.mt, opts$percentMT_upperQuantile))

png("plots/QC/post_QCviolin.png")
VlnPlot(object=SO, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

## Normalize and Scale Data ##

SO <- NormalizeData(SO)
SO <- FindVariableFeatures(object = SO, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(x = SO)
SO <- ScaleData(object = SO, features = all.genes)

## Save Seurat Object ##

saveRDS(SO, io$out.file)
