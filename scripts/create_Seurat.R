args <- commandArgs()

help <- function(){
    cat("accessibility_profiles.R :
Make profiles plots around promoters. Output will be save in the parsed directory in the covPath. The promter and body regions (gene body w#ithout promoter) will be saved in the outDir.\n")
    cat("Usage: \n")
    cat("--covThresh   : Minimum library size                                        [required]\n")
    cat("--featThresh  : Minimum number of expressed features                        [required]\n")
    cat("--top50       : Maximum fraction of reads accounting for top 50 features    [required]\n")
    cat("--featureLQ   : Feature lower quantile                                      [required]\n")
    cat("--featureUQ   : Feature upper quantile                                      [required]\n")
    cat("--CountUQ     : Count upper quantile                                        [required]\n")
    cat("--MTUQ        : percent MT upper quantile                                   [required]\n")
    cat("--integrate   : TRUE or FALSE determining whether to integrate              [required]\n")
    cat("--res         : resolution for clustering                                   [required]\n")
    cat("\n")
    q()
}

io   <- list()
opts <- list()

# Save values of each argument
if( !is.na(charmatch("--help",args)) || !is.na(charmatch("--help",args)) ){
    help()
} else {
    coverage_threshold      <- sub( '--covThresh=', '', args[grep('--covThresh=', args)] )
    features_threshold      <- sub( '--featThresh=', '', args[grep('--featThresh=', args)] )
    top50_threshold         <- sub( '--top50=', '', args[grep('--top50=', args)] )
    Feature_lowerQuantile   <- sub( '--featureLQ=', '', args[grep('--featureLQ=', args)])
    Feature_upperQuantile   <- sub( '--featureUQ=', '', args[grep('--featureUQ=', args)])
    Count_upperQuantile     <- sub( '--CountUQ=', '', args[grep('--CountUQ=', args)])
    percentMT_upperQuantile <- sub( '--MTUQ=', '',args[grep('--MTUQ=',args)])
    integrateTF             <- sub( '--integrate=', '',args[grep('--integrate=',args)])
    res                     <- sub( '--res=', '', args[grep('--res=',args)])

}




##coverage_threshold = snakemake@params[['coverage_threshold']]

#setwd("../")
#coverage_threshold <- 1e5    # Minimum library size (coverage)
#features_threshold = snakemake@params[['features_threshold']]
#features_threshold <- 1000   # Minimum number of expressed features

#top50_threshold = snakemake@params[['top50_threshold']]
#top50_threshold <- 0.75      # Maximum fraction of reads accounting for the top 50 features
#Feature_lowerQuantile = snakemake@params[['Feature_lowerQuantile']]
#Feature_lowerQuantile <- .01
#Feature_upperQuantile = snakemake@params[['Feature_upperQuantile']]
#Feature_upperQuantile <- .99
#Count_upperQuantile = snakemake@params[['Count_upperQuantile']]
#Count_upperQuantile <- .99
#percentMT_upperQuantile = snakemake@params[['percentMT_upperQuantile']]
#percentMT_upperQuantile <- .85
#integrateTF = snakemake@params[['integrateTF']]
#integrateTF = FALSE
#res = 1.0

opts <- list()
opts$coverage_threshold <- as.numeric(coverage_threshold)

opts$features_threshold <- as.numeric(features_threshold)

opts$top50_threshold <- as.numeric(top50_threshold)

opts$Feature_lowerQuantile <- as.numeric(Feature_lowerQuantile)

opts$Feature_upperQuantile <- as.numeric(Feature_upperQuantile)

opts$Count_upperQuantile <- as.numeric(Count_upperQuantile)

opts$percentMT_upperQuantile <- as.numeric(percentMT_upperQuantile)

opts$integrate <- integrateTF

opts$res <- as.numeric(res)

## I/O ##
io <- list()
io$in.gene_metadata   <- "data/gene_metadata.tsv"
io$in.sample_metadata <- "data/counts/sample_metadata.tsv"
io$in.raw_counts      <- "data/counts/raw_counts_.filt.tsv"
io$out.file           <- "data/seurat/SeuratObject.rds"
io$out.file.CCred     <- "data/seurat/SeuratObject_CCred.rds"
io$dataDir            <- "data/seurat"
io$plotDir            <- "plots/seurat"

if( !( file.exists(io$dataDir)) )  {
    dir.create( io$dataDir, FALSE, TRUE )  
}

if(!( file.exists( io$plotDir ) ) ) {
    dir.create( io$plotDir, FALSE, TRUE )  
}

library(Seurat)
library(data.table)
library(purrr)
library(scater)
library(ggplot2)
library(cowplot)
library(png)
library(dplyr)
library(stringr)

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

## read in counts data and associated metadata ##
counts <- fread_df(io$in.raw_counts) %>% 
    tibble::column_to_rownames("ens_id")
counts$Genes <- NULL

feature_metadata           <- fread_df(io$in.gene_metadata)
rownames(feature_metadata) <- feature_metadata$ens_id
stopifnot(rownames(counts) == feature_metadata$ens_id)

counts           <- as.matrix(counts)
rownames(counts) <- feature_metadata$gene

# should  all be equal now since made from the same table
#if( ! all.equal(rownames(counts), rownames(feature_metadata)) ){
#    print("rownames in metadata not in the same order as counts")
#    genes            <- rownames(feature_metadata[rownames(feature_metadata) %in% rownames(counts),])
#    feature_metadata <- feature_metadata[genes,]
#    counts           <- counts[rownames(feature_metadata),]
#}

## Create Seurat Object ##
#rownames(counts) <- feature_metadata[rownames(counts),3]
SO <- CreateSeuratObject(counts = counts, min.cells = 0,
                         min.features = 0)

## save unique gene ids  
feature_metadata$gene_unique <- rownames(SO)

ctr <- 0
SO@meta.data$origin <- str_extract(rownames(SO@meta.data), "[^_]+")

SO[["percent.mt"]] <- PercentageFeatureSet(object = SO, pattern = "^MT-")
#SO[["percent.mt"]] <- PercentageFeatureSet(object = SO, features = mt)

png(paste0(io$plotDir, "/pre_QCviolin.png"))
VlnPlot(object=SO, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

filtSO.list <- list()

if (opts$integrate) {
    SO.list <- SplitObject(SO, split.by = "origin")
    for (subSO in SO.list) {
        ctr <- ctr + 1

############################
## Filter by library size ##
############################
        libsize.drop <- subSO$nCount_RNA < opts$coverage_threshold
        sum(libsize.drop)

        subSO$sample <- rownames(subSO[[]])

        libsize.drop_dt <- data.table(
            sample=subSO$sample, 
            size=subSO$nCount_RNA, 
            color=c("black","red")[as.numeric(libsize.drop)+1]
        ) %>% setkey(size) %>% .[,col:=size] %>% .[,sample:=factor(sample,levels=sample)]

        p1 <- ggplot(libsize.drop_dt, aes(x=sample, y=size)) +
            geom_bar(stat='identity', position="dodge", fill="#3CB54E") +
            geom_hline(yintercept=opts$coverage_threshold, colour="black", linetype="dashed") +
            scale_fill_gradient(low="red", high="green") +
            labs(y="Library size") +
                                        #barplot_theme() +
            xlab(paste("Threshold", opts$coverage_threshold))+
            theme(
                legend.position = "none",
                axis.title.x = element_text(size=rel(1.8)),
                axis.title.y = element_text(size=rel(1.8)),
                                        # axis.text.x = element_text(colour="black", color=foo$color, angle=90, size=10, vjust=0.5, hjust=1.0)
                axis.text.x = element_blank(),
                axis.ticks.x = element_blank()
            )
        print(p1)
        save_plot(paste0(io$plotDir, "/rna_libsize_barplot", as.character(ctr), ".pdf"), p1)

######################################
## Filter by number of expressed genes
######################################
        feature.drop <- subSO$nFeature_RNA < opts$features_threshold
        sum(feature.drop)

        feature.drop_dt <- data.table(
            sample=subSO$sample,
            features = subSO$nFeature_RNA, 
            color = c("black","red")[as.numeric(feature.drop)+1]
        ) %>% setkey(features) %>% .[,col:=features] %>% .[,sample:=factor(sample,levels=sample)]

        p2 <- ggplot(feature.drop_dt, aes(x=sample, y=features)) +
            geom_bar(stat='identity', position="dodge", fill="#3CB54E") +
            geom_hline(yintercept=opts$features_threshold, colour="black", linetype="dashed") +
                                        # scale_fill_gradient(low="red", high="green") +
            labs(y="Expressed genes") +
                                        #barplot_theme() +
            xlab(paste("Threshold", opts$coverage_threshold))+
            theme(
                legend.position = "none",
                axis.title.x = element_text(size=rel(1.8)),
                axis.title.y = element_text(size=rel(1.8)),
                                        # axis.text.x = element_text(colour="black", color=foo$color, angle=90, size=10, vjust=0.5, hjust=1.0)
                axis.text.x = element_blank(),
                axis.ticks.x = element_blank()
            )
        print(p2)
        save_plot(paste0(io$plotDir, "/rna_features_barplot", as.character(ctr), ".pdf"), p2)

#########
## Subset
#########

        subSO <- subset(x = subSO,
             subset = nFeature_RNA >= quantile(subSO$nFeature_RNA,
                                               opts$Feature_lowerQuantile)
             & nFeature_RNA <= quantile(subSO$nFeature_RNA, opts$Feature_upperQuantile)
             & nCount_RNA <= quantile(subSO$nCount_RNA, opts$Count_upperQuantile)
             & nCount_RNA >= opts$coverage_threshold
             & percent.mt <= quantile(subSO$percent.mt, opts$percentMT_upperQuantile))



#SO <- subset(x = SO, subset = nFeature_RNA >= quantile(SO$nFeature_RNA, opts$Feature_lowerQuantile) & nFeature_RNA <= quantile(SO$nFeature_RNA, opts$Feature_upperQuantile) & nCount_RNA <= quantile(SO$nCount_RNA, opts$Count_upperQuantile) & nCount_RNA >= opts$coverage_threshold & percent.mt <= quantile(SO$percent.mt, opts$percentMT_upperQuantile))


##############################
## Normalize and Scale Data ##
##############################
## need to add anchoring approach when we have multiple samples
        subSO        <- NormalizeData(subSO)
        subSO        <- FindVariableFeatures(object = subSO, selection.method = "vst", nfeatures = 2000)
        filtSO.list <- c(filtSO.list, subSO)
    }

} else {    
    ## Calculate quality metrics ##
    #mt <- feature_metadata[grep("^MT-", feature_metadata$gene), "ens_id"]
    ## Calculate quality metrics ##
    SO[["percent.mt"]] <- PercentageFeatureSet(object = SO, pattern = "^MT-")
    #SO[["percent.mt"]] <- PercentageFeatureSet(object = SO, features = mt)

    png(sub("$", "/pre_QCviolin.png", io$plotDir))
    VlnPlot(object=SO, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
    dev.off()

############################
## Filter by library size ##
############################
    libsize.drop <- SO$nCount_RNA < opts$coverage_threshold
    sum(libsize.drop)

    SO$sample <- rownames(SO[[]])

    libsize.drop_dt <- data.table(
        sample=SO$sample, 
        size=SO$nCount_RNA, 
        color=c("black","red")[as.numeric(libsize.drop)+1]
    ) %>% setkey(size) %>% .[,col:=size] %>% .[,sample:=factor(sample,levels=sample)]

    p1 <- ggplot(libsize.drop_dt, aes(x=sample, y=size)) +
        geom_bar(stat='identity', position="dodge", fill="#3CB54E") +
        geom_hline(yintercept=opts$coverage_threshold, colour="black", linetype="dashed") +
        scale_fill_gradient(low="red", high="green") +
        labs(y="Library size") +
                                        #barplot_theme() +
        xlab(paste("Threshold", opts$coverage_threshold))+
        theme(
            legend.position = "none",
            axis.title.x = element_text(size=rel(1.8)),
            axis.title.y = element_text(size=rel(1.8)),
                                        # axis.text.x = element_text(colour="black", color=foo$color, angle=90, size=10, vjust=0.5, hjust=1.0)
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank()
        )
    print(p1)
    save_plot(sub("$", "/rna_libsize_barplot.pdf", io$plotDir), p1)

######################################
## Filter by number of expressed genes
######################################
    feature.drop <- SO$nFeature_RNA < opts$features_threshold
    sum(feature.drop)

    feature.drop_dt <- data.table(
        sample=SO$sample,
        features = SO$nFeature_RNA, 
        color = c("black","red")[as.numeric(feature.drop)+1]
    ) %>% setkey(features) %>% .[,col:=features] %>% .[,sample:=factor(sample,levels=sample)]

    p2 <- ggplot(feature.drop_dt, aes(x=sample, y=features)) +
        geom_bar(stat='identity', position="dodge", fill="#3CB54E") +
        geom_hline(yintercept=opts$features_threshold, colour="black", linetype="dashed") +
                                        # scale_fill_gradient(low="red", high="green") +
        labs(y="Expressed genes") +
                                        #barplot_theme() +
        xlab(paste("Threshold", opts$coverage_threshold))+
        theme(
            legend.position = "none",
            axis.title.x = element_text(size=rel(1.8)),
            axis.title.y = element_text(size=rel(1.8)),
                                        # axis.text.x = element_text(colour="black", color=foo$color, angle=90, size=10, vjust=0.5, hjust=1.0)
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank()
        )
    print(p2)
    save_plot(sub("$", "/rna_features_barplot.pdf", io$plotDir), p2)

#########
## Subset
#########

    SO <- subset(x = SO,
             subset = nFeature_RNA >= quantile(SO$nFeature_RNA,
                                               opts$Feature_lowerQuantile)
             & nFeature_RNA <= quantile(SO$nFeature_RNA, opts$Feature_upperQuantile)
             & nCount_RNA <= quantile(SO$nCount_RNA, opts$Count_upperQuantile)
             & nCount_RNA >= opts$coverage_threshold
             & percent.mt <= quantile(SO$percent.mt, opts$percentMT_upperQuantile))



#SO <- subset(x = SO, subset = nFeature_RNA >= quantile(SO$nFeature_RNA, opts$Feature_lowerQuantile) & nFeature_RNA <= quantile(SO$nFeature_RNA, opts$Feature_upperQuantile) & nCount_RNA <= quantile(SO$nCount_RNA, opts$Count_upperQuantile) & nCount_RNA >= opts$coverage_threshold & percent.mt <= quantile(SO$percent.mt, opts$percentMT_upperQuantile))

    png(sub("$", "/post_QCviolin.png", io$plotDir))
    VlnPlot(object=SO, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
    dev.off()

##############################
## Normalize and Scale Data ##
##############################
## need to add anchoring approach when we have multiple samples
    SO        <- NormalizeData(SO)
    SO        <- FindVariableFeatures(object = SO, selection.method = "vst", nfeatures = 2000)
}

### Integrating (if necessary) ###

if (opts$integrate) {
    for (i in 1:length(filtSO.list)) {
        filtSO.list[[i]] <- NormalizeData(filtSO.list[[i]], verbose = FALSE)
        filtSO.list[[i]] <- FindVariableFeatures(filtSO.list[[i]], selection.method = "vst", 
                                             nfeatures = 2000, verbose = FALSE)
    }
    SO.anchors <- FindIntegrationAnchors(object.list = filtSO.list, dims = 1:30, k.filter=50)
    SO.integrated <- IntegrateData(anchorset = SO.anchors, dims = 1:30)
    SO <- SO.integrated
}

png(paste0(io$plotDir, "/post_QCviolin.png"))
VlnPlot(object=SO, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

## cell cycle scoring
cc_genes  <- read.table("data/regev_lab_cell_cycle_genes.txt")
s.genes   <- cc_genes$V1[1:43]
g2m.genes <- cc_genes$V1[44:97]

#s.genes.ens   <- feature_metadata[feature_metadata$gene %in% s.genes,"ens_id"]
#g2m.genes.ens <- feature_metadata[feature_metadata$gene %in% g2m.genes,"ens_id"]


SO        <- CellCycleScoring(SO, s.features = s.genes, g2m.features = g2m.genes)
CC_SO        <- ScaleData(object = SO, features = rownames(SO), vars.to.regress=c("S.Score","G2M.Score"))

SO        <- ScaleData(object = SO, features = rownames(SO))

CC_SO     <- RunPCA(CC_SO, approx = FALSE, verbose = FALSE)
SO        <- RunPCA(SO, approx = FALSE, verbose = FALSE)

# Examine and visualize PCA results a few different ways
                                        #print(x = SO[["pca"]], dims = 1:5, nfeatures = 5)
if(dim(SO$pca)[2] >=6){
    png(sub("$", "/post_dimHeatPCA.png", io$plotDir))
    DimHeatmap(object = SO
             , dims = 1:6, cells = 500, balanced = TRUE)
    dev.off()
}else{
    png(sub("$", "/post_dimHeatPCA.png", io$plotDir))
    DimHeatmap(object = SO
             , dims = 1:dim(SO$pca)[2], cells = 500, balanced = TRUE)
    dev.off()
}

# Same as above but for cell cycle

if(dim(CC_SO$pca)[2] >=6){
    png(sub("$", "/post_dimHeatPCA_CCreduced.png", io$plotDir))
    DimHeatmap(object = CC_SO
             , dims = 1:6, cells = 500, balanced = TRUE)
    dev.off()
}else{
    png(sub("$", "/post_dimHeatPCA_CCreduced.png", io$plotDir))
    DimHeatmap(object = CC_SO
             , dims = 1:dim(SO$pca)[2], cells = 500, balanced = TRUE)
    dev.off()
}

## may want to add an option for ndims
png(sub("$", "/post_elbowPlot.png", io$plotDir))
ElbowPlot(object = SO, ndims = 60)
dev.off()

SO        <- FindNeighbors(SO, reduction = "pca", dims = 1:40)
SO        <- FindClusters(object = SO, resolution = opts$res)
SO        <- RunUMAP(SO, reduction = "pca", dims = 1:40)

png(sub("$", "/post_SO_UMAP.png", io$plotDir))
DimPlot(SO, reduction = "umap", group.by = "seurat_clusters")
dev.off()

png(sub("$", "/CCPhase_UMAP.png", io$plotDir))
DimPlot(SO, reduction = "umap", group.by = "Phase")
dev.off()

SO <- FindVariableFeatures(SO)
var_genes <- VariableFeatures(SO)
write.table(var_genes, file = "data/seurat/post_var_genes.tsv", sep = "\t", row.names=F, col.names=F, quote=F)

png(sub("$", "/origin_UMAP.png", io$plotDir))
DimPlot(SO, reduction = "umap", group.by = "origin")
dev.off()

    
## Same as above but for cell cycle
png(sub("$", "/post_elbowPlot_CCreduced.png", io$plotDir))
ElbowPlot(object = CC_SO, ndims = 60)
dev.off()

CC_SO        <- FindNeighbors(CC_SO, reduction = "pca", dims = 1:40)
CC_SO        <- FindClusters(object = CC_SO, resolution = opts$res)
CC_SO        <- RunUMAP(CC_SO, reduction = "pca", dims = 1:40)

png(sub("$", "/post_SO_UMAP_CCreduced.png", io$plotDir))
DimPlot(CC_SO, reduction = "umap", group.by = "seurat_clusters")
dev.off()

png(sub("$", "/CCPhase_UMAP_CCreduced.png", io$plotDir))
DimPlot(CC_SO, reduction = "umap", group.by = "Phase")
dev.off()

CC_SO <- FindVariableFeatures(CC_SO)
var_genes <- VariableFeatures(CC_SO)
write.table(var_genes, file = "data/seurat/post_var_genes_CCreduced.tsv", sep = "\t", row.names=F, col.names=F, quote=F)

png(sub("$", "/origin_UMAP_CCreduced.png", io$plotDir))
DimPlot(CC_SO, reduction = "umap", group.by = "origin")
dev.off()

## Plot Variable Genes ##

top10 <- head(VariableFeatures(SO))
plot1 <- VariableFeaturePlot(SO)
plot2 <- LabelPoints(plot=plot1, points=top10, repel = TRUE)

png(sub("$", "/post_var_genes_scatter.png", io$plotDir))
plot2
dev.off()

## Plot Variable Genes for CC reduced ##

top10 <- head(VariableFeatures(CC_SO))
plot1 <- VariableFeaturePlot(CC_SO)
plot2 <- LabelPoints(plot=plot1, points=top10, repel = TRUE)

png(sub("$", "/post_var_genes_scatter_CCreduced.png", io$plotDir))
plot2
dev.off()

## Find DE genes

SO.markers <- FindAllMarkers(SO, logfc.threshold = 0.2)
top10      <- SO.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
write.table(SO.markers, file = "data/seurat/post_DEgenes.tsv", sep = "\t", row.names=F, quote=F)

png(sub("$", "/post_DEheatmap.png", io$plotDir))
DoHeatmap(SO, features = top10$gene) + NoLegend()
dev.off()

## Find DE genes for CC reduced ##

CC_SO.markers <- FindAllMarkers(CC_SO, logfc.threshold = 0.2)
top10      <- CC_SO.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
write.table(CC_SO.markers, file = "data/seurat/post_DEgenes_CCreduced.tsv", sep = "\t", row.names=F, quote=F)

png(sub("$", "/post_DEheatmap_CCreduced.png", io$plotDir))
DoHeatmap(CC_SO, features = top10$gene) + NoLegend()
dev.off()

## Save Seurat Object ##
saveRDS(SO, io$out.file)

## Save CC reduced Seurat Object ##
saveRDS(CC_SO, io$out.file.CCred)
