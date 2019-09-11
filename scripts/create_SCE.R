library(SingleCellExperiment)
library(data.table)
library(purrr)
library(scater)
library(scran)
library(ggplot2)
library(cowplot)

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
opts$coverage_threshold <- 2e5    # Minimum library size (coverage)
opts$features_threshold <- 1000   # Minimum number of expressed features
opts$top50_threshold <- 0.75      # Maximum fraction of reads accounting for the top 50 features

# Lenient thresholds
#opts$coverage_threshold <- 1e4    # Minimum library size (coverage)
#opts$features_threshold <- 500   # Minimum number of expressed features
#opts$top50_threshold <- 0.75      # Maximum fraction of reads accounting for the top 50 features

opts$MT_threshold <- 0.25         # Maximum fraction of reads mapping to mithocondrial genes

## I/O ##
io <- list()
io$in.gene_metadata <- "Box/My_NMT-seq_work/Stephen/code/MCF7_matched-selected/features/genes_hg19/gene_metadata.tsv.gz"
io$in.sample_metadata <- "Box/My_NMT-seq_work/Stephen/code/MCF7_matched-selected/rna/counts_hg19/sample_metadata.tsv"
io$in.raw_counts <- "Box/My_NMT-seq_work/Stephen/code/MCF7_matched-selected/rna/counts_hg19/raw_counts_2019-08-25_23hr20.tsv"
io$out.file <- "Box/My_NMT-seq_work/Stephen/code/MCF7_matched-selected/rna/parsed_hg19/SingleCellExperiment.rds"
io$out.sample_metadata <- "Box/My_NMT-seq_work/Stephen/code/MCF7_matched-selected/sample_metadata.tsv"



## read in counts data and associated metadata ##

counts <- fread_df(io$in.raw_counts) %>% 
  tibble::column_to_rownames("ens_id") %>% 
  as.matrix()

sample_metadata <- fread_df(io$in.sample_metadata) 
rownames(sample_metadata) <- sample_metadata$id_rna

cells <- sample_metadata$id_rna %>% 
  .[. %in% colnames(counts)]

sample_metadata <- sample_metadata[cells, ]
counts <- counts[, cells]

feature_metadata <- fread_df(io$in.gene_metadata) 

# Define mithocondrial genes
mt <- feature_metadata$gene[feature_metadata$chr == "MT"]

# rename duplicated genes 
feature_metadata$gene <- make.unique(feature_metadata$gene)

rownames(feature_metadata) <- feature_metadata$ens_id


## Parse data ##

genes <- rownames(feature_metadata[rownames(feature_metadata) %in% rownames(counts),])
feature_metadata <- feature_metadata[genes,]
counts <- counts[rownames(feature_metadata),]


rownames(counts) <- NULL
rownames(feature_metadata) <- NULL
## Create SCEset object ##
  
sce <- SingleCellExperiment(assays = list(counts = counts), 
                            rowData = feature_metadata, 
                            colData = sample_metadata)

rownames(sce) <- genes

# Calculate quality metrics
sce <- calculateQCMetrics(sce, feature_controls=list(Mt=rownames(sce) %in% mt))

sce$gene_strand <- sce$strand
sce$strand <- NULL



# Filter by library size
libsize.drop <- sce$total_counts < opts$coverage_threshold
sum(libsize.drop)

libsize.drop_dt <- data.table(
  sample=colnames(sce), 
  size=sce$total_counts, 
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
save_plot("Box/My_NMT-seq_work/Stephen/code/MCF7_matched-selected/plots/rna_libsize_hg19.pdf", p1)



# Filter by number of expressed genes
feature.drop <- sce$total_features_by_counts < opts$features_threshold
feature.drop_dt <- data.table(
  sample = colnames(sce), 
  features = sce$total_features_by_counts, 
  color = c("black","red")[as.numeric(feature.drop)+1]
) %>% setkey(features) %>% .[,col:=features] %>% .[,sample:=factor(sample,levels=sample)]

p2 <- ggplot(feature.drop_dt, aes(x=sample, y=features)) +
  geom_bar(stat='identity', position="dodge", fill="#3CB54E") +
  geom_hline(yintercept=opts$features_threshold, colour="black", linetype="dashed") +
  # scale_fill_gradient(low="red", high="green") +
  labs(y="Total number of expressed genes") +
  #barplot_theme() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=rel(1.8)),
    # axis.text.x = element_text(colour="black", color=foo$color, angle=90, size=10, vjust=0.5, hjust=1.0)
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )
print(p2)
save_plot("Box/My_NMT-seq_work/Stephen/code/MCF7_matched-selected/plots/rna_features_hg19.pdf", p2)

# Samples with large proportion of mithocondrial genes
mt.drop <- sce$pct_counts_Mt > opts$MT_threshold*100

# Remove cells that do not pass QC

#drop.samples <- colnames(sce)[( libsize.drop | feature.drop | top50.drop | mt.drop)]

drop.samples <- colnames(sce)[( libsize.drop | feature.drop |  mt.drop)]
length(drop.samples)
#[1] 16

# Update sample metadata
setDT(sample_metadata)
sample_metadata[,pass_rnaQC:=ifelse(id_rna%in%drop.samples,FALSE,TRUE)]

# Save updated metadata
fwrite(sample_metadata, io$out.sample_metadata, sep="\t", na="NA")

# Re-calculate QC statistics
sce <- sce[,!colnames(sce) %in% drop.samples]
sce <- calculateQCMetrics(sce)


# # Normalisation and log transformation ##
# Transcript counts are now normalised based on size factors using the convolution approach from the scran package.
# Lowly expressed genes are removed before normalisation but they are included afterwards, since they are interesting for some analysis.




sce_filt <- sce

# Temporarily remove the lowly expressed genes
rowData(sce_filt)$gene_strand <- rowData(sce_filt)$strand
rowData(sce_filt)$start_bp <- rowData(sce_filt)$start
rowData(sce_filt)$end_bp <- rowData(sce_filt)$end

rowData(sce_filt)$start <- NULL
rowData(sce_filt)$end <- NULL
rowData(sce_filt)$strand <- NULL
sce_filt$start <- NULL
sce_filt$end <- NULL
sce_filt$strand <- NULL
# rowData(sce)$strand <- NULL
# rowData(sce)$gene_strand <- NULL
sce_filt <- sce_filt[!(rowMeans(counts(sce)) <= 1 | rowData(sce)$pct_dropout_by_counts > 90),]

# Compute size factors without the lowly expressed genes

sf = computeSumFactors(sce_filt, positive=TRUE, sf.out=T, sizes=seq(21, 77, 5))

# qplot(sf, sce_filt$total_counts, log="xy", ylab="Library size (mapped reads)", xlab="Size factor")
p3 <- ggplot(data.frame(sf=log(sf), counts=log(sce_filt$total_counts))) +
  geom_point(aes(x=sf,y=counts)) +
  labs(y="Library size (log)", x="Size factor (log)") +
  theme_bw() +
  theme(
    axis.title = element_text(colour="black", size=15),
    axis.text = element_text(colour="black", size=12)
  )

save_plot("Box/My_NMT-seq_work/Stephen/code/MCF7_matched-selected/plots/rna_size_factors.pdf", p3)

# Normalise and log transform with the lowly expressed genes
sizeFactors(sce) <- sf; sce$sizeFactor <- sf
sizeFactors(sce_filt) <- sf; sce_filt$sizeFactor <- sf
sce <- normalize(sce, exprs_values="counts")
sce_filt <- normalize(sce_filt, exprs_values="counts")

# Update quality metrics
sce = calculateQCMetrics(sce)


foo <- data.frame(sd=apply(exprs(sce),1,sd), mean=apply(exprs(sce),1,mean))
p4 <- ggplot(foo, aes(x=mean, y=sd)) +
  geom_point() +
  stat_smooth() +
  scale_color_manual(values=c("black","red")) +
  xlab('Mean') + ylab('Standard deviation')
print(p4)

save_plot("Box/My_NMT-seq_work/Stephen/code/MCF7_matched-selected/plots/rna_mean_var.pdf", p4)
## Save SCESet object ##

dir.create(dirname(io$out.file), recursive = TRUE)
saveRDS(sce,io$out.file)
