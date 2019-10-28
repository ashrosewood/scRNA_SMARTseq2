annoFile = snakemake@params[['anno']]
#annoFile <- "/home/groups/CEDAR/anno/biomaRt/hg19.Ens_87.biomaRt.geneAnno.Rdata"

biotypes <- snakemake@params[['biotypes']]
#biotypes <- "protein_coding"

countsFile <- snakemake@input[['countsFile']]
#countsFile <- "../data/counts/raw_counts_.tsv"

mito <- snakemake@params[['mito']]
#mito <- 0

##----------load counts------------#
print("Loading counts table")
print(countsFile)

## must be a tsv or txt tab sep file

counts           <- read.table(file=countsFile, sep = "\t", header = T, skip=1)
names(counts)    <- sub("Geneid", "ens_id", names(counts))
rownames(counts) <- counts$ens_id

if( length(grep("samples.hisat2.|.bam",names(counts))) > 0 ){
    names(counts) <- gsub("samples.hisat2.|_output.bam",'',names(counts))
}

# take chrom name
counts$Chr    <- sub("\\;.*", "", paste(counts$Chr))

# take the first Etart
counts$Start  <- sub("\\;.*", "", paste(counts$Start))

# take the last End
counts$End    <- sub(".*\\;", "", paste(counts$End))

# take the first strand
counts$Strand <- sub("\\;.*", "", paste(counts$Strand))


##----------load anno------------#
print("Loading annotation table")
print(annoFile)

## load
anno             <- get(load(file=annoFile))

## check if all ensembl ids in the counts are in the anno file
test_list        <- anno$external_gene_name
names(test_list) <- anno$ensembl_gene_id
counts$Genes     <- "Not_in_anno"

if( sum(! paste(counts$ens_id) %in% paste(anno$ensembl_gene_id) ) > 0 ){
    print( "not all ensembl ids in counts in annoFile" )
    for (i in seq(1:length(counts$Genes))) {
        if (counts$ens_id[i] %in% names(test_list)) {
            counts$Genes[i] <- test_list[ paste(counts$ens_id[i]) ]
        }
    }
}else{
    iv           <- match( paste(counts$ens_id) , paste(anno$ensembl_gene_id) )
    counts$Genes <- anno[iv,"external_gene_name"]
}

if(strsplit(biotypes, split='\\,')[[1]]!=""){   
    anno.sub   <- anno[paste(anno$gene_biotype) %in% strsplit(biotypes, split='\\,')[[1]] ,]
    counts.sub <- counts[paste(counts$ens_id) %in% unique(paste(anno.sub$ensembl_gene_id)) , ]
}else{
    print("no biotypes provided")
    counts.sub <- counts
}

## Note this will only work for orgnaisms that's mito gene start with "MT-", mouse is "mt"
if(mito==1){
    print("tossing MT- genes")
    counts.sub <- counts.sub[grep("^MT-", paste(counts.sub$Genes), invert=TRUE), ]
}

#---------------clean up the counts names------------------#
# take chrom name
counts.sub$Chr    <- sub("\\;.*", "", paste(counts.sub$Chr))

# take the first Etart
counts.sub$Start  <- sub("\\;.*", "", paste(counts.sub$Start))

# take the last End
counts.sub$End    <- sub(".*\\;", "", paste(counts.sub$End))

# take the first strand
counts.sub$Strand <- sub("\\;.*", "", paste(counts.sub$Strand))

Anno <- counts.sub[, c("ens_id", "Genes", "Chr", "Start", "End", "Strand")]
names(Anno) <- c("ens_id", "gene", "chr", "start", "end", "strand")

write.table(Anno, "data/gene_metadata.tsv", sep = "\t", na = "NA", row.names= FALSE, col.names = TRUE)

counts.sub$Chr    <- NULL
counts.sub$Start  <- NULL
counts.sub$End    <- NULL
counts.sub$Strand <- NULL
counts.sub$Length <- NULL

write.table(counts.sub, file=sub(".tsv", ".filt.tsv", countsFile), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

print(paste("getwd:", getwd()))

                                        #files <- list.files(sub("$", "/samples/hisat2", dirname(getwd())), full.names=T)
files <- list.files("samples/hisat2", full.names=T)
files <- sub("^", paste0(getwd(), "/"), files)
print(head(files))

sample_meta <- data.frame(rna_file=files
                          ,id_rna=sub("_output.bam", "", basename(files)))
print(head(sample_meta))

write.table(sample_meta, "data/counts/sample_metadata.tsv", sep = "\t", na = "NA", row.names= FALSE, col.names = TRUE)
