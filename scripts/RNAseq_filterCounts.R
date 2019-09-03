annoFile = snakemake@params[['anno']]

biotypes <- snakemake@params[['biotypes']]

countsFile <- snakemake@input[['countsFile']]

mito <- snakemake@params[['mito']]

counts <- read.table(file=countsFile, sep = "\t", header = T, skip = 1)

##----------load counts------------#
print("Loading counts table")
print(countsFile)

## must be a tsv or txt tab sep file

counts <- read.table(file=countsFile, sep = "\t", header = T)

##----------load anno------------#
print("Loading annotation table")
print(annoFile)

## load
anno             <- get(load(file=annoFile))

## check if all ensembl ids in the counts are in the anno file
test_list        <- anno$external_gene_name
names(test_list) <- anno$ensembl_gene_id
counts$Genes     <- "Not_in_anno"

if( sum(! paste(counts$Geneid) %in% paste(anno$ensembl_gene_id) ) > 0 ){
    print( "not all ensembl ids in counts in annoFile" )
    for (i in seq(1:length(counts$Genes))) {
        if (counts$Geneid[i] %in% names(test_list)) {
            counts$Genes[i] <- test_list[ paste(counts$Geneid[i]) ]
        }
    }
}else{
    iv           <- match( paste(counts$Geneid) , paste(anno$ensembl_gene_id) )
    counts$Genes <- anno[iv,"external_gene_name"]
}

if(strsplit(biotypes, split='\\,')[[1]]!=""){
    anno.sub   <- anno[paste(anno$gene_biotype) %in% strsplit(biotypes, split='\\,')[[1]] ,]
    counts.sub <- counts[paste(counts$Genes) %in% unique(paste(anno.sub$external_gene_name)) , ]
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

if( length(grep("samples.hisat2.|_output.sam",names(counts.sub)) ) > 0 ){
    names(counts.sub) <- gsub("samples.hisat2.|_output.sam", "", names(counts.sub) )
}

# take chrom name
counts.sub$Chr    <- sub("\\;.*", "", paste(counts.sub$Chr))

# take the first Etart
counts.sub$Start  <- sub("\\;.*", "", paste(counts.sub$Start))

# take the last End
counts.sub$End    <- sub(".*\\;", "", paste(counts.sub$End))

# take the first strand
counts.sub$Strand <- sub("\\;.*", "", paste(counts.sub$Strand))


write.table(counts.sub, file=sub(".txt", ".filt.txt", countsFile), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
