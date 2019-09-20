library(data.table)
library(purrr)
library(rtracklayer)
## Download gtf annotation from Ensembl ##

## use gtf that we already have downloaded in 
#gtf_url <- "ftp://ftp.ensembl.org/pub/grch37/release-87/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.chr.gtf.gz"
#(outfile <- paste0("/home/groups/CEDAR/woodfin/projects/NMT-seq/20190523_NM/anno/hg19/", basename(gtf_url)))
#dir.create(dirname(outfile), recursive = TRUE)

#download.file(gtf_url, outfile)

#outfile <- "/home/groups/CEDAR/anno/gtf/hg19_ens87.chr.gtf"

## also generate gene metadata ##
#chrs <- sub("^", "chr", c("X", "Y", "MT", 1:22))

#gtf <- rtracklayer::import(outfile) %>%
#  as.data.frame() %>%
#  setDT() %>%
#  .[type == "gene" & gene_biotype == "protein_coding" & seqnames %in% chrs, 
#    .(ens_id = gene_id,
#      gene = gene_name, 
#      chr = seqnames, 
#      start, 
#      end, 
#      strand)] 
  
#gene_file <- "data/gene_metadata.tsv"

#fwrite(gtf, gene_file, sep = "\t", na = "NA")

#gtf$score <- 1000
#write.table(gtf[,c("chr", "start", "end", "gene", "score", "strand")]
#          , sub(".tsv", ".bed", gene_file), sep = "\t", quote=FALSE, row.names=FALSE, col.names=FALSE)


############
#gtf_url <- "ftp://ftp.ensembl.org/pub/grch37/release-97/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz"
#(outfile <- paste0("/home/groups/CEDAR/woodfin/projects/NMT-seq/20190523_NM/anno/hg19/", basename(gtf_url)))
#dir.create(dirname(outfile), recursive = TRUE)

#download.file(gtf_url, outfile)


## also generate gene metadata ##
#chrs <- sub("^", "chr", c("X", "Y", "MT", 1:22))

#chrs <- c("X", "Y", "MT", 1:22)

#gtf <- rtracklayer::import(outfile)


#gtf <- as.data.frame(gtf) %>%
#  setDT() %>%
#  .[type == "gene" & gene_biotype == "protein_coding" & seqnames %in% chrs, 
#    .(ens_id = gene_id,
#      gene = gene_name, 
#      chr = seqnames, 
#      start, 
#      end, 
#      strand)] 
  
#gene_file <- paste0(dirname(outfile), "/gene_hg19.87_metadata.tsv")

#fwrite(gtf, gene_file, sep = "\t", na = "NA")


############
outfile <- "/home/groups/CEDAR/anno/CellRanger/hg19/refdata-cellranger-hg19-3.0.0/genes/genes.gtf"


## also generate gene metadata ##
chrs <- sub("^", "chr", c("X", "Y", "MT", 1:22))

chrs <- c("X", "Y", "MT", 1:22)

gtf <- rtracklayer::import(outfile)


gtf <- as.data.frame(gtf) %>%
  setDT() %>%
  .[type == "gene" & gene_biotype == "protein_coding" & seqnames %in% chrs, 
    .(ens_id = gene_id,
      gene = gene_name, 
      chr = seqnames, 
      start, 
      end, 
      strand)] 
  
write.table(gtf, file="data/gene_hg19.cellRanger_metadata.tsv", sep = "\t")

#library(GenomicFeatures)

#gtf <- gtf[!gtf$chr %in% c("MT", "Y")]

#tss <- promoters(as(gtf, "GRanges"), upstream=2000, downstream=2000)

#tss <- as.data.table(tss)
#names(tss) <- sub("seqnames", "chr", names(tss))
