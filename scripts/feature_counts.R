library(data.table)
library(purrr)
library(Rsubread)

fwrite_tsv <- partial(fwrite, sep = "\t", na = "NA")


## in/out ##
io <- list()
io$anno_infile <- "/home/groups/CEDAR/anno/CellRanger/hg19/refdata-cellranger-hg19-3.0.0/genes/genes.gtf"
#io$anno_infile <- "Box/My_NMT-seq_work/Stephen/code/MCF7_matched-selected/features/genes/Homo_sapiens.GRCh38.97.gtf.gz"

io$samples_indir <- "samples/hisat2"
io$outdir <- "data/counts"

## Options ##
opts <- list()
opts$nthreads <- 2
opts$allowMultiOverlap = FALSE        # should a read be allowed to be assigned to more than one feature?
opts$isPairedEnd = TRUE              # paired end data? If TRUE, fragments will be counted instead of individual reads
opts$requireBothEndsMapped = FALSE    # should both ends from the same fragment be required to be aligned?
opts$strandSpecific = 0               # 0=unstranded, 1=stranded, 2=reversely stranded
opts$countChimericFragments = FALSE   # should we allow a fragment to have its two reads mapped to different chromosomes?
opts$countMultiMappingReads <- FALSE  # logical indicating if multi-mapping reads/fragments should be counted
opts$isGTFAnnotationFile = TRUE

opts$generate_metadata <- TRUE




#######################
## Run featureCounts ##
#######################

files <- dir(io$samples_indir, full = T, pattern = ".bam$") %>%
  .[file.exists(.)]

fc <- featureCounts(files = files, 
                    annot.inbuilt = NULL, 
                    annot.ext = io$anno_infile, 
                    isGTFAnnotationFile = opts$isGTFAnnotationFile, 
                    useMetaFeatures = TRUE, 
                    allowMultiOverlap = opts$allowMultiOverlap, 
                    isPairedEnd = opts$isPairedEnd, 
                    requireBothEndsMapped = opts$requireBothEndsMapped, 
                    nthreads = opts$nthreads, 
                    strandSpecific = opts$strandSpecific, 
                    countChimericFragments = opts$countChimericFragments)

counts <- as.data.frame(fc$counts) %>% 
  setDT(keep.rownames = "ens_id")

old_names <- colnames(counts)[!colnames(counts) %in% "ens_id"]
new_names <- basename(old_names) %>%
  tools::file_path_sans_ext() # %>%
#  gsub("M7C2.|.output", "", .) %>%
#  gsub("_L001.*", "", .)

setnames(counts, old_names, new_names)
setcolorder(counts, gtools::mixedsort(colnames(counts)))


##################
## Save results ##
##################

dir.create(io$outdir, recursive = TRUE)
#now <- format(Sys.time(), "%Y-%m-%d_%Hhr%M" )
outfile <- paste0(io$outdir, "/raw_counts_", ".tsv") 


fwrite_tsv(counts, outfile)

if (opts$generate_metadata) {
  meta <- data.table(rna_file = files, id_rna = new_names)
  fwrite_tsv(meta, paste0(io$outdir, "/sample_metadata.tsv"))
}
