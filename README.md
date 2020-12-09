# scRNA_SMARTSEQ2_transcriptomeMapping

This pipeline is for the processing of paired-end single-cell RNA sequencing. The pipeline takes individual FASTQs for each cell as input and proceeds to trim, map, filter and count the reads. Then using this counts table it produces a Seurat object with the data. Standard outputs include counts tables, metadata files, QC plots, and the Seurat object.

Trimming
--------
Pipeline uses TrimGalore to trim the samples as well as run fastQC on the reads. Outputs gzipped trimmed FASTQs, a trimming report, and fastqc reports.

Fastqscreen
-----------
Pipeline runs FastQ screen on the trimmed reads and produces even more QC reports.

Mapping
-------
This pipeline uses HISAT2 to map the reads to a genome of choice. Location of HISAT tool, genome index and GTF file needs to be specified in the config file. Rule is currently set up to run HISAT2 over 12 threads, but that can be changed in the align.smk rule. HISAT2 outputs SAM files by default, so an extra step is added to convert the SAM files into BAM files.

Feature Count
-------------
Uses featureCounts to count the mapped reads. Location of GTF file needs to be specified in the config file.

Filter Counts
-------------
This is a custom R script to filter the counts. Options to filter out mitochondrial genes and filter by biotype need to be specified in the config file. Also outputs a sample metadata file connecting the sample name, output bam, and cell name. Also outputs a gene metadata file containing chromosome location information, ensembl ID, and gene name for all genes in dataset.

Seurat
------
This takes the filtered count table and creates a Seurat Object using it. Produces QC plots as well as filters the Seurat object based on parameters that need to be set in the config file. There is also an option to integrate different sample together in the Seurat object. It will then cluster the dataset and produce UMAPs as well as output differential expression heatmaps. It also produces a second Seurat object which has the effects of cell cycle removed and does the same clustering and differential expression analysis on this second object.

