## must recreate all envs
## you will want to remove all packages and envs after export to prevent your home directory from hitting 10G and locking you out

#################
## rules/align_rmdp.smk 
#################
conda env create -f envs/trimG.yaml
conda activate trimG.yaml
conda env export --no-builds > envs/trimG.yaml
conda deactivate
conda remove --name trimG.yaml --all

## not working make from scratch
conda create -f envs/Fastqscreen.yaml --name screen

conda create --name Fastqscreen -c bioconda fastq-screen
conda activate Screen
conda env export --no-builds > envs/fastqscreen.yaml
conda deactivate
conda remove Screen

conda env create -f featureCounts.yaml
conda activate Fcounts
conda env export --no-builds > envs/featureCounts.yaml
conda deactivate 
conda remove --name Fcounts --all

#############
## rules/analysis.smk
#############
conda create -n smartseq2 r-seurat bioconductor-singlecellexperiment r-data.table bioconductor-rsubread r-purrr bioconductor-rtracklayer r-ggplot2 r-cowplot bioconductor-scater r-png
conda activte smartseq2
conda env export --no-builds > envs/scAnalysis.yaml
conda deactivate
conda remove --name smartseq2 --all


conda create env envs/NMT_Create_Seurat.yaml
## failed
conda activate NMT_Create_Seurat
conda env export --no-builds > envs/NMT_Create_Seurat.new.yaml
conda deactivate
conda remove --name NMT_create_Seurat --all

grep prefix envs/* | sed "s/\:.*//g" | while read line; do
    sed -i "s/prefix.*//g" $line
done

## make a new hisat2 samtools env for mapping
conda create --prefix=/home/groups/CEDAR/woodfin/projects/CHC/wong_mohammad/20200309_CHC_PDAC_scRNA_myEnvs/scNMT_transcriptomeMapping/HiSat2
conda activate /home/groups/CEDAR/woodfin/projects/CHC/wong_mohammad/20200309_CHC_PDAC_scRNA_myEnvs/scNMT_transcriptomeMapping/HiSat2
conda install -c bioconda hisat2
conda install -c bioconda samtools
conda env export --no-builds > envs/Mapping.yaml
conda deactivate
sed -i "s/prefix.*//g" envs/Mapping.yaml
