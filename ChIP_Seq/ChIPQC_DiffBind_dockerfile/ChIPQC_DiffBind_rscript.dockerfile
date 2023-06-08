FROM bioconductor/bioconductor_docker:latest
# create workdir
RUN mkdir -p /home/analysis/custom_R_installation_dir
# install ballgown library from Bioconductor
RUN R --vanilla -e 'BiocManager::install(c("DESeq2","DiffBind","ChIPQC","rtracklayer","TxDb.Hsapiens.UCSC.hg38.knownGene"), lib="/home/analysis/custom_R_installation_dir")'
# install data.table and openxlsx libraries
RUN install2.r \
    --error \
    --libloc /home/analysis/custom_R_installation_dir \
    data.table \
    stringr
# switch to workdir
WORKDIR /home/analysis
# copy R scripts
COPY ChIPQC.R .
COPY DiffBind.R .
