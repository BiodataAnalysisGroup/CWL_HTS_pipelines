FROM bioconductor/bioconductor_docker:latest
# create workdir
RUN mkdir -p /home/analysis/custom_R_installation_dir
# install ballgown library from Bioconductor
RUN R --vanilla -e 'BiocManager::install(c("DESeq2","apeglm","BiocParallel","sva","RUVSeq"), lib="/home/analysis/custom_R_installation_dir")'
# install data.table and openxlsx libraries
RUN install2.r \
    --error \
    --libloc /home/analysis/custom_R_installation_dir \
    data.table \
    dplyr
# switch to workdir
WORKDIR /home/analysis
# copy rscript
COPY DESeq2.R .
