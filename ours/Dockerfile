FROM r-base:4.0.2

MAINTAINER "Carlos Fernandez-Lozano"

RUN install2.r \
    BiocManager \
    sva \
    pkgconfig \
    biospear \
&& rm -rf /tmp/downloaded_packages/ /tmp/*.rds

# Bioconductor
RUN R -e "BiocManager::install('PCAtools')"

COPY Pred1ctors.R /Pred1ctors.R
COPY models /models

CMD ["Rscript", "/Pred1ctors.R", "/data/GRCh37ERCC_refseq105_genes_count.csv", "/output/predictions.csv"]
