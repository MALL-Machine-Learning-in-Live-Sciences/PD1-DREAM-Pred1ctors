FROM r-base:4.0.2

MAINTAINER "Carlos Fernandez-Lozano"

RUN apt-get update && apt-get install -y --no-install-recommends \
	libgdal-dev \
	curl \
&& rm -rf /var/lib/apt/lists/*

RUN install2.r \
    BiocManager \
    sva \
    pkgconfig \
    biospear \
&& rm -rf /tmp/downloaded_packages/ /tmp/*.rds

# Bioconductor
RUN R -e "BiocManager::install('PCAtools')"

RUN R -e "BiocManager::install('sva')"

RUN R -e "BiocManager::install('mixOmics')"
RUN R -e "BiocManager::install('survcomp')"

RUN install2.r \
    biospear \
&& rm -rf /tmp/downloaded_packages/ /tmp/*.rds

RUN install2.r \
    data.table \
&& rm -rf /tmp/downloaded_packages/ /tmp/*.rds

RUN R -e "BiocManager::install('impute')"

COPY Pred1ctors_os_v3.R /Pred1ctors_os_v3.R
COPY models /models
COPY utils.r /utils.r

CMD ["Rscript", "/Pred1ctors_os_v3.R", "/data/GRCh37ERCC_refseq105_genes_count.csv","/data/clinical_data.csv" ,"/output/predictions.csv"]
