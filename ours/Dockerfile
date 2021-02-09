FROM r-base:4.0.2

MAINTAINER "Carlos Fernandez-Lozano"

RUN apt-get update && apt-get install -y --no-install-recommends \
        libgdal-dev \
        curl \
&& rm -rf /var/lib/apt/lists/*

HEALTHCHECK --interval=10s --timeout=10s --start-period=60s --retries=3 \
        CMD curl -k --fail http://localhost:3838

# Packages
RUN install2.r --error \
    data.table \
    BiocManager \
    minqa \
    nloptr \
    RcppEigen \
    profileModel \
    pbkrtest \
    lme4 \
    brglm \
    gtools \
    car \
    BradleyTerry2 \
    doMC \
    caret \
&& rm -rf /tmp/downloaded_packages/ /tmp/*.rds

RUN install2.r \
    caretEnsemble \
    CHAID \
    class \
    Cubist \
    dplyr \
    e1071 \
    foreach \
    frbs \
    gam \
    import \
    inTrees \
    ipred \
    kernlab \
    MASS \
    Matrix \
    mboost \
    mda \
    mgcv \
    party \
    plyr \
    protoclass \
    proxy \
    randomForest \
    ranger \
    rocc \
    rpart \
    rpartScore \
&& rm -rf /tmp/downloaded_packages/ /tmp/*.rds

# Bioconductor
RUN R -e "BiocManager::install('NOISeq')"
RUN R -e "BiocManager::install('PCAtools')"

RUN install2.r \
    sva \
    pkgconfig \
    biospear \
&& rm -rf /tmp/downloaded_packages/ /tmp/*.rds

COPY Pred1ctors.R /Pred1ctors.R

CMD ["Rscript", "/Pred1ctors.R", "/data/GRCh37ERCC_refseq105_genes_count.csv", "/output/predictions.csv"]
