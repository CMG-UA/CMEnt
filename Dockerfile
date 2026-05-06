FROM rocker/r2u:24.04
LABEL maintainer="vaslem"
ENV DEBIAN_FRONTEND="noninteractive" \
    TZ="Europe/Brussels" \
    HOME=/home/root \
    LANG=en_US.UTF-8 \
    LC_ALL=en_US.UTF-8

COPY . /CMEnt
WORKDIR /CMEnt
SHELL ["/bin/bash", "-c"]

# Install system dependencies
RUN apt update && apt install -y r-cran-devtools r-cran-remotes r-cran-optparse r-cran-biocmanager
RUN Rscript -e "pak::pkg_install('R6', upgrade = TRUE)"
RUN Rscript -e "devtools::install()"
# Pre-install array annotation packages used when opening saved viewer outputs.
RUN Rscript -e "BiocManager::install(c('IlluminaHumanMethylation27kanno.ilmn12.hg19', 'IlluminaHumanMethylation450kanno.ilmn12.hg19', 'IlluminaHumanMethylationEPICanno.ilm10b4.hg19', 'IlluminaHumanMethylationEPICv2anno.20a1.hg38'), ask = FALSE, update = FALSE)"
RUN apt clean && rm -rf /var/lib/apt/lists/*

ENTRYPOINT ["inst/bin/run_cment.R"]
