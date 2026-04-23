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
RUN Rscript -e "install.packages(c('devtools', 'remotes'))"
RUN Rscript -e "devtools::install()"
RUN Rscript -e "install.packages(c('optparse', 'BiocManager'))"
# Pre-install array annotation packages used when opening saved viewer outputs.
RUN Rscript -e "BiocManager::install(c('IlluminaHumanMethylation27kanno.ilmn12.hg19', 'IlluminaHumanMethylation450kanno.ilmn12.hg19', 'IlluminaHumanMethylationEPICanno.ilm10b4.hg19', 'IlluminaHumanMethylationEPICv2anno.20a1.hg38'), ask = FALSE, update = FALSE)"

ENTRYPOINT ["inst/bin/run_cment.R"]
