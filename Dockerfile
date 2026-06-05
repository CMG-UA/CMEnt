# syntax=docker/dockerfile:1.4
FROM ubuntu:24.04 AS builder
ARG SAMTOOLS_VERSION=1.23.1

# install samtools/htslib/bgzip/tabix, to get the latest version (1.22.1 instead of 1.13)
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    ca-certificates \
    curl \
    bzip2 \
    build-essential \
    libcurl4-openssl-dev \
    libncurses-dev \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    lsb-release \
    libssl-dev && \
    rm -rf /var/lib/apt/lists/*

RUN mkdir -p /tmp/samtools_install && \
    cd /tmp/samtools_install && \
    curl -fsSL "https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2" | tar jxf - && \
    cd "samtools-${SAMTOOLS_VERSION}" && \
    ./configure --prefix=/opt/samtools && \
    make -j"$(nproc)" && \
    make install && \
    cd htslib-* && \
    make -j"$(nproc)" && \
    make install && \
    rm -rf /tmp/samtools_install

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

COPY --link --from=builder /opt/samtools/bin/* /usr/local/bin/
COPY --link --from=builder /opt/samtools/lib/* /usr/local/lib/
RUN mkdir -p "${HOME}" && \
    printf '%s\n' 'options(CMEnt.auto_install_dep_if_missing = TRUE)' > "${HOME}/.Rprofile"

# Install system dependencies
RUN apt update && apt install -y r-cran-devtools r-cran-remotes r-cran-optparse r-cran-biocmanager
# Pre-install some annotation packages
RUN Rscript -e "BiocManager::install(c('IlluminaHumanMethylation27kanno.ilmn12.hg19', 'IlluminaHumanMethylation450kanno.ilmn12.hg19', 'IlluminaHumanMethylationEPICanno.ilm10b4.hg19', 'IlluminaHumanMethylationEPICv2anno.20a1.hg38'), ask = FALSE, update = FALSE)"
RUN Rscript -e "BiocManager::install(c('JASPAR2024', 'BSgenome.Hsapiens.UCSC.hg38'), ask = FALSE, update = FALSE)"
RUN Rscript -e "pak::pkg_install('R6', upgrade = TRUE)"
RUN Rscript -e "devtools::install()"
RUN apt clean && rm -rf /var/lib/apt/lists/*

ENTRYPOINT ["./run_cment.R"]
