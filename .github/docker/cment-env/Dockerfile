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
ENV BIOCONDUCTOR_VERSION=3.23 \
    DEBIAN_FRONTEND="noninteractive" \
    TZ="Europe/Brussels" \
    HOME=/home/root \
    LANG=en_US.UTF-8 \
    LC_ALL=en_US.UTF-8 \
    PATH="/usr/lib/google-cloud-sdk/bin:${PATH}"
WORKDIR /home/root
# Setting bash as the default shell
SHELL ["/bin/bash", "-c"]

ENV R_COMPILE_AND_INSTALL_PACKAGES=never

RUN apt-get update && apt-get install -y --no-install-recommends \
    curl \
    git \
    libcurl4-openssl-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libfribidi-dev \
    libharfbuzz-dev \
    libjpeg-dev \
    libpng-dev \
    libssl-dev \
    libtiff5-dev \
    libxml2-dev \
    pandoc \
    ghostscript \
    texlive-latex-base \
    texlive-fonts-recommended \
    texlive-fonts-extra \
    texlive-latex-extra

# Install tabix (through samtools) from the builder stage
COPY --from=builder /opt/samtools /opt/samtools
COPY --link --from=builder /opt/samtools/bin/* /usr/local/bin/
COPY --link --from=builder /opt/samtools/lib/* /usr/local/lib/

COPY .github/docker/cment-env/docker.Rprofile .Rprofile
ENV R_PROFILE_USER=/home/root/.Rprofile
RUN mkdir -p "$HOME/.cache/R/ExperimentHub" && \
    for i in {1..5}; do \
        Rscript -e "install.packages('BiocManager'); BiocManager::install(version='${BIOCONDUCTOR_VERSION}', ask=FALSE)" && \
        break || sleep 30; \
    done && \
    rm -rf /tmp/downloaded_packages

COPY .github/docker/cment-env/cran-requirements.txt cran-requirements.txt
COPY .github/docker/cment-env/bioc-requirements.txt bioc-requirements.txt
COPY .github/docker/cment-env/install-requirements.R install-requirements.R

RUN Rscript install-requirements.R
# Install additional packages that are not in the requirements files but are needed for testing and development
RUN Rscript -e "install_deps(packages=c('covr', 'xml2', 'rcmdcheck', 'mockery', 'knitr', 'rmarkdown', 'rcmdcheck', 'devtools', 'magick'), bioc_packages=c('BiocCheck', 'IlluminaHumanMethylation450kanno.ilmn12.hg19', 'IlluminaHumanMethylationEPICanno.ilm10b4.hg19', 'BSgenome.Hsapiens.UCSC.hg19'))"

RUN apt-get clean \
    && rm -rf /var/lib/apt/lists/*