FROM bioconductor/bioconductor_docker:RELEASE_3_21

ENV DEBIAN_FRONTEND=noninteractive
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
    nodejs \
    npm \
    pandoc \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /requirements

COPY .github/docker/cment-env/cran-requirements.txt cran-requirements.txt
COPY .github/docker/cment-env/bioc-requirements.txt bioc-requirements.txt
COPY .github/docker/cment-env/install-requirements.R install-requirements.R

RUN Rscript install-requirements.R
RUN Rscript -e "install.packages(c('covr', 'xml2', 'rcmdcheck', 'mockery', 'knitr', 'rmarkdown', 'rcmdcheck', 'devtools', 'magick'), repos = 'https://cloud.r-project.org')"
RUN Rscript -e "BiocManager::install(c('IlluminaHumanMethylation450kanno.ilmn12.hg19', 'IlluminaHumanMethylationEPICanno.ilm10b4.hg19', 'BSgenome.Hsapiens.UCSC.hg19'), ask = FALSE, update = FALSE)"