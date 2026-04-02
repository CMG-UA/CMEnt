FROM rocker/r2u:24.04
LABEL maintainer="vaslem"
ENV DEBIAN_FRONTEND="noninteractive" \
    TZ="Europe/Brussels" \
    HOME=/home/root \
    LANG=en_US.UTF-8 \
    LC_ALL=en_US.UTF-8

COPY . /DMRsegal
WORKDIR /DMRsegal
SHELL ["/bin/bash", "-c"]

# Install system dependencies
RUN Rscript -e "devtools::install()"

ENTRYPOINT ["inst/bin/find_dmrs_from_seeds.R"]

