# syntax=docker/dockerfile:1.4
ARG BASE_IMAGE
FROM ${BASE_IMAGE} AS final

FROM continuumio/miniconda3:latest AS builder
RUN conda install -y -c conda-forge -c bioconda combined-pvalues bedtools && \
    conda clean -afy

FROM final

COPY .github/docker/cment-env-notebooks/cran-requirements.txt cran-requirements.txt
COPY .github/docker/cment-env-notebooks/bioc-requirements.txt bioc-requirements.txt
COPY .github/docker/cment-env-notebooks/install-requirements.R install-requirements.R
COPY --from=builder /opt/conda /opt/conda
ENV PATH="/opt/conda/bin:${PATH}"
RUN Rscript install-requirements.R && \
    rm -rf /tmp/downloaded_packages
