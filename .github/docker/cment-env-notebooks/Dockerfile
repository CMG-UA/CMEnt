# syntax=docker/dockerfile:1.4
ARG BASE_IMAGE
FROM ${BASE_IMAGE} AS final

FROM continuumio/miniconda3:latest AS builder
ARG COMBP_SHA=c18e137ef42e32d9ec6acdde7562d2338db0be14
RUN conda install -y -c conda-forge -c bioconda bedtools && \
    python -m pip install --no-cache-dir "https://github.com/brentp/combined-pvalues/archive/${COMBP_SHA}.tar.gz" && \
    conda clean -afy

FROM final

COPY .github/docker/cment-env-notebooks/cran-requirements.txt cran-requirements.txt
COPY .github/docker/cment-env-notebooks/bioc-requirements.txt bioc-requirements.txt
COPY .github/docker/cment-env-notebooks/install-requirements.R install-requirements.R
COPY --from=builder /opt/conda /opt/conda
ENV PATH="/opt/conda/bin:${PATH}"
RUN Rscript install-requirements.R && \
    rm -rf /tmp/downloaded_packages
