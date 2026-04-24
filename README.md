<!-- README.md is generated from README.Rmd. Please edit that file -->

# CMEnt <img src="man/figures/logo.png" align="right" height="139"/>

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![R-CMD-check](https://github.com/CMG-UA/CMEnt/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/CMG-UA/CMEnt/actions/workflows/R-CMD-check.yaml)
[![codecov](https://codecov.io/gh/CMG-UA/CMEnt/graph/badge.svg?token=58D9F1XJHY)](https://codecov.io/gh/CMG-UA/CMEnt)
<!-- badges: end -->

## Overview

CMEnt is an R package for identifying Differentially Methylated
Regions (DMRs) from pre-computed seeds, commonly identified
Differentially Methylated Positions (DMPs). The package implements a
correlation-based approach to expand regions around significant DMPs,
considering both statistical significance and biological relevance of
methylation changes.

## Features

- Correlation-based region expansion from genomic seeds
- DMRs scoring based on 5-fold SVM predictions
- Motif discovery within DMRs
- DMRs motif-based interactions
- lookup of DMRs motifs matching binding factors on JASPAR database
- Support for 450K, EPIC, EPICv2,
  [mouse](https://github.com/chiaraherzog/IlluminaMouseMethylationanno.12.v1.mm10)
  array platforms **and** NGS datasets in the form of bed (or tabix)
  files
- Support for genome builds hg19, hg38, hs1, mm10 and mm39
- Visualization functions for DMRs and methylation profiles

## Considerations

- CMEnt is designed to work with pre-computed seeds (DMPs) and does
  not perform DMP identification itself. Users should first identify
  DMPs using their preferred method before using CMEnt for region
  expansion.
- CMEnt is not preprocessing input methylation files (e.g. no
  normalization, no filtering, no probe removal), apart from conversion
  to different formats for efficiency, and optional confounders
  adjustment. It is recommended to preprocess beta files using
  established pipelines (e.g. minfi, meffil) before using CMEnt for
  DMR identification.

## Installation

You can install the development version of CMEnt from GitHub.
Currently the repository is private, so you need to use the organization
read-only access authentication key:

    # install.packages("devtools") # nolint
    devtools::install_github("CMG-UA/CMEnt")

Another option is to use the Docker image available on DockerHub, which
contains a pre-installed version of CMEnt along with all
dependencies. You can pull the image using the following command:

    docker pull vlemonidis/cment:latest

and run it using:

    docker run --rm vlemonidis/cment:latest --help

To launch the packaged example output in `CMEntViewer`, publish the
Shiny port from the container and bind the app to `0.0.0.0`:

``` bash
docker run --rm -p 3838:3838 \
  vlemonidis/cment:latest \
  launchCMEntViewer \
  --output_prefix /CMEnt/inst/extdata/example_output \
  --launch_browser FALSE \
  --host 0.0.0.0 \
  --port 3838
```

Then open `http://localhost:3838` in your browser.

If you want to view outputs created on your machine, mount them into
the container. Keeping the R library and cache on Docker volumes avoids
re-installing on-demand packages and cached annotations on every run:

``` bash
docker volume create cment-r-lib
docker volume create cment-r-cache

docker run --rm -p 3838:3838 \
  -v "$PWD:/work" -w /work \
  -v cment-r-lib:/usr/local/lib/R/site-library \
  -v cment-r-cache:/home/root/.cache/R \
  vlemonidis/cment:latest \
  launchCMEntViewer \
  --output_prefix /work/path/to/output_prefix \
  --launch_browser FALSE \
  --host 0.0.0.0 \
  --port 3838
```

## Usage

Here’s a basic example of how to use CMEnt:

    suppressPackageStartupMessages(library(CMEnt))
    beta <- loadExampleInputDataChr5And11("beta")
    dmps <- loadExampleInputDataChr5And11("dmps")
    pheno <- loadExampleInputDataChr5And11("pheno")
    array_type <- loadExampleInputDataChr5And11("array_type")
    # Find DMRs using parallel processing
    dmrs <- findDMRsFromSeeds(
        beta = beta,
        seeds = dmps,
        pheno = pheno,
        sample_group_col = "Sample_Group",
        array = array_type,  # Specify array platform
        min_cpgs = 3,
        min_seeds = 2,
        njobs = 2  # Use parallel processing
    )

    # View summary of DMR statistics
    summary(dmrs)

## Citation

If you use CMEnt in your research, please cite:

    To cite package 'CMEnt' in publications use:

      Lemonidis et al (2026). CMEnt: Characterization of Methylation using positional Entanglement. R package version 0.99.0.
      University of Antwerp, Antwerp, Belgium.

    A BibTeX entry for LaTeX users is:

      @Manual{
        title = {From DMPs to DMRs: CMEnt, Characterization of Methylation using positional Entanglement},
        author = {Vasileios Lemonidis, Joe Ibrahim, Joris R Vermeesch, Timon Vandamme, Guy Van Camp, Ken op de Beeck},
        year = {2026},
        note = {R package version 0.99.0},
        url = {https://github.com/CMG-UA/CMEnt},
      }

## License

This project is licensed under the GPL-2 License - see the
[LICENSE.md](LICENSE.md) file for details.

## Authors

- **Vasileios Lemonidis** - *Initial work* - \[University of Antwerp\]

## Acknowledgments

- Center of Medical Genetics, University of Antwerp
- [Genetics of Cancer and Hearing Loss Research
  Group](https://www.linkedin.com/company/genetics-of-cancer-and-hearing-loss)
