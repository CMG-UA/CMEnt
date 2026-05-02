read_pkgs <- function(path) {
  if (!file.exists(path)) {
    return(character())
  }

  x <- readLines(path, warn = FALSE)
  x <- trimws(x)
  x <- x[nzchar(x)]
  x <- x[!startsWith(x, "#")]
  sort(unique(x))
}

cran <- read_pkgs("cran-requirements.txt")
bioc <- read_pkgs("bioc-requirements.txt")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

if (length(cran) > 0L) {
  install.packages(
    cran,
    repos = "https://cloud.r-project.org"
  )
}

if (length(bioc) > 0L) {
  BiocManager::install(
    bioc,
    ask = FALSE,
    update = FALSE
  )
}