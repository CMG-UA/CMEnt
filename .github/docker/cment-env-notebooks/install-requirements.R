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

install_deps(packages = cran, bioc_packages = bioc)