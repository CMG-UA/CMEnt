options(Ncpus = max(1, parallel::detectCores() - 1))
# Ensure binary packages are preferred
options(HTTPUserAgent = sprintf(
    "R/%s R (%s)",
    getRversion(),
    paste(getRversion(), R.version["platform"], R.version["arch"], R.version["os"])
))
options(warn = 1)
options(error = function() {
    calls <- sys.calls()
    if (length(calls) >= 2L) {
        sink(stderr())
        on.exit(sink(NULL))
        cat("Backtrace:\n")
        calls <- rev(calls[-length(calls)])
        for (i in seq_along(calls)) {
            cat(i, ": ", deparse(calls[[i]], nlines = 1L), "\n", sep = "")
        }
    }
    if (!interactive()) {
        q(status = 1)
    }
})

try_wrapper <- function(expr, attempts = 3, delay = 5) {
    expr_sub <- substitute(expr)
    for (attempt in seq_len(attempts)) {
        success <- FALSE
        value <- tryCatch({
            value <- eval(expr_sub, envir = parent.frame())
            success <- TRUE
            value
        }, error = function(e) {
            if (attempt == attempts) {
                stop("Error after ", attempts, " attempts: ", conditionMessage(e), call. = FALSE)
            } else {
                message("Attempt ", attempt, " failed: ", conditionMessage(e), ". Retrying...")
                Sys.sleep(delay)
                NULL
            }
        })
        if (success) return(value)
    }
}


install_deps <- function(packages = NULL,
                         bioc_packages = NULL,
                         git_packages = NULL,
                         attempts = 3) {
    if (is.null(packages)) packages <- c()
    if (!is.null(bioc_packages)) {
        bioc_packages <- as.list(bioc_packages)
        packages <- c(packages, "BiocManager")
        bioc_packages <- append(bioc_packages, "remotes", after = 2)
    }

    # Cache installed packages for performance
    installed_pkgs <- installed.packages()[, "Package"]

    to_install <- packages[!(packages %in% installed_pkgs)]
    if (length(to_install) > 0) {
        try_wrapper({
            install.packages(to_install, dependencies = TRUE)
        }, attempts = attempts)
        installed_pkgs <- installed.packages()[, "Package"]
    }
    errored <- to_install[!(to_install %in% installed_pkgs)]
    if (length(errored) > 0) {
        stop("The following packages could not be installed: ", paste(errored, collapse = ", "))
    }

    if (!is.null(bioc_packages)) {
        options(repos = BiocManager::repositories())

        bioc_packages_names <- sapply(bioc_packages, function(x) if (is.list(x)) x[[1]] else x)
        bioc_to_install <- bioc_packages[
            !(bioc_packages_names %in% installed_pkgs) | sapply(bioc_packages, function(x) {
                is.list(x) && isTRUE(x[["force"]])
            })
        ]
        if (length(bioc_to_install) > 0) {
            without_extra_args <- c()
            for (pkg in bioc_to_install) {
                if (is.list(pkg)) {
                    try_wrapper({
                        # Keep the target package name explicit; modifyList() drops unnamed entries.
                        install_args <- c(
                            list(pkgs = pkg[[1]]),
                            utils::modifyList(
                                list(ask = FALSE, update = FALSE),
                                pkg[-1]
                            )
                        )
                        do.call(BiocManager::install, install_args)
                    }, attempts = attempts)
                } else {
                    without_extra_args <- c(without_extra_args, pkg)
                }
            }
            if (length(without_extra_args) > 0) {
                try_wrapper({
                    BiocManager::install(without_extra_args, ask = FALSE, update = FALSE)
                }, attempts = attempts)
            }
            installed_pkgs <- installed.packages()[, "Package"]
        }
        errored <- bioc_packages_names[!(bioc_packages_names %in% installed_pkgs)]
        if (length(errored) > 0) {
            stop("The following Bioconductor packages could not be installed: ", paste(errored, collapse = ", "))
        }
    }

    if (!is.null(git_packages)) {
        if (!"BiocManager" %in% installed_pkgs) {
            try_wrapper({
                install.packages("BiocManager", dependencies = TRUE)
            }, attempts = attempts)
            installed_pkgs <- installed.packages()[, "Package"]
        }

        # Ensure remotes can resolve Bioconductor dependencies when installing from GitHub.
        options(repos = BiocManager::repositories())

        if (!"remotes" %in% installed_pkgs) {
            try_wrapper({
                install.packages("remotes", dependencies = TRUE)
            }, attempts = attempts)
            installed_pkgs <- installed.packages()[, "Package"]
        }

        ch_packages <- sub("@.*$", "", git_packages)
        ch_packages <- sub(".*/", "", ch_packages)
        git_to_install <- git_packages[!(ch_packages %in% installed_pkgs)]
        if (length(git_to_install) > 0) {
            for (git_pkg in git_to_install) {
                try_wrapper({
                    remotes::install_github(
                        git_pkg,
                        dependencies = NA,
                        upgrade = "never",
                        build_vignettes = FALSE
                    )
                }, attempts = attempts)
            }
            installed_pkgs <- installed.packages()[, "Package"]
        }

        errored <- ch_packages[!(ch_packages %in% installed_pkgs)]
        if (length(errored) > 0) {
            stop("The following GitHub packages could not be installed: ", paste(errored, collapse = ", "))
        }
    }
}
