#' DMRSegal utility functions
#' @keywords internal

#' Get beta value column names and indices
#'
#' @param beta.file Path to beta value file
#' @param beta.col.names Column names to match
#' @param is.tabix Whether the file is a tabix file
#' @return List with beta column information
#' @keywords internal
.get.beta.col.names.and.inds <- function(beta.file, beta.col.names=NULL, is.tabix=FALSE) {
    if (endsWith(beta.file, 'gz'))
        conn <- gzfile(beta.file, 'r')
    else
        conn <- file(beta.file, 'r')
    file.beta.col.names <- scan(conn, sep = '\t',
                               what = character(),
                               nlines = 1,
                               quiet = TRUE)
    close(conn)
    
    if (is.tabix) {
        file.beta.col.names <- file.beta.col.names[7:length(file.beta.col.names)]
    }
    
    if (is.null(beta.col.names)) {
        beta.col.names <- file.beta.col.names[-1]
        cols.inds <- seq(2, length(beta.col.names))
    } else {
        cols.inds <- match(beta.col.names, file.beta.col.names)
        cols.inds <- cols.inds[!is.na(cols.inds), drop=FALSE]
        if (length(cols.inds) == 0) {
            stop("Beta file does not contain any phenotype rownames as column name. ",
                 "First 5 supplied phenotype rownames: ",
                 paste(beta.col.names[1:min(5, length(beta.col.names))], 
                       collapse=','))
        }
        beta.col.names <- file.beta.col.names[cols.inds]
    }
    
    ret <- list(
        beta.col.names = beta.col.names,
        beta.col.inds = cols.inds,
        file.beta.col.names = file.beta.col.names[-1]
    )
    invisible(ret)
}

#' Subset beta values from file
#'
#' @param beta.file Path to beta value file
#' @param sites Sites to subset
#' @param beta.row.names Row names
#' @param beta.col.names Column names
#' @return Matrix of beta values
#' @keywords internal
.subset.beta <- function(beta.file,
                        sites,
                        beta.row.names = NULL,
                        beta.col.names = NULL) {
    if (is.null(beta.row.names))
        beta.row.names <- unlist(data.table::fread(
            file = beta.file,
            select = 1,
            header = TRUE
        ))
    
    ret <- .get.beta.col.names.and.inds(beta.file, beta.col.names)
    beta.col.names <- ret$beta.col.names
    cols.inds <- ret$beta.col.inds
    
    sites.inds <- which(beta.row.names %in% sites)
    sites.inds.steps <- diff(c(-1, sites.inds)) - 1
    
    if (endsWith(beta.file, 'gz'))
        conn <- gzfile(beta.file, 'r')
    else
        conn <- file(beta.file, 'r')
    
    beta.sites <- lapply(sites.inds.steps, function(step) {
        l <- scan(
            conn,
            what = character(),
            skip = step,
            nlines = 1,
            quiet = TRUE
        )
        as.numeric(l[cols.inds])
    })
    close(conn)
    
    beta.sites <- do.call(rbind, beta.sites)
    rownames(beta.sites) <- sites
    colnames(beta.sites) <- beta.col.names
    beta.sites
}

#' Test connectivity between CpG sites
#'
#' @param site1.beta Beta values for first site
#' @param site2.beta Beta values for second site
#' @param sample.groups Sample group factor
#' @param max.pval Maximum p-value threshold
#' @param casecontrol Case-control vector
#' @param min.delta_beta Minimum delta beta
#' @param extreme.verbosity Verbose output flag
#' @return List with connectivity test results
#' @keywords internal
.test.connectivity <- function(site1.beta,
                             site2.beta,
                             sample.groups,
                             max.pval,
                             casecontrol = NULL,
                             min.delta_beta = 0,
                             extreme.verbosity = FALSE) {
    pval <- 0
    delta_beta <- NULL
    max.pval.corrected <- max.pval / length(unique(sample.groups))
    
    for (g in levels(sample.groups)) {
        if (sum(sample.groups == g) < 3)
            next
        
        op <- options(warn = 2)$warn
        corr.ret <- try(psych::corr.test(
            site1.beta[sample.groups == g],
            site2.beta[sample.groups == g]
        ))
        options(warn = op)
        
        if (inherits(corr.ret, "try-error")) {
            if (extreme.verbosity) {
                message("Error in correlation test for group ", g)
            }
            return(list(FALSE, pval, delta_beta,
                       failing = g,
                       reason = 'error occurred'))
        }
        
        r <- max(pval, corr.ret$p)
        if (is.null(r) || is.na(r)) {
            return(list(FALSE, pval, delta_beta,
                       failing = g,
                       reason = 'na pval'))
        }
        
        pval <- r
        if (pval > max.pval.corrected) {
            return(list(FALSE, pval, delta_beta,
                       failing = g,
                       reason = 'pval>max.pval (corrected)'))
        }
    }
    
    if (!is.null(casecontrol) && (min.delta_beta > 0)) {
        if (length(casecontrol) != length(site2.beta)) {
            return(list(FALSE, pval, delta_beta,
                       failing = 'all',
                       reason = 'casecontrol length mismatch'))
        }
        
        delta_beta <- abs(mean(site2.beta[casecontrol]) - 
                         mean(site2.beta[!casecontrol]))
        if (delta_beta < min.delta_beta) {
            return(list(FALSE, pval, delta_beta,
                       failing = 'all',
                       reason = 'delta_beta < min_delta_beta'))
        }
    }
    
    list(TRUE, pval, delta_beta)
}