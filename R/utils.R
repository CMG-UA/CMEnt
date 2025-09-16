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
    if (endsWith(beta.file, 'gz')) {
        conn <- gzfile(beta.file, 'r')
    } else {
        conn <- file(beta.file, 'r')
    }
    file.beta.col.names <- scan(conn, sep = '\t',
                               what = character(),
                               nlines = 1,
                               quiet = TRUE)
    close(conn)
    
    if (is.tabix) {
        file.beta.col.names <- file.beta.col.names[7:length(file.beta.col.names)]
    }
    
    if (is.null(beta.col.names)) {
        cols.inds <- seq(2, length(file.beta.col.names))
        beta.col.names <- file.beta.col.names[-1]
    } else {
        cols.inds <- match(beta.col.names, file.beta.col.names)
        cols.inds <- cols.inds[!is.na(cols.inds), drop=FALSE]
        if (length(cols.inds) == 0) {
            stop("Beta file does not contain any phenotype rownames as column name. ",
                 "First 5 supplied phenotype rownames: ",
                 paste(beta.col.names[seq_len(min(5, length(beta.col.names)))],
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

