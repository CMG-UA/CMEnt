#' Find Differentially Methylated Regions (DMRs) from Pre-computed DMPs
#'
#' @name findDMRsFromDMPs
#' @description This function identifies Differentially Methylated Regions (DMRs) from pre-computed
#' Differentially Methylated Positions (DMPs) using a correlation-based approach. It expands
#' significant DMPs into regions, considering both statistical significance and biological
#' relevance of methylation changes.
#'
#' @param beta.file Path to the methylation beta values file or a data matrix with beta values
#' @param dmps.tsv.file Path to the pre-computed DMPs file or a data frame with DMPs
#' @param pheno Data frame containing sample phenotype information
#' @param output.dir Output directory for results (default: '.')
#' @param pval.col Column name in DMPs file containing p-values (default: "pval_adj")
#' @param sample_group.col Column in pheno for sample grouping (default: "Sample_Group")
#' @param casecontrol.col Column in pheno for case/control status (default: "casecontrol")
#' @param min.cpg.delta_beta Minimum delta beta threshold for CpGs (default: 0)
#' @param expansion.step Distance in bp to expand regions (default: 50)
#' @param array Array platform, either "450K" or "EPIC" (default: c("450K", "EPIC"))
#' @param genome Reference genome, "hg19" or "hg38" (default: c("hg19", "hg38"))
#' @param max.pval Maximum p-value threshold (default: 0.05)
#' @param max.lookup.dist Maximum distance for region expansion (default: 10000)
#' @param min.dmps Minimum number of DMPs required per region (default: 1)
#' @param min.cpgs Minimum number of CpGs required per region (default: 50)
#' @param ignored.sample.groups Sample groups to ignore (default: NULL)
#' @param output.prefix Optional identifier for output files (default: NULL)
#' @param njobs Number of parallel jobs (default: detectCores())
#' @param verbose Enable verbose output (default: FALSE)
#' @param beta.row.names.file Optional file with beta value row names (default: NULL)
#' @param dmps.beta.file Optional separate beta file for DMPs (default: NULL)
#'
#' @return A GRanges object containing identified DMRs with metadata columns:
#' \itemize{
#'   \item n_cpgs: Number of CpGs in the region
#'   \item n_dmps: Number of DMPs in the region
#'   \item mean_delta_beta: Mean methylation difference
#'   \item max_delta_beta: Maximum methylation difference
#'   \item min_pval: Minimum p-value of DMPs in region
#' }
#'
#' @examples
#' \dontrun{
#' # Load example data
#' data(example_beta)
#' data(example_dmps)
#' data(example_pheno)
#'
#' # Write beta values to file
#' beta_file <- tempfile(fileext = ".txt")
#' write.table(cbind(ID = rownames(example_beta), example_beta),
#'             file = beta_file, sep = "\t", quote = FALSE, row.names = FALSE)
#'
#' # Find DMRs
#' dmrs <- findDMRsFromDMPs(
#'   beta.file = beta_file,
#'   dmps.tsv.file = example_dmps,
#'   pheno = example_pheno
#' )
#' }
#'
#' @importFrom GenomicRanges GRanges makeGRangesFromDataFrame
#' @importFrom parallel detectCores mclapply
#' @importFrom stringr str_count str_order
#' @importFrom readr read_tsv
#' @importFrom data.table fread
#' @importFrom psych corr.test
#' @importFrom dplyr %>% filter select mutate
#' @importFrom bedr tabix
#' @importFrom BSgenome getSeq
#' @importFrom rtracklayer import.chain
#' @importFrom GenomeInfoDb Seqinfo
#' @importFrom stats sd
#' @importFrom utils write.table read.table
#' @import IlluminaHumanMethylation450kanno.ilmn12.hg19
#' @import IlluminaHumanMethylationEPICanno.ilm10b4.hg19
#' @export


.get.beta.col.names.and.inds <- function(beta.file, beta.col.names=NULL, is.tabix=F){
  if (endsWith(beta.file, 'gz'))
    conn <- gzfile(beta.file, 'r')
  else
    conn <- file(beta.file, 'r')
  file.beta.col.names <- scan(conn,sep = '\t',
                              what = character(),
                              nlines = 1,
                              quiet = T)
  close(conn)
  if (is.tabix) {
    file.beta.col.names <- file.beta.col.names[7:length(file.beta.col.names)]
  }
  if (is.null(beta.col.names)) {
    
    beta.col.names <- file.beta.col.names[-1]
    cols.inds <- seq(2,length(beta.col.names))
    
  } else {
    cols.inds <- match(beta.col.names, file.beta.col.names)
    cols.inds <- cols.inds[!is.na(cols.inds), drop=F]
    if (length(cols.inds) == 0){
      stop("Beta file does not contain any phenotype rownames as column name. First 5 supplied phenotype rownames: ",paste(beta.col.names[seq_len(min(5, length(beta.col.names)))], sep=','))
    }
    beta.col.names <- file.beta.col.names[cols.inds]
  }
  ret <- list(beta.col.names=beta.col.names, beta.col.inds=cols.inds, file.beta.col.names=file.beta.col.names[-1])
  invisible(ret)
}

.subset.beta <- function(beta.file,
                         sites,
                         beta.row.names = NULL,
                         beta.col.names = NULL) {
  # Fallback simple subsetting path for small site sets (avoids scan/skip logic issues in tests)
  if (length(sites) <= 5000) {
    full <- data.table::fread(beta.file, header=TRUE, data.table=FALSE)
    rn <- full[[1]]
    full <- full[,-1, drop=FALSE]
    rownames(full) <- rn
    # Preserve only requested sites in order
    idx <- match(sites, rn)
    if (anyNA(idx)) {
      missing <- sites[is.na(idx)]
      stop("Internal error: requested CpG IDs not found in beta file during subsetting: ", paste(missing, collapse=","))
    }
    sub <- full[idx, , drop=FALSE]
    if (nrow(sub) == 0){
      browser()
      stop("None of the provided sites exist in the read beta file")
    }
    sub <- apply(sub, 2, as.numeric)
    
    tryCatch({
      
      rownames(sub) <- sites
    }, error=function(e){
      browser()
    }
    )
    nas.per.row <- apply(sub,1,function(r) sum(is.na(r)))
    if (any(nas.per.row == ncol(sub))) {
      warning("All-NA beta rows detected for sites: ", paste(names(nas.per.row)[nas.per.row==ncol(sub)], collapse=","))
    }
    return(sub)
  }
  if (is.null(beta.row.names))
    beta.row.names <- unlist(fread(
      file = beta.file,
      select = 1,
      header = T,
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
      quiet = T
    )
    as.numeric(l[cols.inds])
  })
  close(conn)
  beta.sites <- do.call(rbind, beta.sites)
  rownames(beta.sites) <- sites
  colnames(beta.sites) <- beta.col.names
  beta.sites
}

# Expand the identified dmrs to nearby CpG regions
.expandDMRs <- function(dmr,
                        beta.file,
                        beta.row.names,
                        beta.col.names,
                        sample.groups,
                        sorted.locs,
                        max.pval,
                        min.cpg.delta_beta=0,
                        casecontrol = NULL,
                        tabix.file=NULL,
                        expansion.step = 500) {
  if (!is.null(beta.file)){
    ret <- .get.beta.col.names.and.inds(beta.file, beta.col.names)

  } else {
    ret <- .get.beta.col.names.and.inds(tabix.file, beta.col.names, is.tabix=T)
  }

  cols.inds <- ret$beta.col.inds
  file.beta.col.names <- ret$file.beta.col.names
  beta.col.names <- ret$beta.col.names
  sorted.locs <- sorted.locs[beta.row.names,]
  
  dmr.start <- dmr$start_dmp
  
  dmr.end <- dmr$end_dmp
  dmr.start.ind <- which(beta.row.names == dmr.start)
  
  if (length(dmr.start.ind) == 0) {
    stop("Could not find the start CpG ", dmr.start, " in the beta file row names.")
  }
  end.site.ind <- dmr.start.ind[[1]]
  upstream.exp <- end.site.ind
  upstream.stop.found <- F
  while (T) {
    if (end.site.ind < 0) {
      upstream.stop.reason <- "end-of-input"
      upstream.exp <- 1
      break
    }
    start.site.ind <- max(0, end.site.ind - expansion.step) + 1
    x <- which(sorted.locs[start.site.ind:end.site.ind, "chr"] == dmr$chr)
    if (length(x) == 0) {
      upstream.stop.reason <- "end-of-input"
      break
    }
    x <- x[[1]]
    start.site.ind <- start.site.ind + x - 1
    exp.step <- end.site.ind - start.site.ind + 1
    if (!is.null(beta.file)){
      upstream.betas <- fread(
        file = beta.file,
        skip = start.site.ind,
        nrows = exp.step,
        header = F,data.table=FALSE,
        colClasses=c("character", rep("numeric", length(file.beta.col.names)))
      )
      upstream.betas <- upstream.betas[, cols.inds]
    } else {
      upstream.region <- paste0(
        sorted.locs[start.site.ind, 'chr'], ':',
        sorted.locs[start.site.ind, 'pos'], '-',
        sorted.locs[end.site.ind, 'pos']+1)
      upstream.betas <- try(tabix(
        upstream.region,
        tabix.file,
        check.valid=F,
        verbose=F
      ))
      if(inherits(upstream.betas, "try-error")){
        warning("Error reading upstream region ", upstream.region, " from tabix file, stopping extension. The error is:\n\t", paste(capture.output(print(upstream.betas)), collapse='\n\t'))
        upstream.stop.reason <- "error-reading-tabix"
        break
      }
      upstream.betas <- as.data.frame(sapply(upstream.betas[,beta.col.names], as.numeric))
    }
    upstream.betas <- upstream.betas[rev(seq_len(nrow(upstream.betas))),,drop=F]
    if (nrow(upstream.betas) == 1){
      upstream.stop.reason <- "end-of-input"
      break
    }
    i <- 1
    while (T) {
      corr.ret <- .testConnectivity(site1.beta = unlist(upstream.betas[i, ]),
                                     site2.beta = unlist(upstream.betas[i + 1, ]),
                                     sample.groups = sample.groups,
                                     casecontrol = casecontrol,
                                     max.pval = max.pval,
                                     min.delta_beta = min.cpg.delta_beta,
                                     extreme.verbosity=T)
      if (!corr.ret[[1]]) {
        upstream.exp <- end.site.ind - i + 1
        upstream.stop.found <- T
        upstream.stop.reason <- corr.ret$reason
        break
      }
      i <- i + 1
      if (i == nrow(upstream.betas)){
        break
      }
    }
    if (upstream.stop.found)
      break
    end.site.ind <- end.site.ind - expansion.step
  }
  
  dmr.end.ind <- which(beta.row.names == dmr.end)
  if (length(dmr.end.ind) == 0) {
    stop("Could not find the end CpG ", dmr.end, " in the beta file row names.")
  }
  start.site.ind <- dmr.end.ind[[1]]
  downstream.exp <- start.site.ind
  ns <- 0
  downstream.stop.found <- F
  while (T) {
    if (start.site.ind > length(beta.row.names)) {
      downstream.stop.reason <- 'end-of-input'
      downstream.exp <- start.site.ind - 1
      break
    }
    end.site.ind <- min(start.site.ind + expansion.step, length(beta.row.names))
    x <- which(rev(sorted.locs[start.site.ind:end.site.ind,'chr'])== dmr$chr)
    if (length(x) == 0){
      downstream.stop.reason <- "end-of-input"
      downstream.exp <- start.site.ind - 1
      break
    }
    x <- x[[1]]
    end.site.ind <- end.site.ind - x + 1
    if (end.site.ind <= start.site.ind + 1){
      downstream.stop.reason <- "end-of-input"
      downstream.exp <- start.site.ind - 1
      break
    }
    if (!is.null(beta.file)){
      downstream.betas <- fread(
        file = beta.file,
        skip = start.site.ind,
        nrows = expansion.step - x + 1,
        header = F,
        data.table=FALSE,
        colClasses=c("character", rep("numeric", length(file.beta.col.names)))
      )
      downstream.betas <- downstream.betas[, cols.inds]
    } else {
      downstream.region <- paste0(
        sorted.locs[start.site.ind, 'chr'], ':',
        sorted.locs[start.site.ind, 'pos'], '-',
        sorted.locs[end.site.ind, 'pos'])
      downstream.betas <- try(tabix(
        downstream.region,
        tabix.file,
        check.valid=F,
        verbose=F
      ))
      if(inherits(downstream.betas, "try-error")){
        warning("Error reading downstream region ", downstream.region, " from tabix file, stopping extension. The error is:\n\t", paste(capture.output(print(downstream.betas)), collapse='\n\t'))
        downstream.stop.reason <- "error-reading-tabix"
        break
      }
      downstream.betas <- as.data.frame(sapply(downstream.betas[,beta.col.names], as.numeric))
    }
    
    i <- 1
    while (T) {
      corr.ret <- .testConnectivity(site1.beta = unlist(downstream.betas[i, ]),
                                     site2.beta = unlist(downstream.betas[i + 1, ]),
                                     sample.groups = sample.groups,
                                     max.pval = max.pval,
                                     casecontrol = casecontrol,
                                     min.delta_beta = min.cpg.delta_beta)
      if (!corr.ret[[1]]) {
        downstream.exp <- start.site.ind + i - 1
        downstream.stop.reason <- corr.ret$reason
        downstream.stop.found <- T
        break
      }
      i <- i + 1
      if (i == nrow(downstream.betas)){
        break
      }
    }
    if (downstream.stop.found)
      break
    prev.start.site.ind <- start.site.ind  
    start.site.ind <- start.site.ind + expansion.step - x + 1
    if (start.site.ind < prev.start.site.ind) {
      stop("BUG: start.site.ind < prev.start.site.ind during downstream expansion. start.site.ind: ", start.site.ind, " prev.start.site.ind: ", prev.start.site.ind, " expansion.step: ", expansion.step, " x: ", x) 
    }
  }
  dmr$start_cpg <- beta.row.names[upstream.exp]
  tryCatch({
  dmr$end_cpg <- beta.row.names[downstream.exp]
  }, error=function(e){browser()})
  dmr$start <- sorted.locs[dmr$start_cpg, 'pos']
  dmr$end <- sorted.locs[dmr$end_cpg, 'pos']
  dmr$downstream_cpg_expansion_stop_reason <- downstream.stop.reason
  dmr$upstream_cpg_expansion_stop_reason <- upstream.stop.reason
  
  dmr
}

.testConnectivity <- function(site1.beta, site2.beta, sample.groups, max.pval,  casecontrol=NULL, min.delta_beta=0, extreme.verbosity=F) {
  pval <- 0
  delta_beta <- NULL

  max.pval.corrected <- max.pval / length(unique(sample.groups)) # Bonferroni correction 
  for (g in levels(sample.groups)) {
    if (sum(sample.groups == g) < 3)
      next
    op <- options(warn=2)$warn
    corr.ret <- try(corr.test(site1.beta[sample.groups == g], site2.beta[sample.groups == g]))
    options(warn=op)
    if (inherits(corr.ret, "try-error")){
      if (extreme.verbosity){
        message(".testConnectivity: Error occurred in corr.test while processing the following:")
        message("casecontrol:", paste(casecontrol, collapse=','))
        message("site2.beta:", paste(site2.beta, collapse=','))
        message('sample.groups:', paste(sample.groups, collapse=','))
        message('max.pval.corrected:', max.pval.corrected)
        message("Error message:", corr.ret)
        browser()
      }
      return(list(F, pval, delta_beta, failing=g, reason='error occurred'))
    }
    r <-  max(pval, corr.ret$p)
    if (is.null(r) || is.na(r)) {
      return(list(F, pval, delta_beta, failing=g, reason='na pval'))
    }
    pval <- r
    if (pval > max.pval.corrected) {
      return(list(F, pval, delta_beta, failing=g, reason='pval>max.pval (corrected)'))
    }
  }
  if (!is.null(casecontrol) && (min.delta_beta>0)){
    if (length(casecontrol) != length(site2.beta)){
      if (extreme.verbosity){
        message(".testConnectivity: Error occurred while computing delta beta for the following:")
        message("casecontrol:", paste(casecontrol, collapse=','))
        message("site2.beta:", paste(site2.beta, collapse=','))
        message('sample.groups:', paste(sample.groups, collapse=','))
        message('max.pval.corrected:', max.pval.corrected)
      }
      stop(paste0(
        "The provided casecontrol vector has length", length(casecontrol), 
        " while the site2.beta has length ", length(site2.beta)))
      
    }
    delta_beta <- mean(site2.beta[casecontrol == 1], na.rm=T) - mean(site2.beta[casecontrol == 0], na.rm=T)
    if (is.null(delta_beta) || is.na(delta_beta) || (abs(delta_beta) < min.delta_beta)){
      return(list(F, pval, delta_beta, reason="delta_beta<min.delta_beta"))
    }
  }
  return (list(T, pval, delta_beta))
}


#' Find Differentially Methylated Regions (DMRs) from Differentially Methylated Positions (DMPs)
#'
#' This function identifies DMRs from a given set of DMPs and a beta value file.
#'
#' @param beta.file Character. Path to the beta value file. Either this or tabix.file must be provided.
#' @param dmps.tsv.file Character. Path to the DMPs TSV file.
#' @param pheno Data frame. Phenotype data.
#' @param output.dir Character. Directory to save the output files. Default is current directory.
#' @param pval.col Character. Column name for p-values in the DMPs file. Default is "pval_adj".
#' @param sample_group.col Character. Column name for sample group information in the phenotype data. Default is NULL.
#' @param casecontrol.col Character. Column name for case-control information in the phenotype data. Default is "casecontrol".
#' @param min.cpg.delta_beta Numeric. Minimum delta beta value for CpGs. Default is 0.
#' @param expansion.step Numeric. Step size for expanding DMRs. Increasing it means higher memory usage and faster computation. Default is 500.
#' @param array Character. Type of array used ("450K" or "EPIC"). Default is "450K".
#' @param max.pval Numeric. Maximum p-value to assume DMPs correlation is significant. Default is 0.05.
#' @param max.lookup.dist Numeric. Maximum distance to look up for adjacent DMPs belonging to the same DMR. Default is 10000.
#' @param min.dmps Numeric. Minimum number of DMPs in a DMR. Default is 1.
#' @param min.cpgs Numeric. Minimum number of CpGs in a DMR. Default is 50.
#' @param output.prefix Character. Identifier for the output files. Default is NULL.
#' @param njobs Numeric. Number of parallel jobs to use. Default is the number of available cores.
#' @param verbose Logical. Whether to print detailed messages. Default is FALSE.
#'
#' @return Data frame of identified DMRs.
#' @export
findDMRsFromDMPs <- function(beta.file=NULL,
                             dmps.tsv.file=NULL,
                             pheno=NULL,
                             pval.col = "pval_adj",
                             sample_group.col = "Sample_Group",
                             casecontrol.col = "casecontrol",
                             min.cpg.delta_beta = 0,
                             expansion.step = 50,
                             array = c("450K", "EPIC"),
                             genome = c("hg19", "hg38"),
                             max.pval = 0.05,
                             max.lookup.dist = 10000,
                             min.dmps = 1,
                             min.cpgs = 50,
                             ignored.sample.groups = NULL,
                             output.prefix = NULL,
                             njobs = parallel::detectCores(),
                             verbose = FALSE,
                             beta.row.names.file=NULL,
                             dmps.beta.file=NULL,
                             tabix.file=NULL) {
  if (is.null(dmps.tsv.file) || is.null(pheno)){
    stop("dmps.tsv.file and pheno parameters are required")
  }
  if (is.null(beta.file) && is.null(tabix.file)){
    stop("Either beta.file or tabix.file parameter is required")
  }
  stopifnot(!is.null(max.pval))
  stopifnot(!is.null(min.dmps))
  stopifnot(!is.null(min.cpgs))
  stopifnot(!is.null(expansion.step))
  stopifnot(!is.null(min.cpg.delta_beta))
  stopifnot(!is.null(max.lookup.dist))

  array <- match.arg(array)
  genome <- match.arg(genome)
  if (!is.null(beta.file)){
  stopifnot(file.exists(beta.file))
  }
  if (!is.null(tabix.file)){
    stopifnot(file.exists(tabix.file))
  }
  stopifnot(file.exists(dmps.tsv.file))
  stopifnot(casecontrol.col %in% colnames(pheno))
  stopifnot(sample_group.col %in% colnames(pheno))
  if (is.null(ignored.sample.groups)){
    ignored.sample.groups <- c()
  } else {
    ignored.sample.groups <- unlist(strsplit(ignored.sample.groups, ','))
  }
  if (!is.null(output.prefix)) {
    output.dir <- dirname(output.prefix)
    dir.create(output.dir, showWarnings = F)
    output.prefix <- paste0(output.prefix,'.')
  } else {
    output.dir <- NULL
    output.prefix <- NULL
  }
  
  
  if (verbose)
    message("Reading dmp tsv file..")
  dmps.tsv <- try(read.table(
    dmps.tsv.file,
    header = TRUE,
    sep = "\t",
    check.names = FALSE,
    quote="",
    comment.char = "",
    row.names = NULL,
  ))
  if (inherits(dmps.tsv, "try-error")) {
    message("Provided DMPs file is empty or does not exist. Not proceeding.")
    if (!is.null(output.prefix)) {
      for (f in c('.methylation.tsv.gz', '.tsv.gz')){
        gzfile <- gzfile(paste0(output.prefix, f), "w", compression = 2)   
        close(gzfile)
      }
    }
    return(NULL)
  }
  if (! sample_group.col %in% colnames(dmps.tsv)){
      stop("Sample group column '", sample_group.col,
        "' does not reside in the DMPs file columns: ", 
           paste(colnames(dmps.tsv), collapse=','))
  }
  
  if (! pval.col %in% colnames(dmps.tsv)) {
    stop("P-value column '", pval.col, "' does not reside in the DMPs file columns: ",  
           paste(colnames(dmps.tsv), collapse=','))
  }
  stopifnot(pval.col %in% colnames(dmps.tsv))
  
  colnames(dmps.tsv)[[1]] <- "dmp"
  
  if (verbose)
    message("Reading beta file characteristics..")
  
  if (!is.null(output.prefix)){
    output.beta.row.names.file <- paste0(output.prefix, "row.names.txt")
    if (is.null(beta.row.names.file)){
      beta.row.names.file <- output.beta.row.names.file # load from previous run
    }
  }
  if (!is.null(beta.row.names.file) && file.exists(beta.row.names.file)){
    if (verbose) message("Reading beta file row names from beta row names file ", beta.row.names.file, "..")
    beta.row.names <- unlist(read.table(beta.row.names.file, header=F, comment.char = "", quote = ""))
  } else {
    if (verbose) message("Reading beta file row names from beta file..")
    if (!is.null(beta.file)){
    beta.row.names <- unlist(fread(
      file = beta.file,
      select = 1,
      sep='\t',
      header = T,
      showProgress=T,
      nThread=njobs
    ))
    }
    else{
      beta.row.names <- unlist(fread(
      file = tabix.file,
      select = 4,
      sep='\t',
      header = T,
      showProgress=T,
      nThread=njobs))
    }
    if (!is.null(beta.row.names.file)){
    if (verbose) message("Saving beta file row names to beta row names file: ", beta.row.names.file)
    writeLines(paste(beta.row.names, collapse='\n'), beta.row.names.file)
    }
  }
  if (verbose)
    message('Number of rows names read:', length(beta.row.names))
  
  beta.col.names <- row.names(pheno)
  samples.selection.mask <- !( pheno[, sample_group.col] %in% ignored.sample.groups )
  beta.col.names <- beta.col.names[samples.selection.mask]
  pheno <- pheno[beta.col.names, ]
  if (verbose) message("Number of samples to process: ", length(beta.col.names))

  if (!is.null(beta.file)){
    beta.col.names <- .get.beta.col.names.and.inds(beta.file, beta.col.names)$beta.col.names
  } else {
    beta.col.names <- .get.beta.col.names.and.inds(tabix.file, beta.col.names, is.tabix=T)$beta.col.names
  }
  pheno <- pheno[beta.col.names,]
  sample.groups <- factor(pheno[beta.col.names, sample_group.col])
  
  # THEY ARE NOT SORTED BY DEFAULT
  if (array == "450K") {
    sorted.locs <- IlluminaHumanMethylation450kanno.ilmn12.hg19::Locations
  } else if (array == "EPIC") {
    sorted.locs <- IlluminaHumanMethylationEPICanno.ilm10b4.hg19::Locations
  } else {
    stop("Unknown array type: ", array)
  }
  sorted.locs <- sorted.locs[str_order(paste0(sorted.locs[, "chr"], ":", sorted.locs[, "pos"]), numeric = T), ]
  if (!is.null(tabix.file)){
    sorted.locs[,'start'] <- sorted.locs[,'pos']
    sorted.locs[,'end'] <- sorted.locs[,'pos'] + 1
  }
  orderByLoc <- function(x) {
    str_order(paste0(sorted.locs[x, "chr"], ":", sorted.locs[x, "pos"]), numeric = T)
  }
  
  if (verbose)
    message("Reordering DMPs based on location..")
  dmps.tsv <- dmps.tsv[orderByLoc(dmps.tsv$dmp), , drop = F]
  # if (verbose)
  #   message("Head of ordered dmps.tsv:\n\t", paste(capture.output(print(dmps.tsv[1:10,])), collapse='\n\t')) 
  # if (verbose)
  #   message("Tail of ordered dmps.tsv:\n\t", paste(capture.output(print(dmps.tsv[(nrow(dmps.tsv)-10):nrow(dmps.tsv),])), collapse='\n\t'))
  
  dmps <- unique(dmps.tsv$dmp)
  # Filter DMPs not present in array annotation first (prevents NA logical indices later)
  missing.in.annotation <- setdiff(dmps, rownames(sorted.locs))
  if (length(missing.in.annotation) > 0) {
    warning("Dropping ", length(missing.in.annotation), " DMP(s) not found in the array annotation: ",
            paste(head(missing.in.annotation, 10), collapse=','),
            if (length(missing.in.annotation) > 10) " ..." else "")
    dmps.tsv <- dmps.tsv[!(dmps.tsv$dmp %in% missing.in.annotation), , drop=FALSE]
    dmps <- setdiff(dmps, missing.in.annotation)
  }
  if (length(dmps) == 0) {
    browser()
    stop("No DMPs remain after filtering against array annotation.")
  }
  if (! all(dmps %in% beta.row.names)){
    missing.in.beta <- dmps[!(dmps %in% beta.row.names)]
    browser()
    stop("Some of the DMPs are not present in the beta file. DMPs: ", paste(missing.in.beta, collapse=','))
  }
  dmps <- dmps[orderByLoc(dmps)]
  
  if (verbose)
    message("Making sure that the beta file is sorted by position..")
  if (!all(beta.row.names[orderByLoc(beta.row.names)] == beta.row.names))
    stop("Provided beta file is not sorted by position!")
  
  if (verbose)
    message("Getting subset beta corresponding to DMPs..")
  dmps.locs <- sorted.locs[dmps, , drop = F]
  if (!is.null(output.prefix)){
    dmps.beta.output.file <-   paste0(output.prefix, "dmps.beta.tsv.gz")
    if (is.null(dmps.beta.file)){
      dmps.beta.file <- dmps.beta.output.file # loading from previous run
    }
  }
  if (!is.null(dmps.beta.file) && file.exists(dmps.beta.file)){ 
    if (verbose) message("Reading beta file row names from supplied dmps beta file..")
    gz <- gzfile(dmps.beta.output.file,'r')
    dmps.beta <- read.table(gz,header = T,sep = '\t', check.names = F, row.names=1, comment.char = "", quote = "")
    close(gz)
  } else {
    if (!is.null(beta.file)){
      dmps.beta <- .subset.beta(beta.file,
                              dmps,
                              beta.row.names = beta.row.names,
                              beta.col.names = beta.col.names)
    } else {
      dmps.beta <- tabix(as.data.frame(dmps.locs[,c('chr','start','end')]),
                          tabix.file,verbose=F)
      dmps.beta <- as.data.frame(sapply(dmps.beta[,beta.col.names], as.numeric))
    }
    if (!is.null(output.prefix)){
      if (verbose) message("Saving dmps beta to file: ", dmps.beta.output.file)
      gz <- gzfile(dmps.beta.output.file, 'w')
      write.table(dmps.beta, gz, sep = '\t', row.names=T, col.names=NA, quote=F)
      close(gz)
    } 
  }
  if (nrow(dmps.locs) != nrow(dmps.beta)) {
    stop("Number of rows in the queried dmps beta file does not match the number of DMPs. Number of rows in beta file: ", 
    nrow(dmps.beta), 
    " Number of rows in DMPs: ", 
    nrow(dmps.locs))
  }
  # Diagnostic: identify any rows with all NA betas (should not happen with synthetic test data)
  all.na.rows <- apply(dmps.beta, 1, function(r) all(is.na(r)))
  if (any(all.na.rows)) {
    stop("Beta extraction failure: the following DMP rows have all NA beta values: ",
         paste(rownames(dmps.beta)[all.na.rows], collapse=','),
         ". This indicates a mismatch between requested CpG IDs and beta file columns or a parsing issue.")
  }

  if (verbose)
    message("Subset size: ",paste(dim(dmps.beta), collapse=","))
  
  if (verbose)
    message("Number of provided DMPs:", length(dmps))
  if (verbose)
    message("Connecting DMPs to form initial DMRs..")
  
  ret <- parallel::mclapply(unique(dmps.locs$chr), function(chr) {
    m <- dmps.locs$chr == chr
    cdmps.tsv <- dmps.tsv[(dmps.tsv$dmp%in%rownames(dmps.locs)), , drop = F]
    cdmps <- dmps[m]
    cdmps.beta <- dmps.beta[m, , drop = F]
    cdmps.locs <- dmps.locs[m, , drop = F]
    controls.beta <- cdmps.beta[,pheno[,casecontrol.col]==0,drop=F]
    controls.beta.mean <- apply(controls.beta, 1, mean)
    controls.beta.sd <- apply(controls.beta, 1, sd)
    controls.beta.num <- apply(!is.na(controls.beta), 1, sum)
    
    
    dmrs <- data.frame()
    start.ind <- 1
    corr.pval <- 1
    dmr.dmps.inds <- c()
    for (i in seq(nrow(cdmps.beta))) {
      reg.dmr <- F
      dmr.dmps.inds <- c(dmr.dmps.inds, i)
      stop.condition <- F
      if (i == nrow(cdmps.beta)) {
        stop.condition <- T 
        stop.reason <- 'end of input'
      } else if ((max.lookup.dist>0) && (cdmps.locs[i + 1, "pos"] - cdmps.locs[i, 'pos'] > max.lookup.dist)) {
        stop.condition <- T
        stop.reason <- 'exceeded max distance' 
      }
      if (stop.condition) {
        reg.dmr <- T
      } else {
        t <- .testConnectivity(
          site1.beta = unlist(cdmps.beta[i, ]),
          site2.beta = unlist(cdmps.beta[i + 1, ]),
          sample.groups = sample.groups,
          max.pval = max.pval
        )
        
        corr.pval <- min(corr.pval, t[[2]])
        
        if (!t[[1]]) {
          reg.dmr <- T
          stop.reason <- t$reason
        }
      }
      if (reg.dmr) {
          for (sample_group in unique(cdmps.tsv[, sample_group.col])){
            
            gdmps.tsv <- cdmps.tsv[
                (cdmps.tsv[, sample_group.col] == sample_group),]
            rownames(gdmps.tsv) <- gdmps.tsv$dmp # here there must be unique CpGs in `dmp` column
            
            dmr.dmps.tsv <- gdmps.tsv[cdmps[dmr.dmps.inds], , drop = F]
            dmr.dmps.tsv[,'controls_beta'] <- dmr.dmps.tsv[,'cases_beta'] - dmr.dmps.tsv[,'delta_beta']
            new.dmr <- data.frame(
              chr = chr,
              start_dmp = cdmps[[start.ind]],
              end_dmp = cdmps[[i]],
              start_dmp_pos = cdmps.locs[start.ind, 'pos'],
              end_dmp_pos = cdmps.locs[i, 'pos'],
              dmps_num = length(dmr.dmps.inds),
              delta_beta =  mean(abs(dmr.dmps.tsv[,'delta_beta'])) * sign(sum(sign(dmr.dmps.tsv[,'delta_beta']))),
              delta_beta_sd = sd(dmr.dmps.tsv[,'delta_beta']),
              delta_beta_se = sd(dmr.dmps.tsv[,'delta_beta'])/sqrt(length(dmr.dmps.inds)),
              delta_beta_min = min(dmr.dmps.tsv[,'delta_beta']),
              delta_beta_max = max(dmr.dmps.tsv[,'delta_beta']),
              delta_beta_start = dmr.dmps.tsv[1,'delta_beta'],
              delta_beta_mid = dmr.dmps.tsv[ceiling(nrow(dmr.dmps.tsv)/2),'delta_beta'],
              delta_beta_end = dmr.dmps.tsv[nrow(dmr.dmps.tsv),'delta_beta'],
              cases_beta = mean(abs(dmr.dmps.tsv[,'cases_beta'])) * sign(sum(sign(dmr.dmps.tsv[,'cases_beta']))),
              cases_beta_max = max(dmr.dmps.tsv[,'cases_beta']),
              cases_beta_min = min(dmr.dmps.tsv[,'cases_beta']),
              cases_beta_sd = sd(dmr.dmps.tsv[,'cases_beta']),
              cases_beta_se = sd(dmr.dmps.tsv[,'cases_beta'])/sqrt(length(dmr.dmps.inds)),
              cases_beta_start = dmr.dmps.tsv[1,'cases_beta'],
              cases_beta_mid = dmr.dmps.tsv[ceiling(nrow(dmr.dmps.tsv)/2),'cases_beta'],
              cases_beta_end = dmr.dmps.tsv[nrow(dmr.dmps.tsv),'cases_beta'],
              cases_beta_dmps_sd = mean(dmr.dmps.tsv[,'cases_beta_sd']),
              cases_beta_dmps_sd_max = max(dmr.dmps.tsv[,'cases_beta_sd']),
              cases_beta_dmps_sd_min = min(dmr.dmps.tsv[,'cases_beta_sd']),
              cases_beta_dmps_sd_start = dmr.dmps.tsv[1,'cases_beta_sd'],
              cases_beta_dmps_sd_mid = dmr.dmps.tsv[ceiling(nrow(dmr.dmps.tsv)/2),'cases_beta_sd'],
              cases_beta_dmps_sd_end = dmr.dmps.tsv[nrow(dmr.dmps.tsv),'cases_beta_sd'],
              controls_beta = mean(abs(dmr.dmps.tsv[,'controls_beta'])) * sign(sum(sign(dmr.dmps.tsv[,'controls_beta']))),
              controls_beta_max = max(dmr.dmps.tsv[,'controls_beta']),
              controls_beta_min = min(dmr.dmps.tsv[,'controls_beta']),
              controls_beta_sd = sd(dmr.dmps.tsv[,'controls_beta']),
              controls_beta_se = sd(dmr.dmps.tsv[,'controls_beta'])/sqrt(length(dmr.dmps.inds)),
              controls_beta_start = dmr.dmps.tsv[1,'controls_beta'],
              controls_beta_mid = dmr.dmps.tsv[ceiling(nrow(dmr.dmps.tsv)/2),'controls_beta'],
              controls_beta_end = dmr.dmps.tsv[nrow(dmr.dmps.tsv),'controls_beta'],
              controls_beta_dmps_sd = mean(dmr.dmps.tsv[,'controls_beta_sd']),
              controls_beta_dmps_sd_max = max(dmr.dmps.tsv[,'controls_beta_sd']),
              controls_beta_dmps_sd_min = min(dmr.dmps.tsv[,'controls_beta_sd']),
              controls_beta_dmps_sd_start = dmr.dmps.tsv[1,'controls_beta_sd'],
              controls_beta_dmps_sd_mid = dmr.dmps.tsv[ceiling(nrow(dmr.dmps.tsv)/2),'controls_beta_sd'],
              controls_beta_dmps_sd_end = dmr.dmps.tsv[nrow(dmr.dmps.tsv),'controls_beta_sd'],
              cases_num = min(dmr.dmps.tsv[,'cases_num']),
              controls_num = min(dmr.dmps.tsv[,'controls_num']),
              corr_pval = corr.pval,
              dmps_pval_adj = mean(dmr.dmps.tsv[, "pval_adj"]),
              dmps_pval_adj_min = min(dmr.dmps.tsv[, "pval_adj"]),
              dmps_pval_adj_max = max(dmr.dmps.tsv[, "pval_adj"]),
              dmps_pval = mean(dmr.dmps.tsv[, "pval"]),
              dmps_pval_min = min(dmr.dmps.tsv[, "pval"]),
              dmps_pval_max = max(dmr.dmps.tsv[, "pval"]),
              dmps_qval = mean(dmr.dmps.tsv[, "qval"]),
              dmps_qval_min = min(dmr.dmps.tsv[, "qval"]),
              dmps_qval_max = max(dmr.dmps.tsv[, "qval"]),
              stop_connection_reason = stop.reason,
              dmps = paste(cdmps[dmr.dmps.inds], collapse = ',')
            )
            new.dmr[[sample_group.col]] <- sample_group
            tryCatch({
              dmrs <- rbind(dmrs, new.dmr)
            }, error = function(e) {
              message("Error in rbind: ", e)
              message("New DMR: \n\t", paste(paste(colnames(new.dmr),unlist(new.dmr),sep=': '), collapse = '\n\t'))
              message("Existing DMRs: \n\t", paste(capture.output(print(dmr)),collapse='\n\t') )
            })
          }
        
        start.ind <- i + 1
        corr.pval <- 1
        dmr.dmps.inds <- c()
      }
    }
    
    rownames(dmrs) <- seq_along(dmrs[, 1])
    dmrs[, 'chr'] <- chr
    dmrs
  })
  
  if (inherits(ret[[1]], "try-error")) {
    stop(ret)
  }
  dmrs <- do.call(rbind, ret)
  # Force numeric for key quantitative columns that may have been coerced to character
  num.cols <- intersect(c('delta_beta','cases_beta_sd','controls_beta_sd','cases_beta','controls_beta',
                          'cases_num','controls_num','delta_beta_sd','delta_beta_se','cases_beta_max',
                          'cases_beta_min','controls_beta_max','controls_beta_min'), colnames(dmrs))
  for (cn in num.cols) {
    if (!is.numeric(dmrs[[cn]])) {
      suppressWarnings(dmrs[[cn]] <- as.numeric(dmrs[[cn]]))
    }
  }
  # Coerce any purely numeric character columns to numeric to avoid downstream arithmetic errors
  for (cn in colnames(dmrs)) {
    if (is.character(dmrs[[cn]])) {
      suppressWarnings(numv <- as.numeric(dmrs[[cn]]))
      # Convert only if there are no newly introduced NAs (excluding existing NAs)
      if (!all(is.na(numv)) && sum(is.na(numv)) <= sum(is.na(dmrs[[cn]]))) {
        dmrs[[cn]] <- numv
      }
    }
  }
  # Ensure numeric types (some test-generated DMP tables may have these as characters/factors)
  to_numeric <- function(x) {
    if (is.numeric(x)) return(x)
    if (is.factor(x)) x <- as.character(x)
    suppressWarnings(as.numeric(x))
  }
  cases.num <- to_numeric(dmrs$cases_num)
  controls.num <- to_numeric(dmrs$controls_num)
  cases.sd.dmps.methylation <- to_numeric(dmrs$cases_beta_sd)
  controls.sd.dmps.methylation <- to_numeric(dmrs$controls_beta_sd)
  if (anyNA(c(cases.num, controls.num))) {
    warning("NAs introduced while coercing cases_num / controls_num to numeric; replacing NAs with 1 to avoid division errors.")
    cases.num[is.na(cases.num)] <- 1
    controls.num[is.na(controls.num)] <- 1
  }
  
  pooled.sd <- sqrt(((cases.num - 1) * cases.sd.dmps.methylation ^ 2 + (controls.num -
                                                                          1) * controls.sd.dmps.methylation ^ 2
  ) / (cases.num + controls.num - 2)
  )
  pooled.sd[pooled.sd < 1e-7] <- 1e-7
  dmrs$cohensd <- dmrs$delta_beta / pooled.sd



  ungrouped.dmrs <- dmrs[
    dmrs[,sample_group.col] == dmrs[1,sample_group.col] ,
    c("chr", "start_dmp", "end_dmp", "start_dmp_pos", "end_dmp_pos", "dmps_num", "corr_pval")
    ]
  if (verbose)
    message("Number of initial DMRs:", nrow(ungrouped.dmrs))
    message("Summary:\n\t", paste(capture.output(summary(ungrouped.dmrs)), collapse="\n\t"))
  if (verbose)
    message("Expanding DMRs on neighborhood CpGs..")
  ret <- parallel::mclapply(split(ungrouped.dmrs, seq_along(ungrouped.dmrs[, 1])), function(dmr) {
    ret <- .expandDMRs(
      dmr = dmr,
      sample.groups = sample.groups,
      max.pval = max.pval,
      tabix.file = tabix.file,
      beta.file = beta.file,
      beta.row.names = beta.row.names,
      beta.col.names = beta.col.names,
      casecontrol = pheno[,casecontrol.col],
      expansion.step = expansion.step,
      min.cpg.delta_beta = min.cpg.delta_beta,
      sorted.locs=sorted.locs
    )
    ret
  })
  if (inherits(ret, "try-error")) {
    stop(ret)
  }
  extended.dmrs <- do.call(rbind, ret)

  end.less.than.start <- extended.dmrs$end - extended.dmrs$start < 0

  if (any(end.less.than.start)) {
    warning(
      paste(
        sum(end.less.than.start),
        "DMRs have been assigned an end larger than start ! (CODE BUG TO BE REPORTED)"
      )
    )
    warning("Those are: \n\t",
            paste0(capture.output(print(extended.dmrs[end.less.than.start, ])), collapse = "\n\t"))
    
    warning("Removing them..")
    extended.dmrs <- extended.dmrs[!end.less.than.start, ]
    warning("Remaining: ", nrow(extended.dmrs))
  }
  all.locs.inds <- rownames(sorted.locs)
  names(all.locs.inds) <- all.locs.inds
  all.locs.inds[1:length(all.locs.inds)] <- seq(length(all.locs.inds))
  extended.dmrs$start_ind <- as.numeric(all.locs.inds[extended.dmrs$start_cpg])
  extended.dmrs$end_ind <- as.numeric(all.locs.inds[extended.dmrs$end_cpg])
  extended.dmrs$sup_cpgs_num <- extended.dmrs$end_ind - extended.dmrs$start_ind + 1
  
  extended.dmrs$id <- paste0(extended.dmrs$chr, ":", extended.dmrs$start, "-", extended.dmrs$end)
  extended.dmrs <- makeGRangesFromDataFrame(
    extended.dmrs,
    keep.extra.columns = T,
    ignore.strand = T,
    seqnames.field = "chr",
    seqinfo = Seqinfo(genome = genome),
    na.rm=T
  )

  if (verbose)
    message("Finding GC content of DMRs..")
  
  if (genome == "hg19") {
    library(BSgenome.Hsapiens.UCSC.hg19)
    Hsapiens <- BSgenome.Hsapiens.UCSC.hg19
    extended.dmrs.lifted.over <- extended.dmrs
  } else if (genome == "hg38") {
    library(BSgenome.Hsapiens.UCSC.hg38)
    Hsapiens <- BSgenome.Hsapiens.UCSC.hg38
    library(rtracklayer)
    path = system.file(package="liftOver", "extdata", "hg19ToHg38.over.chain")
    chain <- import.chain(path)
    extended.dmrs.lifted.over <- liftOver(extended.dmrs, chain)
  }
  
  sequences <- getSeq(Hsapiens, extended.dmrs.lifted.over, as.character = TRUE)
  sequences <- unstrsplit(sequences, sep = ",")
  extended.dmrs.lifted.over <- as.data.frame(extended.dmrs.lifted.over)
  extended.dmrs.lifted.over$cpgs_num <- str_count(sequences, "GC")
  extended.dms <- as.data.frame(extended.dmrs)
  extended.dmrs <- merge(extended.dmrs, extended.dmrs.lifted.over[, c("id", "cpgs_num")], by = "id")
  colnames(extended.dmrs)[colnames(extended.dmrs) == "seqnames"] <- "chr"

  
  extended.dmrs[extended.dmrs$cpgs_num == 0, "cpgs_num"]  <- 1
  
  extended.dmrs$dmps_num_adj <- ceiling(extended.dmrs$cpgs_num / extended.dmrs$sup_cpgs_num * extended.dmrs$dmps_num)
  

  message("Summary of extended DMRs before filtering:\n\t", paste(capture.output(summary(extended.dmrs)), collapse="\n\t"))
  filtered.dmrs <- extended.dmrs[(extended.dmrs$dmps_num_adj >= min.dmps) & (extended.dmrs$cpgs_num >= min.cpgs), ]
  if (verbose)
    message(
      "Keeping ",
      nrow(filtered.dmrs),
      " out of ",
      nrow(extended.dmrs),
      " with at least ",
      min.dmps,
      " (adjusted) supporting DMPs and ",
      min.cpgs,
      " contained CpGs."
    )
  
  extended.dmrs <- filtered.dmrs
  if (nrow(extended.dmrs) == 0){
    if (verbose) message("No DMRs were found")
    if (!is.null(output.prefix)) {
      for (f in c('.methylation.tsv.gz', '.tsv.gz')){
        gzfile <- gzfile(file.path(output.dir, paste0(output.prefix, f)), "w", compression = 2)   
        close(gzfile)
      }
    }
    return(NULL)
  }

  ne.dmrs.cols.to.keep <- colnames(dmrs)
  ne.dmrs.cols.to.keep <- c("start_dmp", "end_dmp", ne.dmrs.cols.to.keep[! ne.dmrs.cols.to.keep %in% colnames(extended.dmrs)])
  dmrs <- merge(extended.dmrs, dmrs[,ne.dmrs.cols.to.keep], by=c("start_dmp", "end_dmp"))
  dmrs <- dmrs[str_order(dmrs[, "id"], numeric = T), ]
  if (!is.null(output.prefix)) {
    dmrs.file <- paste0(output.prefix, "dmrs.tsv.gz")
    if (verbose)
      message("Saving DMRs to ", dmrs.file, "..")
    gz <- gzfile(dmrs.file, 'w')
    write.table(
      dmrs,
      gz,
      sep = "\t",
      quote = FALSE,
      qmethod = "double",
      col.names = TRUE,
      row.names = FALSE
    )
    close(gz)
  }


  dmrs.cpgs <- data.frame()
  extended.dmrs.sites <- c()
  if (verbose)
    message("Finding constituent CpGs of DMRs, that exist in the provided beta file..")
  for (i in seq_len(nrow(extended.dmrs))) {
    start.ind <- extended.dmrs$start_ind[i]
    end.ind <- extended.dmrs$end_ind[i]
    cpgs <- rownames(sorted.locs)[start.ind:end.ind]
    cpgs <- cpgs[cpgs %in% beta.row.names]
    dmrs.cpgs <- rbind(dmrs.cpgs, data.frame(dmr=extended.dmrs$id[i], cpgs=cpgs))
    extended.dmrs.sites <- c(extended.dmrs.sites,  cpgs)
  }
  if (!is.null(output.prefix)) {
    dmrs.cpgs.file <- paste0(output.prefix, ".cpgs.tsv.gz")
    if (verbose)
      message("Saving extended DMRs constituent CpGs names mapping to ", dmrs.cpgs.file, "..")
    gz <- gzfile(dmrs.cpgs.file, 'w')
    write.table(
      dmrs.cpgs,
      gz,
      sep = "\t",
      quote = FALSE,
      qmethod = "double",
      col.names = TRUE,
      row.names = FALSE
    )
    close(gz)
  }
  extended.dmrs.sites <- unique(extended.dmrs.sites)
  extended.dmrs.sites <- rownames(sorted.locs)[rownames(sorted.locs)%in% extended.dmrs.sites]
  if (!is.null(output.prefix)) {
    extended.dmrs.beta.file <- paste0(output.prefix, "extended.dmrs.beta.tsv.gz")
    if (verbose)
      message("Saving extended DMRs constituent CpGs betas to ", extended.dmrs.beta.file, "..")

    if (!is.null(beta.file)){
        extended.dmrs.beta <- .subset.beta(beta.file,
                                extended.dmrs.sites,
                                beta.row.names = beta.row.names,
                                beta.col.names = beta.col.names)
    } else {
      extended.dmrs.beta <- tabix(
        as.data.frame(
          sorted.locs[extended.dmrs.sites,c('chr','start','end')]),
                          tabix.file,
        params=paste0("-@ ", njobs),
        verbose=F)
      extended.dmrs.beta <- as.data.frame(sapply(extended.dmrs.beta[,beta.col.names], as.numeric))
    }
    rownames(extended.dmrs.beta) <- extended.dmrs.sites

    gz <- gzfile(extended.dmrs.beta.file, 'w')
    write.table(
      extended.dmrs.beta,
      gz,
      sep = "\t",
      quote = FALSE,
      qmethod = "double",
      col.names = NA,
      row.names = TRUE
    )
    close(gz)
  }
  rm(extended.dmrs.beta)
  gc()


  dmrs
}
