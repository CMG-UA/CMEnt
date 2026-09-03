#' Find DMPs from methylation array beta values using limma
#'
#' This helper identifies differentially methylated positions (DMPs) from
#' methylation array beta values using limma, returning the same seed-style
#' output columns as [findDMPsBSSeq()].
#' In contrast to [dmpFinder()] from minfi, this function supports flexible covariate inclusion
#' and returns a consistent output format with the BS-seq DMPs for downstream compatibility with CMEnt's region-finding functions.
#'
#' @param beta A beta input supported by [getBetaHandler()], such as a numeric
#'   matrix/data frame or a beta file path. minfi `MethylSet`/`RatioSet` inputs
#'   are also accepted and converted to beta values.
#' @param samplesheet A data frame or file path to a tab-delimited sample sheet.
#' @param samplesheet_sep Separator for samplesheet files. Default is tab.
#' @param sample_group_col Column in `samplesheet` containing group labels.
#' @param id_col Column in `samplesheet` containing sample IDs. row.names can also be used by specifying `id_col = "row.names"`.
#' @param array Array platform passed to [getSortedGenomicLocs()] when
#'   `sorted_locs` is not supplied.
#' @param genome Optional genome passed to [getSortedGenomicLocs()] when
#'   `sorted_locs` is not supplied. If `NULL`, inferred from `beta` when
#'   possible, otherwise from `array`.
#' @param sorted_locs Optional data frame of probe locations with row names as
#'   site IDs and `chr` plus `start` or `pos` columns.
#' @param njobs Number of jobs used by [getBetaHandler()] when reading beta
#'   files.
#' @param chr Chromosomes to retain, `"auto"` for chr1-chr22, or `"all"` for
#'   chr1-chr22 plus chrX and chrY.
#' @param case_group Group label to treat as case. If `NULL`, the first group in
#'   `sample_group_col` is used.
#' @param covariates Optional covariate column names, or a comma-separated
#'   string, to include in the limma model.
#' @param output_file Optional tab-delimited output path. Files ending in `.gz`
#'   are gzipped.
#'
#' @return A data frame with columns `chr`, `start`, `end`, `site_id`, `pval`,
#'   `qval`, `delta_beta`, and `score`.
#'
#' @examples
#' beta <- matrix(
#'     c(
#'         0.10, 0.12, 0.11, 0.82, 0.80, 0.81,
#'         0.74, 0.72, 0.73, 0.20, 0.22, 0.21,
#'         0.40, 0.42, 0.41, 0.43, 0.44, 0.45
#'     ),
#'     nrow = 3,
#'     byrow = TRUE,
#'     dimnames = list(paste0("cg", 1:3), paste0("s", 1:6))
#' )
#' samplesheet <- data.frame(
#'     Sample_ID = colnames(beta),
#'     Sample_Group = rep(c("control", "case"), each = 3)
#' )
#' sorted_locs <- data.frame(
#'     chr = "chr1",
#'     start = c(100L, 200L, 300L),
#'     row.names = rownames(beta)
#' )
#' dmps <- findDMPsArray(
#'     beta = beta,
#'     samplesheet = samplesheet,
#'     sample_group_col = "Sample_Group",
#'     array = "450K",
#'     genome = "hg19",
#'     sorted_locs = sorted_locs,
#'     case_group = "case",
#'     chr = "all"
#' )
#' head(dmps)
#'
#' @importFrom stats as.formula p.adjust
#' @export
findDMPsArray <- function(
    beta,
    samplesheet,
    samplesheet_sep = "\t",
    sample_group_col = NULL,
    id_col = "Sample_ID",
    array = NULL,
    genome = NULL,
    sorted_locs = NULL,
    njobs = getOption("CMEnt.njobs", 1L),
    chr = "auto",
    case_group = NULL,
    covariates = NULL,
    output_file = NULL
) {
    .assertPackagesInstalled(
        pkg_names = "limma",
        context = "findDMPsArray()",
        reason = "Array DMP calling is delegated to limma."
    )

    if (is_bsseq(beta)) {
        stop("BSseq inputs should be analyzed with findDMPsBSSeq() so coverage counts are retained.")
    }
    if (inherits(beta, "MethylSet") || inherits(beta, "RatioSet")) {
        .assertPackagesInstalled("minfi", "findDMPsArray()", "minfi objects must be converted to beta values.")
        beta <- minfi::getBeta(beta)
    }
    sample_group_col <- .requireSampleGroupCol(sample_group_col, "findDMPsArray()")
    if (is.null(array)) {
        stop("array must be provided for array DMP calling.", call. = FALSE)
    }
    array <- match.arg(array, c("450K", "27K", "EPIC", "EPICv2", "Mouse"))
    genome <- .resolveBuildDMRsGenome(beta = beta, array = array, genome = genome)

    beta_handler <- getBetaHandler(
        beta = beta,
        array = array,
        genome = genome,
        sorted_locs = sorted_locs,
        njobs = njobs
    )
    sample_ids <- beta_handler$getBetaColNames()
    if (anyDuplicated(sample_ids) > 0L) {
        stop("Duplicate sample IDs found in beta column names.")
    }

    pheno <- .readDMPsSampleSheet(samplesheet, samplesheet_sep)
    pheno <- .prepareDMPsPheno(
        pheno = pheno,
        sample_ids = sample_ids,
        sample_group_col = sample_group_col,
        id_col = id_col,
        case_group = case_group,
        covariates = covariates
    )
    covariates <- attr(pheno, "covariates")

    chr <- .normalizeDMPsChromosomes(chr)
    beta_mat <- as.matrix(beta_handler$getBeta(chr = chr))
    if (nrow(beta_mat) == 0L) {
        stop("No loci remain after chromosome filtering.")
    }
    if (anyDuplicated(rownames(beta_mat)) > 0L) {
        stop("Duplicate site IDs found in beta row names.")
    }
    locs <- beta_handler$getBetaLocs()[rownames(beta_mat), , drop = FALSE]

    eps <- 1e-6
    m_values <- log2(pmax(pmin(beta_mat, 1 - eps), eps) / (1 - pmax(pmin(beta_mat, 1 - eps), eps)))
    condition <- pheno$condition
    design_terms <- c("condition", .quoteNonSyntactic(covariates))
    design <- stats::model.matrix(stats::as.formula(paste("~", paste(design_terms, collapse = " + "))), data = pheno)
    coef_name <- "conditioncase"
    if (!coef_name %in% colnames(design)) {
        stop("Could not construct case/control contrast from the provided samplesheet.")
    }
    design_qr <- qr(design)
    if (design_qr$rank < ncol(design)) {
        keep_cols <- sort(design_qr$pivot[seq_len(design_qr$rank)])
        collinear_covs <- setdiff(colnames(design), colnames(design)[keep_cols])
        if (coef_name %in% collinear_covs) {
            stop("The case/control contrast is collinear with the design and cannot be estimated.")
        }
        warning("The following design columns are collinear and will be removed: ", paste(collinear_covs, collapse = ", "))
        design <- design[, keep_cols, drop = FALSE]
    }

    fit <- limma::eBayes(limma::lmFit(m_values, design))
    tab <- limma::topTable(fit, coef = coef_name, number = Inf, sort.by = "none")

    case_idx <- condition == "case"
    control_idx <- condition == "control"
    delta_beta <- rowMeans(beta_mat[, case_idx, drop = FALSE], na.rm = TRUE) -
        rowMeans(beta_mat[, control_idx, drop = FALSE], na.rm = TRUE)
    tab$delta_beta <- delta_beta
    tab$site_id <- rownames(beta_mat)
    tab$chr <- locs$chr
    tab$start <- as.integer(locs$start)

    tab <- tab[
        is.finite(tab$P.Value) &
            is.finite(tab$delta_beta),
        ,
        drop = FALSE
    ]
    tab$qval <- stats::p.adjust(tab$P.Value, method = "BH")

    if (nrow(tab) == 0L) {
        stop("limma returned no usable sites for DMP generation.")
    }
    tab$score <- -log10(tab$P.Value + .Machine$double.xmin) * abs(tab$delta_beta)
    tab$score <- (tab$score - min(tab$score)) / (max(tab$score) - min(tab$score) + .Machine$double.xmin)

    dmps <- data.frame(
        chr = as.character(tab$chr),
        start = as.integer(tab$start),
        end = as.integer(tab$start),
        site_id = tab$site_id,
        pval = tab$P.Value,
        qval = tab$qval,
        delta_beta = tab$delta_beta,
        score = tab$score,
        stringsAsFactors = FALSE
    )
    if (nrow(dmps) > 0L) {
        dmps <- dmps[order(dmps$chr, dmps$start, na.last = TRUE), , drop = FALSE]
    }
    rownames(dmps) <- NULL

    if (!is.null(output_file)) {
        .writeDMPsTable(dmps, output_file)
    }
    dmps
}

.readDMPsSampleSheet <- function(samplesheet, samplesheet_sep) {
    if (is.character(samplesheet) && length(samplesheet) == 1L) {
        pheno <- utils::read.table(samplesheet, header = TRUE, stringsAsFactors = FALSE, sep = samplesheet_sep, check.names = FALSE)
        # If the samplesheet first column is unnamed, assume it's row names and move to row names
        if ((ncol(pheno) > 0L) && (colnames(pheno)[1L] == "")) {
            rownames(pheno) <- pheno[[1L]]
            pheno <- pheno[, -1L, drop = FALSE]
        }
    } else if (is.data.frame(samplesheet)) {
        pheno <- samplesheet
    } else {
        stop("samplesheet argument must be a file path or a data frame.")
    }
    colnames(pheno) <- trimws(colnames(pheno))
    pheno
}

.prepareDMPsPheno <- function(pheno, sample_ids, sample_group_col, id_col, case_group, covariates) {
    if (id_col == "row.names") {
        if (is.null(rownames(pheno))) {
            stop("ID column is specified as 'row.names' but the samplesheet has no row names.")
        }
        pheno[[id_col]] <- rownames(pheno)
    }
    if (!all(c(sample_group_col, id_col) %in% colnames(pheno))) {
        stop("Group column ", sample_group_col, " or ID column ", id_col, " not found in samplesheet.")
    }
    pheno[[id_col]] <- as.character(pheno[[id_col]])
    if (anyDuplicated(pheno[[id_col]]) > 0L) {
        stop("Duplicate sample IDs found in the samplesheet ID column.")
    }
    if (!all(sample_ids %in% pheno[[id_col]])) {
        missing_samples <- setdiff(sample_ids, pheno[[id_col]])
        warning("The following samples are in beta but not in the samplesheet: ", paste(missing_samples, collapse = ", "))
        stop("Not all samples in beta are present in the samplesheet.")
    }
    if (is.null(case_group)) {
        case_group <- unique(pheno[[sample_group_col]])[[1L]]
    } else if (!case_group %in% pheno[[sample_group_col]]) {
        stop("Specified case group ", case_group, " not found in group column ", sample_group_col)
    }

    covariates <- .normalizeDMPsCovariates(covariates)
    missing_covariates <- setdiff(covariates, colnames(pheno))
    if (length(missing_covariates) > 0L) {
        stop("The following covariates are missing from the provided samplesheet: ", paste(missing_covariates, collapse = ","))
    }

    pheno <- pheno[match(sample_ids, pheno[[id_col]]), , drop = FALSE]
    rownames(pheno) <- pheno[[id_col]]
    pheno$condition <- factor(
        pheno[[sample_group_col]] == case_group,
        levels = c(FALSE, TRUE),
        labels = c("control", "case")
    )
    if (length(unique(pheno$condition)) < 2L) {
        stop("Condition has fewer than two levels after applying case_group; cannot fit contrast.")
    }
    attr(pheno, "case_group") <- case_group
    attr(pheno, "covariates") <- covariates
    pheno
}

.normalizeDMPsCovariates <- function(covariates) {
    if (is.null(covariates)) {
        return(NULL)
    }
    if (length(covariates) == 1L && is.character(covariates)) {
        covariates <- unlist(strsplit(covariates, ","))
    }
    covariates <- trimws(as.character(covariates))
    covariates <- covariates[nzchar(covariates)]
    if (length(covariates) == 0L) {
        NULL
    } else {
        covariates
    }
}

.quoteNonSyntactic <- function(x) {
    if (is.null(x)) {
        return(NULL)
    }
    ifelse(make.names(x) == x, x, paste0("`", x, "`"))
}

.normalizeDMPsChromosomes <- function(chr) {
    if (length(chr) == 1L && identical(chr, "auto")) {
        paste0("chr", seq_len(22L))
    } else if (length(chr) == 1L && identical(chr, "all")) {
        c(paste0("chr", seq_len(22L)), "chrX", "chrY")
    } else {
        trimws(unlist(strsplit(chr, ",")))
    }
}
