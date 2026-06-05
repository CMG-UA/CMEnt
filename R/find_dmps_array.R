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
#' @param group_col Column in `samplesheet` containing group labels.
#' @param id_col Column in `samplesheet` containing sample IDs.
#' @param array Array platform passed to [getSortedGenomicLocs()] when
#'   `sorted_locs` is not supplied.
#' @param genome Genome passed to [getSortedGenomicLocs()] when `sorted_locs` is
#'   not supplied.
#' @param sorted_locs Optional data frame of probe locations with row names as
#'   site IDs and `chr` plus `start` or `pos` columns.
#' @param njobs Number of jobs used by [getBetaHandler()] when reading beta
#'   files.
#' @param chr Chromosomes to retain, `"auto"` for chr1-chr22, or `"all"` for
#'   chr1-chr22 plus chrX and chrY.
#' @param case_group Group label to treat as case. If `NULL`, the first group in
#'   `group_col` is used.
#' @param covariates Optional covariate column names, or a comma-separated
#'   string, to include in the limma model.
#' @param fdr_thres FDR threshold for returned DMPs. Default is 0.05.
#' @param output_file Optional tab-delimited output path. Files ending in `.gz`
#'   are gzipped.
#'
#' @return A data frame with columns `chr`, `start`, `end`, `site_id`, `pval`,
#'   `qval`, `delta_beta`, and `score`.
#'
#' @importFrom stats as.formula p.adjust
#' @export
findDMPsArray <- function(
    beta,
    samplesheet,
    samplesheet_sep = "\t",
    group_col = "Sample_Group",
    id_col = "Sample_ID",
    array = c("450K", "27K", "EPIC", "EPICv2", "Mouse"),
    genome = c("hg19", "hg38", "hs1", "mm10", "mm39"),
    sorted_locs = NULL,
    njobs = getOption("CMEnt.njobs", 1L),
    chr = "auto",
    case_group = NULL,
    covariates = NULL,
    fdr_thres = 0.05,
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
        group_col = group_col,
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
    if (qr(design)$rank < ncol(design)) {
        stop("The array DMP model is rank deficient; remove collinear covariates.")
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
    tab$qval <- stats::p.adjust(tab$P.Value, method = "BH")

    tab <- tab[
        is.finite(tab$P.Value) &
            is.finite(tab$qval) &
            is.finite(tab$delta_beta),
        ,
        drop = FALSE
    ]
    if (nrow(tab) == 0L) {
        stop("limma returned no usable sites for DMP generation.")
    }
    tab$score <- -log10(tab$P.Value + .Machine$double.xmin) * abs(tab$delta_beta)
    tab$score <- (tab$score - min(tab$score)) / (max(tab$score) - min(tab$score) + .Machine$double.xmin)

    tab <- tab[tab$qval <= fdr_thres, , drop = FALSE]
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
    } else if (is.data.frame(samplesheet)) {
        pheno <- samplesheet
    } else {
        stop("samplesheet argument must be a file path or a data frame.")
    }
    colnames(pheno) <- trimws(colnames(pheno))
    pheno
}

.prepareDMPsPheno <- function(pheno, sample_ids, group_col, id_col, case_group, covariates) {
    if (!all(c(group_col, id_col) %in% colnames(pheno))) {
        stop("Group column ", group_col, " or ID column ", id_col, " not found in samplesheet.")
    }
    pheno[[id_col]] <- as.character(pheno[[id_col]])
    if (anyDuplicated(pheno[[id_col]]) > 0L) {
        stop("Duplicate sample IDs found in the samplesheet ID column.")
    }
    if (!all(sample_ids %in% pheno[[id_col]])) {
        stop("Not all samples in beta are present in the samplesheet.")
    }
    if (is.null(case_group)) {
        case_group <- unique(pheno[[group_col]])[[1L]]
    } else if (!case_group %in% pheno[[group_col]]) {
        stop("Specified case group ", case_group, " not found in group column ", group_col)
    }

    covariates <- .normalizeDMPsCovariates(covariates)
    missing_covariates <- setdiff(covariates, colnames(pheno))
    if (length(missing_covariates) > 0L) {
        stop("The following covariates are missing from the provided samplesheet: ", paste(missing_covariates, collapse = ","))
    }

    pheno <- pheno[match(sample_ids, pheno[[id_col]]), , drop = FALSE]
    rownames(pheno) <- pheno[[id_col]]
    pheno$condition <- factor(
        pheno[[group_col]] == case_group,
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
