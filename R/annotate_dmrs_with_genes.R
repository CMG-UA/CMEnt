.trimGeneBodiesByPromoters <- function(genes, promoters) {
    if (length(genes) == 0L || length(promoters) == 0L) {
        return(genes)
    }
    genes_unstranded <- genes
    promoters_unstranded <- promoters
    GenomicRanges::strand(genes_unstranded) <- "*"
    GenomicRanges::strand(promoters_unstranded) <- "*"
    promoter_free <- GenomicRanges::setdiff(genes_unstranded, promoters_unstranded, ignore.strand = TRUE)
    if (length(promoter_free) == 0L) {
        return(genes[FALSE])
    }
    hits <- GenomicRanges::findOverlaps(promoter_free, genes, ignore.strand = TRUE)
    if (length(hits) == 0L) {
        return(genes)
    }
    gene_hits <- S4Vectors::subjectHits(hits)
    gene_bodies <- GenomicRanges::pintersect(
        genes[gene_hits],
        promoter_free[S4Vectors::queryHits(hits)],
        ignore.strand = TRUE
    )
    names(gene_bodies) <- names(genes)[gene_hits]
    S4Vectors::mcols(gene_bodies) <- S4Vectors::mcols(genes)[gene_hits, , drop = FALSE]
    gene_bodies
}

.loadGeneAnnotationFeatures <- function(genome, promoter_upstream = 2000,
                                        promoter_downstream = 200,
                                        context = "gene annotation") {
    cache_dir <- getOption(
        "CMEnt.annotation_cache_dir",
        .getOSCacheDir(file.path("R", "CMEnt", "annotation_cache"))
    )
    target_genome <- tolower(genome)
    annotation_source_genome <- if (target_genome == "hs1") "hg38" else target_genome
    annotation_pkgs <- .assertGeneAnnotPkgsInstalled(
        genome = genome,
        context = context
    )
    txdb_pkg <- unname(annotation_pkgs[["txdb"]])
    orgdb_pkg <- unname(annotation_pkgs[["orgdb"]])
    if (annotation_source_genome != target_genome) {
        .log_info(
            "Using ", annotation_source_genome,
            " gene models lifted to ", target_genome, " for annotation.",
            level = 2
        )
    }
    .log_step("Loading gene annotations for ", genome, "...", level = 3)

    # Load TxDb namespace
    if (!isNamespaceLoaded(txdb_pkg)) {
        loadNamespace(txdb_pkg)
    }

    # Load TxDb - the main object has the same name as the package
    txdb <- getExportedValue(txdb_pkg, txdb_pkg)

    # Get genes and promoters
    # Load them from cache if available
    suppressMessages({
        genes_key <- paste0("genes_", genome)
        genes <- if (getOption("CMEnt.use_annotation_cache", TRUE)) {
            .readBiocFileCacheRDS(cache_dir, genes_key)
        } else {
            NULL
        }
        if (!is.null(genes)) {
            .log_info("Loading cached genes from ", genes_key, level = 2)
        } else {
            genes <- GenomicFeatures::genes(txdb)
            # get genes only within standard chromosomes
            std_chroms <- GenomeInfoDb::standardChromosomes(GenomeInfoDb::seqinfo(txdb))
            genes <- genes[as.character(GenomeInfoDb::seqnames(genes)) %in% std_chroms]
            if (annotation_source_genome != target_genome) {
                genes <- .liftOverFromGenomeToGenome(genes, annotation_source_genome, target_genome)
                target_std_chroms <- GenomeInfoDb::standardChromosomes(GenomeInfoDb::Seqinfo(genome = target_genome))
                genes <- genes[as.character(GenomeInfoDb::seqnames(genes)) %in% target_std_chroms]
            }
            tryCatch(
                {
                    .saveBiocFileCacheRDS(genes, cache_dir, genes_key)
                },
                warning = function(w) {
                    .log_warn(
                        "Could not write annotation cache entry '", genes_key,
                        "' (warning: ", conditionMessage(w), "). Continuing without cache persistence."
                    )
                },
                error = function(e) {
                    .log_warn(
                        "Could not write annotation cache entry '", genes_key,
                        "' (error: ", conditionMessage(e), "). Continuing without cache persistence."
                    )
                }
            )
        }
        promoters_key <- paste0("promoters_", genome, "_u", promoter_upstream, "_d", promoter_downstream)
        promoters <- if (getOption("CMEnt.use_annotation_cache", TRUE)) {
            .readBiocFileCacheRDS(cache_dir, promoters_key)
        } else {
            NULL
        }
        if (!is.null(promoters)) {
            .log_info("Loading cached promoters from ", promoters_key, level = 2)
        } else {
            transcripts_by_gene <- GenomicFeatures::transcriptsBy(txdb, by = "gene")
            transcripts_by_gene <- transcripts_by_gene[names(transcripts_by_gene) %in% names(genes)]
            promoters <- GenomicFeatures::promoters(transcripts_by_gene, upstream = promoter_upstream, downstream = promoter_downstream)
            promoters <- stack(promoters)
            if (annotation_source_genome != target_genome) {
                promoters <- .liftOverFromGenomeToGenome(promoters, annotation_source_genome, target_genome)
                target_std_chroms <- GenomeInfoDb::standardChromosomes(GenomeInfoDb::Seqinfo(genome = target_genome))
                promoters <- promoters[as.character(GenomeInfoDb::seqnames(promoters)) %in% target_std_chroms]
            }
            tryCatch(
                {
                    .saveBiocFileCacheRDS(promoters, cache_dir, promoters_key)
                },
                warning = function(w) {
                    .log_warn(
                        "Could not write annotation cache entry '", promoters_key,
                        "' (warning: ", conditionMessage(w), "). Continuing without cache persistence."
                    )
                },
                error = function(e) {
                    .log_warn(
                        "Could not write annotation cache entry '", promoters_key,
                        "' (error: ", conditionMessage(e), "). Continuing without cache persistence."
                    )
                }
            )
        }
    })
    genes <- .trimGeneBodiesByPromoters(genes, promoters)
    .log_success("Gene annotations loaded: ", length(genes), " genes", level = 3)
    list(genes = genes, promoters = promoters, orgdb_pkg = orgdb_pkg)
}

#' Annotate DMRs with Gene Information
#'
#' @description Annotates DMRs with overlapping gene promoters and gene bodies
#' using TxDb annotations. For each DMR, identifies genes whose promoters or
#' gene bodies overlap with the DMR coordinates.
#'
#' @param dmrs Dataframe or GRanges object containing DMR coordinates
#' @param genome Character. Genome version to use for gene annotation. Required for data.frame DMRs; inferred from GRanges seqinfo when available.
#' @param promoter_upstream Integer. Number of base pairs upstream of TSS to
#'   define promoter region (default: 2000)
#' @param promoter_downstream Integer. Number of base pairs downstream of TSS
#'   to define promoter region (default: 200)
#' @param njobs Integer. Number of parallel jobs used to annotate promoter and
#'   gene-body overlaps (default: `getOption("CMEnt.njobs")`)
#' @param site_locs Optional data frame or GRanges with site coordinates used to
#'   compute feature-specific delta beta values.
#' @param site_delta_beta Optional named numeric vector of per-site delta beta values.
#' @param aggfun Function used to aggregate per-site delta beta values.
#'
#' @return The input Dataframe/GRanges object with additional metadata columns:
#' \itemize{
#'   \item in_promoter_of: Character vector of gene symbols with promoters overlapping the DMR (comma-separated)
#'   \item in_gene_body_of: Character vector of gene symbols with gene bodies overlapping the DMR (comma-separated)
#'   \item delta_beta_promoter: Aggregated delta beta of DMR sites overlapping promoters, or NA
#'   \item delta_beta_gene_body: Aggregated delta beta of DMR sites overlapping gene bodies, or NA
#' }
#'
#' @details
#' The function uses genome-appropriate TxDb packages. For `hs1`, CMEnt
#' uses hg38 gene models and lifts them to `hs1` before computing overlaps.
#' Gene symbols are retrieved from the appropriate org.*.eg.db package.
#' Multiple overlapping genes are concatenated with commas.
#'
#' @examples
#' # Annotate DMRs with gene information
#' dmrs <- data.frame(
#'     chr = c("chr1", "chr2"),
#'     start = c(1000000, 2000000),
#'     end = c(1001000, 2001000)
#' )
#' dmrs_annotated <- annotateDMRsWithGenes(dmrs, genome = "hg38")
#'
#' # Use custom promoter definition
#' dmrs_annotated <- annotateDMRsWithGenes(
#'     dmrs,
#'     genome = "hg38",
#'     promoter_upstream = 5000,
#'     promoter_downstream = 1000,
#'     njobs = 2
#' )
#'
#' @export
annotateDMRsWithGenes <- function(dmrs, genome = NULL,
                                  promoter_upstream = 2000,
                                  promoter_downstream = 200,
                                  njobs = getOption("CMEnt.njobs", .defaultNJobs()),
                                  site_locs = NULL,
                                  site_delta_beta = NULL,
                                  aggfun = stats::median) {
    dmrs_df_provided <- is.data.frame(dmrs)
    dmrs <- .convertToGRanges(dmrs, genome)
    genome <- .resolveGRangesGenome(dmrs, "DMRs")

    annotation_features <- .loadGeneAnnotationFeatures(
        genome = genome,
        promoter_upstream = promoter_upstream,
        promoter_downstream = promoter_downstream,
        context = "annotateDMRsWithGenes()"
    )
    genes <- annotation_features$genes
    promoters <- annotation_features$promoters
    orgdb_pkg <- annotation_features$orgdb_pkg

    .log_step("Finding overlaps with promoters and gene bodies...", level = 3)
    annotation_specs <- list(
        list(
            column = "in_promoter_of",
            delta_column = "delta_beta_promoter",
            features = promoters,
            feature_type = "promoter"
        ),
        list(
            column = "in_gene_body_of",
            delta_column = "delta_beta_gene_body",
            features = genes,
            feature_type = "gene_body"
        )
    )
    task <- function(spec) {
        list(
            column = spec$column,
            values = .annotateDMRsWithGeneFeature(
                dmrs = dmrs,
                features = spec$features,
                orgdb_pkg = orgdb_pkg,
                feature_type = spec$feature_type
            )
        )
    }
    if (njobs > 1) {
        bp_param <- .makeBiocParallelParam(
            njobs,
            n_tasks = length(annotation_specs),
            progressbar = FALSE,
            log = getOption("CMEnt.verbose", 1L) >= 1L
        )
        annotation_results <- .safeBiocParallelApply(
            annotation_specs,
            task,
            BPPARAM = bp_param
        )
    } else {
        annotation_results <- lapply(
            annotation_specs,
            task
        )
    }
    for (annotation_result in annotation_results) {
        S4Vectors::mcols(dmrs)[[annotation_result$column]] <- annotation_result$values
    }
    delta_beta_annotations <- .annotateDMRSiteDBByFeature(
        dmrs = dmrs,
        annotation_specs = annotation_specs,
        site_locs = site_locs,
        site_delta_beta = site_delta_beta,
        aggfun = aggfun,
        genome = genome
    )
    for (delta_column in names(delta_beta_annotations)) {
        S4Vectors::mcols(dmrs)[[delta_column]] <- delta_beta_annotations[[delta_column]]
    }
    .log_success("Gene symbols mapped", level = 3)
    if (dmrs_df_provided) {
        dmrs <- .convertToDataFrame(dmrs)
    }
    dmrs
}

#' @keywords internal
#' @noRd
.annotateDMRSiteDBByFeature <- function(
    dmrs, annotation_specs, site_locs,
    site_delta_beta, aggfun, genome
) {
    out <- stats::setNames(
        replicate(length(annotation_specs), rep(NA_real_, length(dmrs)), simplify = FALSE),
        vapply(annotation_specs, `[[`, character(1), "delta_column")
    )
    if (length(dmrs) == 0L || is.null(site_locs) || is.null(site_delta_beta)) {
        return(out)
    }

    site_names <- if (is.data.frame(site_locs)) {
        rownames(site_locs)
    } else {
        names(site_locs)
    }
    site_locs <- .convertSitesToGPos(site_locs, genome)
    if (
        !is.null(site_names) &&
            length(site_names) == length(site_locs) &&
            all(!is.na(site_names)) &&
            all(nzchar(site_names))
    ) {
        names(site_locs) <- site_names
    }
    if (
        is.null(names(site_locs)) ||
            anyNA(names(site_locs)) ||
            any(!nzchar(names(site_locs)))
    ) {
        return(out)
    }

    site_delta_beta <- as.numeric(site_delta_beta)
    if (is.null(names(site_delta_beta))) {
        names(site_delta_beta) <- names(site_locs)[seq_along(site_delta_beta)]
    }
    keep <- names(site_locs) %in% names(site_delta_beta)
    site_locs <- site_locs[keep]
    if (length(site_locs) == 0L) {
        return(out)
    }
    site_names <- names(site_locs)
    site_delta_beta <- site_delta_beta[site_names]

    mcols_names <- colnames(S4Vectors::mcols(dmrs))
    if ("sites" %in% mcols_names) {
        dmr_site_idx <- lapply(
            as.character(S4Vectors::mcols(dmrs)$sites),
            function(ids) {
                idx <- match(.splitCsvValues(ids), site_names)
                unique(idx[!is.na(idx)])
            }
        )
    } else {
        site_hits <- GenomicRanges::findOverlaps(dmrs, site_locs, ignore.strand = TRUE)
        dmr_site_idx <- vector("list", length(dmrs))
        if (length(site_hits) > 0L) {
            hits_by_dmr <- split(
                S4Vectors::subjectHits(site_hits),
                S4Vectors::queryHits(site_hits)
            )
            dmr_site_idx[as.integer(names(hits_by_dmr))] <- lapply(
                hits_by_dmr,
                unique
            )
        }
    }

    aggregate_delta <- function(idx) {
        if (length(idx) == 0L) {
            return(NA_real_)
        }
        vals <- site_delta_beta[idx]
        vals <- vals[is.finite(vals)]
        if (length(vals) == 0L) {
            return(NA_real_)
        }
        val <- tryCatch(
            aggfun(vals, na.rm = TRUE),
            error = function(e) aggfun(vals)
        )
        if (length(val) == 0L || !is.finite(val[1L])) {
            return(NA_real_)
        }
        as.numeric(val[1L])
    }

    for (spec in annotation_specs) {
        feature_hits <- GenomicRanges::findOverlaps(
            site_locs,
            spec$features,
            ignore.strand = TRUE
        )
        if (length(feature_hits) == 0L) {
            next
        }
        feature_site_idx <- unique(S4Vectors::queryHits(feature_hits))
        feature_site_mask <- rep(FALSE, length(site_locs))
        feature_site_mask[feature_site_idx] <- TRUE
        out[[spec$delta_column]] <- vapply(
            dmr_site_idx,
            function(idx) aggregate_delta(idx[feature_site_mask[idx]]),
            numeric(1)
        )
    }
    out
}
