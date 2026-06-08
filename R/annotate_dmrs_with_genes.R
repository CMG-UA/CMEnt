#' Annotate DMRs with Gene Information
#'
#' @description Annotates DMRs with overlapping gene promoters and gene bodies
#' using TxDb annotations. For each DMR, identifies genes whose promoters or
#' gene bodies overlap with the DMR coordinates.
#'
#' @param dmrs Dataframe or GRanges object containing DMR coordinates
#' @param genome Character. Genome version to use for gene annotation. (default: "hg38")
#' @param promoter_upstream Integer. Number of base pairs upstream of TSS to
#'   define promoter region (default: 2000)
#' @param promoter_downstream Integer. Number of base pairs downstream of TSS
#'   to define promoter region (default: 200)
#' @param njobs Integer. Number of parallel jobs used to annotate promoter and
#'   gene-body overlaps (default: `getOption("CMEnt.njobs")`)
#'
#' @return The input Dataframe/GRanges object with additional metadata columns:
#' \itemize{
#'   \item in_promoter_of: Character vector of gene symbols with promoters overlapping the DMR (comma-separated)
#'   \item in_gene_body_of: Character vector of gene symbols with gene bodies overlapping the DMR (comma-separated)
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
annotateDMRsWithGenes <- function(dmrs, genome = "hg38",
                                  promoter_upstream = 2000,
                                  promoter_downstream = 200,
                                  njobs = getOption("CMEnt.njobs", .defaultNJobs())) {
    cache_dir <- getOption(
        "CMEnt.annotation_cache_dir",
        .getOSCacheDir(file.path("R", "CMEnt", "annotation_cache"))
    )
    dmrs_df_provided <- is.data.frame(dmrs)
    dmrs <- .convertToGRanges(dmrs, genome)

    target_genome <- tolower(genome)
    annotation_source_genome <- if (target_genome == "hs1") "hg38" else target_genome
    annotation_pkgs <- .assertGeneAnnotationPackagesInstalled(
        genome = genome,
        context = "annotateDMRsWithGenes()"
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
    .log_step("Loading gene annotations for ", genome, "...", level = 2)

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
        promoters_key <- paste0("promoters_", genome)
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
    .log_success("Gene annotations loaded: ", length(genes), " genes", level = 2)
    .log_step("Finding overlaps with promoters and gene bodies...", level = 2)
    .log_step("Mapping overlapping Entrez IDs to gene symbols...", level = 2)
    annotation_specs <- list(
        list(
            column = "in_promoter_of",
            features = promoters,
            feature_type = "promoter"
        ),
        list(
            column = "in_gene_body_of",
            features = genes,
            feature_type = "gene_body"
        )
    )
    bp_param <- .makeBiocParallelParam(
        njobs,
        n_tasks = length(annotation_specs),
        progressbar = getOption("CMEnt.verbose", 0) > 0L
    )
    annotation_results <- BiocParallel::bplapply(
        annotation_specs,
        function(spec) {
            list(
                column = spec$column,
                values = .annotateDMRsWithGeneFeature(
                    dmrs = dmrs,
                    features = spec$features,
                    orgdb_pkg = orgdb_pkg,
                    feature_type = spec$feature_type
                )
            )
        },
        BPPARAM = bp_param
    )
    for (annotation_result in annotation_results) {
        S4Vectors::mcols(dmrs)[[annotation_result$column]] <- annotation_result$values
    }
    .log_success("Gene symbols mapped", level = 2)
    if (dmrs_df_provided) {
        dmrs <- .convertToDataFrame(dmrs)
    }
    dmrs
}
