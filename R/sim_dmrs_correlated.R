# nolint start: object_name_linter, line_length_linter

#' Simulate DMRs with Correlated Effects.
#'
#' Add simulated DMRs to observed control data. Control data will be split into two (artificial) populations.
#' Extends \link[dmrseq]{simDMRs} function, with additional options to control the correlation structure of methylation changes across CpGs within a DMR
#' and across samples.
#'
#' @param bs A BSseq object containing the control data to which DMRs will be added. All loci must have coverage in every sample.
#' @param num.dmrs The number of DMRs to simulate.
#' @param delta.max0 The maximum absolute change in methylation proportion to apply within a DMR. The actual change for each DMR will be randomly varied around this value.
#' @param use.correlated.effects If TRUE, adds correlated latent deviations across CpGs within
#'   each affected DMR on top of the smooth mean shift, creating more realistic regional
#'   patterns of differential methylation. If FALSE, only the smooth deterministic DMR effect is
#'   applied, without the additional correlated deviations.
#' @param corr.mode Character. `"background"` (default) estimates a neighboring-CpG correlation
#'   distribution from the background signal in `bs` and calibrates DMR-specific correlation
#'   settings from it. `"manual"` uses the supplied `corr.sd` and `sample.sd.frac` directly.
#' @param corr.rate Numeric scalar multiplier applied to the sampled background correlation
#'   targets when `corr.mode = "background"`. Values above 1 strengthen the expected
#'   neighboring-CpG correlation signal, while values between 0 and 1 weaken it.
#' @param corr.sd Controls the magnitude of the correlated effect across CpGs within a DMR when
#'   `corr.mode = "manual"`. Higher values lead to stronger correlation and more coherent DMRs,
#'   while lower values result in weaker correlation and more heterogeneous DMRs.
#' @param corr.length Controls the length scale of the correlation across CpGs within a DMR.
#'   Smaller values lead to more localized correlation (i.e., stronger correlation between
#'   nearby CpGs and weaker correlation between distant CpGs), while larger values lead to
#'   more extended correlation across the DMR.
#' @param sample.sd.frac Controls the magnitude of sample-specific variation relative to the
#'   shared correlated effect when `corr.mode = "manual"`. Higher values make CpG-to-CpG
#'   correlation easier to observe across samples, while lower values keep replicates more
#'   similar to the shared DMR profile.
#' @return A named list containing:
#'   \item{gr.dmrs}{A GRanges object with the locations of the simulated DMRs and metadata
#'   columns `corr_target`, `corr_sd_used`, `sample_sd_frac_used`, and `corr_mode_used`.}
#'   \item{dmr.mncov}{A numeric vector with the mean coverage of each simulated DMR.}
#'   \item{dmr.L}{A numeric vector with the number of CpGs in each simulated DMR.}
#'   \item{bs}{A BSseq object containing the simulated methylation data with the added DMRs.}
#'   \item{delta}{A numeric vector with the maximum methylation change applied to each DMR.}
#' @export
#' @examples
#' \dontrun{
#' # Load example BSseq data
#' data("BS.cancer.ex", package = "bsseqData")
#' keep <- rowSums(as.matrix(bsseq::getCoverage(BS.cancer.ex, type = "Cov")) == 0) == 0
#' # Simulate DMRs with background-calibrated correlation
#' sim_data <- simDMRsCorrelated(BS.cancer.ex[keep, ],
#'     num.dmrs = 100,
#'     delta.max0 = 0.2,
#'     use.correlated.effects = TRUE,
#'     corr.rate = 1.2,
#'     corr.length = 200
#' )
#'
#' # Manual correlation settings remain available
#' sim_data_manual <- simDMRsCorrelated(BS.cancer.ex[keep, ],
#'     num.dmrs = 100,
#'     delta.max0 = 0.2,
#'     corr.mode = "manual",
#'     corr.sd = 0.3,
#'     corr.length = 200,
#'     sample.sd.frac = 0.75
#' )
#' }
simDMRsCorrelated <- function(bs,
                              num.dmrs = 3000,
                              delta.max0 = 0.3,
                              use.correlated.effects = TRUE,
                              corr.mode = c("background", "manual"),
                              corr.rate = 1,
                              corr.sd = 0.30,
                              corr.length = 200,
                              sample.sd.frac = 0.75) {
    default.corr.sd <- 0.30
    default.sample.sd.frac <- 0.75
    corr.mode.missing <- missing(corr.mode)
    nondefault.manual <- (
        !missing(corr.sd) &&
            !isTRUE(all.equal(corr.sd, default.corr.sd))
    ) || (
        !missing(sample.sd.frac) &&
            !isTRUE(all.equal(sample.sd.frac, default.sample.sd.frac))
    )
    corr.mode <- match.arg(corr.mode)
    if (corr.mode.missing && nondefault.manual) {
        corr.mode <- "manual"
        message(
            "Detected non-default manual correlation settings without 'corr.mode'; ",
            "using corr.mode = 'manual' for backward compatibility."
        )
    } else if (!corr.mode.missing && corr.mode == "background" && nondefault.manual) {
        warning(
            "Ignoring non-default 'corr.sd'/'sample.sd.frac' because ",
            "corr.mode = 'background'. Use corr.mode = 'manual' to apply them."
        )
    }
    if (!is.numeric(corr.rate) || length(corr.rate) != 1L ||
        !is.finite(corr.rate) || corr.rate <= 0) {
        stop("'corr.rate' must be a single positive numeric value.")
    }
    if (!is.numeric(corr.sd) || length(corr.sd) != 1L ||
        !is.finite(corr.sd) || corr.sd < 0) {
        stop("'corr.sd' must be a single non-negative numeric value.")
    }
    if (!is.numeric(sample.sd.frac) || length(sample.sd.frac) != 1L ||
        !is.finite(sample.sd.frac) || sample.sd.frac < 0) {
        stop("'sample.sd.frac' must be a single non-negative numeric value.")
    }
    if (!is.numeric(corr.length) || length(corr.length) != 1L ||
        !is.finite(corr.length) || corr.length <= 0) {
        stop("'corr.length' must be a single positive numeric value.")
    }

    # check that all loci have coverage in every sample
    zero.cov <- which(rowSums(as.matrix(bsseq::getCoverage(bs, type = "Cov")) == 0) > 0)
    if (length(zero.cov) > 0) {
        stop(
            "Zero coverage found for at least one sample in ", length(zero.cov),
            " loci. Please filter for loci with coverage at least one in ",
            "all samples before passing to 'simDMRs_correlated'"
        )
    }

    sampleSize <- floor(ncol(bs) / 2)
    message(
        "Simulating DMRs for ", sampleSize, " vs ", ncol(bs) - sampleSize,
        " comparison"
    )

    compute_pair_correlations <- function(m, start.ind, end.ind) {
        if (length(start.ind) == 0L || ncol(m) < 2L) {
            return(numeric(0))
        }
        x <- m[start.ind, , drop = FALSE]
        y <- m[end.ind, , drop = FALSE]
        x.centered <- x - rowMeans(x)
        y.centered <- y - rowMeans(y)
        denom <- sqrt(rowSums(x.centered^2) * rowSums(y.centered^2))
        cors <- rowSums(x.centered * y.centered) / denom
        cors[is.finite(cors)]
    }

    triwt <- function(x, amp = 1, base = 0, width = 1, center = 0,
                      deg = 3, dir = 1) {
        y <- dir * (((width / 2)^deg - abs(x - center)^deg) /
            (width / 2)^deg)^deg * amp + base
        y[abs(x - center) > ceiling(width / 2)] <-
            base[abs(x - center) > ceiling(width / 2)]
        y
    }

    # correlated latent field over CpGs inside a DMR
    # It is generated by sampling from a multivariate normal distribution with an exponential covariance structure,
    # where the correlation between two CpGs decays exponentially with their distance.
    # The 'sigma' parameter controls the overall magnitude of the correlated effect,
    # and the 'length.scale' parameter controls how quickly the correlation decays with distance
    # (smaller values lead to more localized correlation).
    make_correlated_effect <- function(pos, sigma = 0.25, length.scale = 200) {
        n <- length(pos)
        if (n == 1L || sigma <= 0) {
            return(rep(0, n))
        }
        d <- abs(outer(pos, pos, "-"))
        R <- exp(-d / length.scale)
        diag(R) <- diag(R) + 1e-8
        z <- t(chol(R)) %*% rnorm(n)
        as.numeric(sigma * z)
    }

    # apply deterministic smooth effect plus correlated deviation on logit scale
    make_shifted_prob <- function(base.p, diff.hit, pos,
                                  group.field = NULL,
                                  sample.sd = 0,
                                  length.scale = 200,
                                  use.cor = TRUE) {
        eps <- 1e-6
        base.p <- pmin(pmax(base.p, eps), 1 - eps)

        # preserve the original interpretation of Diff.hit as an additive change
        # on the probability scale, then map that change to a logit shift
        target.p <- pmin(pmax(base.p + diff.hit, eps), 1 - eps)
        det.shift <- qlogis(target.p) - qlogis(base.p)

        if (is.null(group.field)) {
            group.field <- rep(0, length(base.p))
        }

        sample.field <- if (use.cor && sample.sd > 0) {
            make_correlated_effect(pos,
                sigma = sample.sd,
                length.scale = length.scale
            )
        } else {
            rep(0, length(base.p))
        }

        plogis(qlogis(base.p) + det.shift + group.field + sample.field)
    }

    calibration_metric <- function(corr.sd, sample.sd.frac, affected.samples,
                                   template.pos, template.cov,
                                   length.scale, pilot.reps = 8L) {
        if (affected.samples < 2L || length(template.pos) < 2L) {
            return(0)
        }
        base.p <- matrix(0.5, nrow = length(template.pos), ncol = affected.samples)
        diff.hit <- rep(0, length(template.pos))
        rep.metrics <- numeric(pilot.reps)
        for (rep.idx in seq_len(pilot.reps)) {
            group.field <- make_correlated_effect(
                template.pos,
                sigma = corr.sd,
                length.scale = length.scale
            )
            sample.sd <- corr.sd * sample.sd.frac
            beta.sim <- matrix(NA_real_, nrow = length(template.pos), ncol = affected.samples)
            for (sample.idx in seq_len(affected.samples)) {
                new.p <- make_shifted_prob(
                    base.p = base.p[, sample.idx],
                    diff.hit = diff.hit,
                    pos = template.pos,
                    group.field = group.field,
                    sample.sd = sample.sd,
                    length.scale = length.scale,
                    use.cor = TRUE
                )
                beta.sim[, sample.idx] <- stats::rbinom(
                    n = length(template.pos),
                    size = template.cov,
                    prob = new.p
                ) / template.cov
            }
            m.sim <- .transformBeta(
                beta.sim,
                pheno = data.frame(sample = seq_len(affected.samples))
            )
            pos.cors <- compute_pair_correlations(
                m.sim,
                start.ind = seq_len(nrow(m.sim) - 1L),
                end.ind = seq_len(nrow(m.sim) - 1L) + 1L
            )
            pos.cors <- pos.cors[pos.cors > 0]
            rep.metrics[rep.idx] <- if (length(pos.cors) == 0L) 0 else mean(pos.cors)
        }
        mean(rep.metrics)
    }

    build_lookup_grid <- function(affected.sample.sizes, template.pos, template.cov,
                                  length.scale, pilot.reps = 8L) {
        if (length(affected.sample.sizes) == 0L) {
            return(list())
        }
        param.grid <- expand.grid(
            corr.sd = c(0.20, 0.25, 0.30, 0.35, 0.40),
            sample.sd.frac = c(0.50, 0.75, 1.00),
            KEEP.OUT.ATTRS = FALSE,
            stringsAsFactors = FALSE
        )
        withr::with_seed(1, {
            setNames(lapply(affected.sample.sizes, function(affected.samples) {
                grid.metrics <- vapply(seq_len(nrow(param.grid)), function(idx) {
                    calibration_metric(
                        corr.sd = param.grid$corr.sd[idx],
                        sample.sd.frac = param.grid$sample.sd.frac[idx],
                        affected.samples = affected.samples,
                        template.pos = template.pos,
                        template.cov = template.cov,
                        length.scale = length.scale,
                        pilot.reps = pilot.reps
                    )
                }, numeric(1))
                cbind(param.grid, metric = grid.metrics)
            }), as.character(affected.sample.sizes))
        })
    }

    fit_background_gamma <- function(pos.cors) {
        if (length(pos.cors) < 20L) {
            return(list(ok = FALSE))
        }
        cors.mean <- mean(pos.cors)
        cors.var <- stats::var(pos.cors)
        if (!is.finite(cors.mean) || !is.finite(cors.var) || cors.mean <= 0 ||
            cors.var <= 1e-8) {
            return(list(ok = FALSE))
        }
        shape <- cors.mean^2 / cors.var
        rate <- cors.mean / cors.var
        if (!is.finite(shape) || !is.finite(rate) || shape <= 0 || rate <= 0) {
            return(list(ok = FALSE))
        }
        list(ok = TRUE, shape = shape, rate = rate)
    }

    background_fit_correlations <- function(beta.mat, chr, pos, group.indices) {
        ord <- order(chr, pos)
        chr.ord <- chr[ord]
        pos.ord <- pos[ord]
        beta.ord <- beta.mat[ord, , drop = FALSE]
        cluster.ord <- bumphunter::clusterMaker(chr.ord, pos.ord, maxGap = 500)
        if (length(chr.ord) < 2L) {
            return(numeric(0))
        }
        pair.end <- which(
            chr.ord[-1L] == chr.ord[-length(chr.ord)] &
                cluster.ord[-1L] == cluster.ord[-length(cluster.ord)]
        ) + 1L
        pair.start <- pair.end - 1L
        if (length(pair.start) == 0L) {
            return(numeric(0))
        }
        pooled <- unlist(lapply(group.indices, function(idx) {
            if (length(idx) < 2L) {
                return(numeric(0))
            }
            m.values <- .transformBeta(
                beta.ord[, idx, drop = FALSE],
                pheno = data.frame(sample = seq_len(length(idx)))
            )
            group.cors <- compute_pair_correlations(
                m.values,
                start.ind = pair.start,
                end.ind = pair.end
            )
            group.cors[group.cors > 0]
        }))
        pooled[is.finite(pooled)]
    }

    pick_background_settings <- function(target, affected.samples, lookup.by.size) {
        grid <- lookup.by.size[[as.character(affected.samples)]]
        if (is.null(grid) || nrow(grid) == 0L) {
            return(list(
                corr.sd = default.corr.sd,
                sample.sd.frac = default.sample.sd.frac,
                metric = 0
            ))
        }
        if (is.na(target)) {
            chosen <- which(
                grid[, "corr.sd"] == default.corr.sd &
                    grid[, "sample.sd.frac"] == default.sample.sd.frac
            )
            if (length(chosen) == 0L) {
                chosen <- which.min(
                    abs(grid[, "corr.sd"] - default.corr.sd) +
                        abs(grid[, "sample.sd.frac"] - default.sample.sd.frac)
                )
            } else {
                chosen <- chosen[[1]]
            }
        } else {
            chosen <- which.min(abs(grid[, "metric"] - target))
        }
        list(
            corr.sd = grid[chosen, "corr.sd"],
            sample.sd.frac = grid[chosen, "sample.sd.frac"],
            metric = grid[chosen, "metric"]
        )
    }

    meth.mat <- as.matrix(bsseq::getCoverage(bs, type = "M"))
    unmeth.mat <- as.matrix(bsseq::getCoverage(bs, type = "Cov")) - meth.mat
    chr <- as.character(GenomeInfoDb::seqnames(bs))
    pos <- start(bs)

    cluster <- bumphunter::clusterMaker(chr, pos, maxGap = 500)
    Indexes <- split(seq(along = cluster), cluster)
    lns <- lengths(Indexes)
    Indexes <- Indexes[lns >= 5 & lns <= 500]

    # sample regions with intermediate methylation values preferentially
    prop.mat <- rowMeans(meth.mat / (meth.mat + unmeth.mat))
    prop.mat <- unlist(lapply(Indexes, function(x) median(prop.mat[x])))

    dmrs.ind <- sample(seq_along(Indexes),
        num.dmrs,
        replace = FALSE,
        prob = pmax(1 - sqrt(2) * abs(0.5 - prop.mat)^0.5, 0)
    )
    dmrs.ind <- Indexes[dmrs.ind]

    fnc <- function(index) {
        gr.dmr <- GenomicRanges::GRanges(
            seqnames = unique(as.character(GenomeInfoDb::seqnames(bs)[index])),
            IRanges::IRanges(
                start = min(start(bs)[index]),
                end = max(start(bs)[index])
            )
        )
        gr.dmr
    }

    gr.dmrs <- suppressWarnings(Reduce("c", lapply(dmrs.ind, fnc)))

    background.targets <- rep(NA_real_, num.dmrs)
    corr.targets <- rep(NA_real_, num.dmrs)
    corr.sd.used <- rep(NA_real_, num.dmrs)
    sample.sd.frac.used <- rep(NA_real_, num.dmrs)
    corr.mode.used <- rep(NA_character_, num.dmrs)
    background.lookup <- list()
    background.fit.ok <- FALSE
    if (use.correlated.effects && corr.mode == "background") {
        beta.mat <- meth.mat / (meth.mat + unmeth.mat)
        group.indices <- list(
            Condition1 = seq_len(sampleSize),
            Condition2 = if ((ncol(bs) - sampleSize) > 0L) {
                seq.int(sampleSize + 1L, ncol(bs))
            } else {
                integer(0)
            }
        )
        pos.cors <- background_fit_correlations(
            beta.mat = beta.mat,
            chr = chr,
            pos = pos,
            group.indices = group.indices
        )
        template.len <- max(5L, as.integer(round(stats::median(lengths(Indexes)))))
        template.spacing <- as.numeric(stats::median(
            unlist(lapply(Indexes, function(idx) diff(pos[idx]))),
            na.rm = TRUE
        ))
        if (!is.finite(template.spacing) || template.spacing <= 0) {
            template.spacing <- 100
        }
        template.pos <- seq.int(
            from = 1L,
            by = max(1L, as.integer(round(template.spacing))),
            length.out = template.len
        )
        template.cov <- max(1L, as.integer(round(stats::median(meth.mat + unmeth.mat))))
        background.lookup <- build_lookup_grid(
            affected.sample.sizes = sort(unique(c(sampleSize, ncol(bs) - sampleSize))),
            template.pos = template.pos,
            template.cov = template.cov,
            length.scale = corr.length
        )
        fit.ret <- fit_background_gamma(pos.cors)
        if (fit.ret$ok) {
            background.targets <- pmin(
                pmax(
                    stats::rgamma(
                        n = num.dmrs,
                        shape = fit.ret$shape,
                        rate = fit.ret$rate
                    ) * corr.rate,
                    0
                ),
                0.99
            )
            background.fit.ok <- TRUE
        } else {
            message(
                "Background neighboring-CpG correlation fit was degenerate; ",
                "using fallback manual settings corr.sd = ", default.corr.sd,
                " and sample.sd.frac = ", default.sample.sd.frac, "."
            )
        }
    }

    Diff <- rep(0, length(bs))
    dmr.mncov <- dmr.L <- deltas <- rep(NA, num.dmrs)

    for (u in seq_len(num.dmrs)) {
        up <- 1 - 2 * (rbinom(1, 1, 0.5) == 1)

        # let effect size change randomly
        delta.max <- delta.max0 + (rbeta(1, 2, 2) - 0.5) / 3
        deltas[u] <- delta.max

        idx <- dmrs.ind[[u]]
        dmr.L[u] <- length(idx)

        prop.mat <- meth.mat[idx, ] / (meth.mat[idx, ] + unmeth.mat[idx, ])

        # change direction if baseline mean is near boundary
        if (up == 1) {
            if (mean(prop.mat) > 1 - delta.max) {
                up <- -1
            }
        } else if (up == -1) {
            if (mean(prop.mat) < delta.max) {
                up <- 1
            }
        }

        # simulated mean as a smooth parabola added or subtracted from baseline mean
        last <- max(pos[idx])
        first <- min(pos[idx])
        width <- last - first

        # widen so first and last CpGs do not have zero difference
        last <- last + 0.2 * width
        first <- first - 0.2 * width
        width <- last - first
        mid <- round((last - first) / 2 + first)

        Diff.hit <- round(
            triwt(pos[idx],
                amp = delta.max,
                base = Diff[idx],
                width = width,
                center = mid,
                deg = 4,
                dir = up
            ),
            4
        )

        mn.cov <- by(
            t(meth.mat[idx, ] + unmeth.mat[idx, ]),
            factor(paste0(
                "Condition",
                c(
                    rep(1, sampleSize),
                    rep(2, ncol(bs) - sampleSize)
                )
            )),
            colMeans
        )
        mn.cov <- rowMeans(cbind(mn.cov[[1]], mn.cov[[2]]))
        dmr.mncov[u] <- mean(mn.cov)

        cov <- meth.mat[idx, ] + unmeth.mat[idx, ]
        prop <- meth.mat[idx, ] / cov

        # randomly choose which condition gets the DMR
        grp <- runif(1) < 0.5
        ss <- ifelse(grp, sampleSize, ncol(bs) - sampleSize)

        dmr.pos <- pos[idx]
        current.corr.sd <- 0
        current.sample.sd.frac <- 0
        if (use.correlated.effects) {
            if (corr.mode == "background") {
                picked <- pick_background_settings(
                    target = if (background.fit.ok) background.targets[u] else NA_real_,
                    affected.samples = ss,
                    lookup.by.size = background.lookup
                )
                current.corr.sd <- picked$corr.sd
                current.sample.sd.frac <- picked$sample.sd.frac
                corr.targets[u] <- if (background.fit.ok) background.targets[u] else picked$metric
                corr.mode.used[u] <- "background"
            } else {
                current.corr.sd <- corr.sd
                current.sample.sd.frac <- sample.sd.frac
                corr.mode.used[u] <- "manual"
            }
        } else {
            corr.mode.used[u] <- "disabled"
        }
        corr.sd.used[u] <- current.corr.sd
        sample.sd.frac.used[u] <- current.sample.sd.frac

        # shared correlated field across samples in the affected group
        # gives region-level coherence
        group.field <- if (use.correlated.effects) {
            make_correlated_effect(dmr.pos,
                sigma = current.corr.sd,
                length.scale = corr.length
            )
        } else {
            rep(0, length(idx))
        }

        # smaller sample-specific field for replicate variation
        sample.sd <- current.corr.sd * current.sample.sd.frac

        for (samp in seq_len(ss)) {
            if (grp) {
                new.p <- make_shifted_prob(
                    base.p = prop[, samp],
                    diff.hit = Diff.hit,
                    pos = dmr.pos,
                    group.field = group.field,
                    sample.sd = sample.sd,
                    length.scale = corr.length,
                    use.cor = use.correlated.effects
                )

                meth.mat[idx, samp] <- rbinom(
                    n = length(idx),
                    size = cov[, samp],
                    prob = new.p
                )
                unmeth.mat[idx, samp] <- cov[, samp] - meth.mat[idx, samp]
            } else {
                j <- sampleSize + samp

                new.p <- make_shifted_prob(
                    base.p = prop[, j],
                    diff.hit = Diff.hit,
                    pos = dmr.pos,
                    group.field = group.field,
                    sample.sd = sample.sd,
                    length.scale = corr.length,
                    use.cor = use.correlated.effects
                )

                meth.mat[idx, j] <- rbinom(
                    n = length(idx),
                    size = cov[, j],
                    prob = new.p
                )
                unmeth.mat[idx, j] <- cov[, j] - meth.mat[idx, j]
            }
        }
    }

    S4Vectors::mcols(gr.dmrs)$corr_target <- corr.targets
    S4Vectors::mcols(gr.dmrs)$corr_sd_used <- corr.sd.used
    S4Vectors::mcols(gr.dmrs)$sample_sd_frac_used <- sample.sd.frac.used
    S4Vectors::mcols(gr.dmrs)$corr_mode_used <- corr.mode.used

    sampnames <- paste0(
        "Condition",
        c(
            rep(1, sampleSize),
            rep(2, ncol(bs) - sampleSize)
        ),
        "_Rep",
        c(
            seq_len(sampleSize),
            seq_len(ncol(bs) - sampleSize)
        )
    )

    colnames(meth.mat) <- colnames(unmeth.mat) <- sampnames

    bsNew <- bsseq::BSseq(
        pos = pos, chr = chr, M = meth.mat,
        Cov = (meth.mat + unmeth.mat),
        sampleNames = sampnames
    )

    sim.dat.red <- list(
        gr.dmrs = gr.dmrs,
        dmr.mncov = dmr.mncov,
        dmr.L = dmr.L,
        bs = bsNew,
        delta = deltas
    )

    sim.dat.red
}

# nolint end
