# Stabilized ICA ---------------------------------------------------------------

.fastica_once <- function(z, k, seed, maxit = 200L) {
  set.seed(as.integer(seed))
  X <- t(z) # samples x genes
  fit <- fastICA::fastICA(
    X, n.comp = k, alg.typ = "parallel", fun = "logcosh",
    alpha = 1, method = "C", row.norm = FALSE,
    maxit = as.integer(maxit), tol = 1e-4, verbose = FALSE
  )
  # fastICA uses X = S A. A is k x genes; transpose to genes x components.
  L <- t(fit$A)
  rownames(L) <- colnames(X)
  colnames(L) <- paste0("IC", seq_len(k))
  for (j in seq_len(ncol(L))) L[, j] <- .orient_by_max_abs(L[, j])
  list(loadings = L)
}

.project_scores_from_loadings <- function(z, loadings) {
  # Stable components are selected mixing vectors from different ICA restarts.
  # Re-estimate sample coordinates by least squares in that stable loading basis.
  X <- t(z) # samples x genes
  scores <- tryCatch(
    t(qr.solve(loadings, t(X))),
    error = function(e) NULL
  )
  if (is.null(scores)) {
    stop("Stable loading basis is rank-deficient; cannot project sample scores.", call. = FALSE)
  }
  rownames(scores) <- rownames(X)
  colnames(scores) <- colnames(loadings)
  scores <- scale(scores)
  scores[!is.finite(scores)] <- 0
  scores
}

.cluster_stability_summary <- function(big, C, cluster, run_of, n_runs_total) {
  out <- vector("list", length(unique(cluster)))
  cluster_ids <- sort(unique(cluster))
  for (ii in seq_along(cluster_ids)) {
    ci <- cluster_ids[ii]
    idx <- which(cluster == ci)
    if (!length(idx)) next

    # Preliminary medoid. Singleton is deliberately NOT assigned stability 1.
    if (length(idx) == 1L) {
      prelim <- idx
    } else {
      subC <- C[idx, idx, drop = FALSE]
      diag(subC) <- NA_real_
      prelim <- idx[which.max(rowMeans(subC, na.rm = TRUE))]
    }

    # One estimate per ICA restart: prevents duplicate components from one run
    # inflating component stability.
    idx_one_per_run <- unlist(lapply(unique(run_of[idx]), function(rr) {
      cand <- idx[run_of[idx] == rr]
      if (length(cand) == 1L) return(cand)
      cand[which.max(C[cand, prelim])]
    }), use.names = FALSE)

    n_runs_represented <- length(unique(run_of[idx_one_per_run]))
    run_fraction <- n_runs_represented / n_runs_total

    if (length(idx_one_per_run) < 2L) {
      stability <- 0
      medoid <- idx_one_per_run[1]
    } else {
      subC <- C[idx_one_per_run, idx_one_per_run, drop = FALSE]
      diag(subC) <- NA_real_
      mean_cor <- rowMeans(subC, na.rm = TRUE)
      medoid <- idx_one_per_run[which.max(mean_cor)]
      vals <- subC[upper.tri(subC)]
      vals <- vals[is.finite(vals)]
      stability <- if (length(vals)) mean(vals) else 0
    }

    out[[ii]] <- list(
      cluster_id = ci,
      medoid_index = medoid,
      stability = stability,
      n_estimates_raw = length(idx),
      n_estimates_balanced = length(idx_one_per_run),
      n_runs_represented = n_runs_represented,
      run_fraction = run_fraction
    )
  }
  out
}

#' Stabilized ICA for one cohort and one k
#'
#' Repeats FastICA, clusters component estimates by absolute loading
#' correlation, chooses a centrotype/medoid, and requires both within-cluster
#' similarity and recurrence across ICA restarts. Multiple estimates from the
#' same restart receive only one vote in the stability calculation.
#'
#' @param x Numeric matrix, genes x samples.
#' @param k Number of ICA components.
#' @param n_runs Number of random restarts (`>=50` recommended for production).
#' @param seed Base random seed.
#' @param min_stability Minimum mean pairwise absolute loading correlation.
#' @param min_run_fraction Minimum fraction of successful restarts represented
#'   by a stable component cluster.
#' @param gene_zscore If TRUE, z-score each gene across samples before ICA.
#' @param maxit Maximum FastICA iterations.
#' @param ... Reserved; outcome-like argument names are rejected.
#' @return A list containing stable loadings, projected sample scores, component
#'   stability metadata and run metadata.
#' @export
run_stabilized_ica <- function(
    x,
    k,
    n_runs = 50L,
    seed = 1L,
    min_stability = 0.55,
    min_run_fraction = 0.80,
    gene_zscore = TRUE,
    maxit = 200L,
    ...
) {
  assert_no_outcome_args(...)
  x <- .validate_expression_matrix(x, "x")
  k <- as.integer(k)
  n_runs <- as.integer(n_runs)
  if (length(k) != 1L || is.na(k) || k < 2L || k >= min(nrow(x), ncol(x))) {
    stop("k must be >=2 and < min(n_genes, n_samples).", call. = FALSE)
  }
  if (n_runs < 5L) stop("n_runs must be >= 5.", call. = FALSE)
  if (!is.finite(min_stability) || min_stability < 0 || min_stability > 1) {
    stop("min_stability must be in [0, 1].", call. = FALSE)
  }
  if (!is.finite(min_run_fraction) || min_run_fraction <= 0 || min_run_fraction > 1) {
    stop("min_run_fraction must be in (0, 1].", call. = FALSE)
  }

  sds <- matrixStats::rowSds(x, na.rm = TRUE)
  x <- x[is.finite(sds) & sds > 0, , drop = FALSE]
  if (nrow(x) <= k) stop("Too few variable genes after zero-SD filtering.", call. = FALSE)

  if (gene_zscore) {
    z <- t(scale(t(x)))
  } else {
    z <- x
  }
  z[!is.finite(z)] <- 0

  gene_names <- rownames(z)
  estimates <- list()
  run_of <- integer()
  source_component <- character()
  ok_runs <- 0L
  failed_runs <- integer()

  for (r in seq_len(n_runs)) {
    res <- tryCatch(.fastica_once(z, k, seed + r, maxit = maxit), error = function(e) NULL)
    if (is.null(res)) {
      failed_runs <- c(failed_runs, r)
      next
    }
    ok_runs <- ok_runs + 1L
    for (j in seq_len(ncol(res$loadings))) {
      estimates[[length(estimates) + 1L]] <- res$loadings[, j]
      run_of <- c(run_of, r)
      source_component <- c(source_component, colnames(res$loadings)[j])
    }
  }

  if (ok_runs < max(5L, ceiling(n_runs * min_run_fraction))) {
    stop(
      "Too few successful ICA runs to satisfy min_run_fraction: ",
      ok_runs, "/", n_runs,
      call. = FALSE
    )
  }

  big <- do.call(cbind, estimates)
  rownames(big) <- gene_names
  colnames(big) <- paste0("run", run_of, "_", source_component)

  # Pearson correlation of loading vectors is used only for stabilization;
  # cross-cohort alignment later uses Spearman on gene loadings.
  B <- scale(big)
  B[!is.finite(B)] <- 0
  C <- abs(crossprod(B) / max(1, nrow(B) - 1L))
  C[C > 1] <- 1
  diag(C) <- 1
  cluster <- stats::cutree(
    stats::hclust(stats::as.dist(1 - C), method = "average"),
    k = k
  )

  summaries <- .cluster_stability_summary(big, C, cluster, run_of, n_runs)
  meta_all <- dplyr::bind_rows(lapply(summaries, function(s) {
    tibble::tibble(
      cluster_id = as.integer(s$cluster_id),
      medoid_index = as.integer(s$medoid_index),
      stability = as.numeric(s$stability),
      n_estimates_raw = as.integer(s$n_estimates_raw),
      n_estimates_balanced = as.integer(s$n_estimates_balanced),
      n_runs_represented = as.integer(s$n_runs_represented),
      run_fraction = as.numeric(s$run_fraction)
    )
  }))
  meta_all$passed <- with(
    meta_all,
    is.finite(stability) & stability >= min_stability &
      is.finite(run_fraction) & run_fraction >= min_run_fraction
  )

  kept <- meta_all[meta_all$passed, , drop = FALSE]
  if (!nrow(kept)) {
    stop("No components passed stability + run-coverage gates.", call. = FALSE)
  }

  L <- big[, kept$medoid_index, drop = FALSE]
  colnames(L) <- paste0("sIC", seq_len(ncol(L)))
  for (j in seq_len(ncol(L))) L[, j] <- .orient_by_max_abs(L[, j])
  kept$component <- colnames(L)
  kept <- kept[, c(
    "component", "cluster_id", "stability", "run_fraction",
    "n_runs_represented", "n_estimates_raw", "n_estimates_balanced", "passed"
  )]

  scores <- .project_scores_from_loadings(z, L)
  names_stability <- stats::setNames(kept$stability, kept$component)
  names_run_fraction <- stats::setNames(kept$run_fraction, kept$component)

  list(
    loadings = L,
    scores = scores,
    score_definition = "least_squares_projection_on_stable_mixing_vectors",
    stability = names_stability,
    run_fraction = names_run_fraction,
    component_meta = tibble::as_tibble(kept),
    all_cluster_meta = tibble::as_tibble(meta_all),
    mean_stability = mean(kept$stability),
    mean_run_fraction = mean(kept$run_fraction),
    stable_fraction = nrow(kept) / k,
    k_requested = k,
    n_kept = nrow(kept),
    n_runs_requested = n_runs,
    n_runs_ok = ok_runs,
    failed_runs = failed_runs,
    genes = gene_names,
    samples = colnames(z)
  )
}

#' Scan a k grid and select a stable decomposition
#'
#' @param x Numeric genes x samples matrix.
#' @param k_grid Candidate k values.
#' @param n_runs ICA restarts per k.
#' @param seed Base seed.
#' @param min_stability Stability threshold.
#' @param min_run_fraction Restart-coverage threshold.
#' @param selection_metric Either `"stability_coverage"` (default; penalizes a
#'   k that retains only a few stable components) or `"mean_stability"`.
#' @param ... Forwarded to [run_stabilized_ica()]; outcome-like names rejected.
#' @return List with `best` fit and a scan table.
#' @export
scan_ica_k <- function(
    x,
    k_grid = c(10L, 15L, 20L, 25L, 30L),
    n_runs = 50L,
    seed = 1L,
    min_stability = 0.55,
    min_run_fraction = 0.80,
    selection_metric = c("stability_coverage", "mean_stability"),
    ...
) {
  assert_no_outcome_args(...)
  x <- .validate_expression_matrix(x, "x")
  selection_metric <- match.arg(selection_metric)
  ks <- sort(unique(as.integer(k_grid)))
  ks <- ks[ks >= 2L & ks < min(nrow(x), ncol(x))]
  if (!length(ks)) stop("No valid k values.", call. = FALSE)

  fits <- list()
  scan <- list()
  for (k in ks) {
    fit <- tryCatch(
      run_stabilized_ica(
        x, k = k, n_runs = n_runs, seed = seed,
        min_stability = min_stability,
        min_run_fraction = min_run_fraction, ...
      ),
      error = function(e) NULL
    )
    if (is.null(fit)) next
    score <- if (selection_metric == "stability_coverage") {
      fit$mean_stability * fit$stable_fraction
    } else {
      fit$mean_stability
    }
    fits[[as.character(k)]] <- fit
    scan[[length(scan) + 1L]] <- tibble::tibble(
      k = k,
      mean_stability = fit$mean_stability,
      mean_run_fraction = fit$mean_run_fraction,
      stable_fraction = fit$stable_fraction,
      n_kept = fit$n_kept,
      n_runs_ok = fit$n_runs_ok,
      selection_score = score
    )
  }
  if (!length(fits)) stop("No successful k in scan.", call. = FALSE)

  scan_tbl <- dplyr::bind_rows(scan)
  ord <- order(
    scan_tbl$selection_score,
    scan_tbl$mean_run_fraction,
    scan_tbl$n_kept,
    -scan_tbl$k,
    decreasing = TRUE
  )
  selected_k <- scan_tbl$k[ord[1]]
  best <- fits[[as.character(selected_k)]]
  scan_tbl$selected <- scan_tbl$k == selected_k
  scan_tbl$selection_metric <- selection_metric
  list(best = best, scan = scan_tbl)
}

#' Summarize stability of retained components
#' @param fit Object from [run_stabilized_ica()] or `scan_ica_k()$best`.
#' @return Component-level stability table.
#' @export
summarize_component_stability <- function(fit) {
  out <- fit$component_meta
  out$k_requested <- fit$k_requested
  out$n_runs_ok <- fit$n_runs_ok
  out
}

#' Re-filter stable components
#' @param fit Stabilized ICA fit.
#' @param min_stability Minimum stability.
#' @param min_run_fraction Minimum restart coverage.
#' @return Filtered fit.
#' @export
select_stable_components <- function(
    fit,
    min_stability = 0.55,
    min_run_fraction = 0.80
) {
  keep <- which(
    is.finite(fit$stability) & fit$stability >= min_stability &
      is.finite(fit$run_fraction) & fit$run_fraction >= min_run_fraction
  )
  if (!length(keep)) stop("No stable components.", call. = FALSE)
  fit$loadings <- fit$loadings[, keep, drop = FALSE]
  fit$scores <- fit$scores[, keep, drop = FALSE]
  fit$stability <- fit$stability[keep]
  fit$run_fraction <- fit$run_fraction[keep]
  fit$component_meta <- fit$component_meta[keep, , drop = FALSE]
  fit$n_kept <- length(keep)
  fit$mean_stability <- mean(fit$stability)
  fit$mean_run_fraction <- mean(fit$run_fraction)
  fit$stable_fraction <- fit$n_kept / fit$k_requested
  fit
}
