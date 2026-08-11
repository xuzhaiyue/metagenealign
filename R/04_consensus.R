# Unit-balanced consensus metagenes -------------------------------------------

#' Orient member loadings to a medoid by signed Spearman correlation
#'
#' @param loadings Named list of loading vectors.
#' @param member_ids Member node IDs.
#' @param medoid_id Optional medoid node ID.
#' @param min_shared_genes Minimum genes for sign correlation.
#' @return List with oriented loadings, medoid and signs.
#' @export
orient_component_signs <- function(
    loadings,
    member_ids,
    medoid_id = NULL,
    min_shared_genes = 50L
) {
  member_ids <- unique(member_ids[member_ids %in% names(loadings)])
  if (!length(member_ids)) stop("No member loadings found.", call. = FALSE)
  if (!is.null(medoid_id) && !medoid_id %in% member_ids) {
    stop("medoid_id is not a member.", call. = FALSE)
  }

  if (is.null(medoid_id)) {
    if (length(member_ids) == 1L) {
      medoid_id <- member_ids[[1]]
    } else {
      M <- matrix(NA_real_, length(member_ids), length(member_ids),
                  dimnames = list(member_ids, member_ids))
      diag(M) <- 1
      for (i in seq_len(length(member_ids) - 1L)) {
        for (j in seq.int(i + 1L, length(member_ids))) {
          r <- .abs_spearman(
            loadings[[member_ids[i]]], loadings[[member_ids[j]]],
            min_shared = min_shared_genes
          )
          M[i, j] <- r; M[j, i] <- r
        }
      }
      score <- rowMeans(M, na.rm = TRUE)
      score[!is.finite(score)] <- -Inf
      medoid_id <- member_ids[which.max(score)]
    }
  }

  signs <- stats::setNames(rep(1, length(member_ids)), member_ids)
  oriented <- lapply(member_ids, function(m) {
    v <- loadings[[m]]
    rho <- .signed_spearman(v, loadings[[medoid_id]], min_shared = min_shared_genes)
    if (is.finite(rho) && rho < 0) {
      signs[[m]] <<- -1
      v <- -v
    }
    v
  })
  names(oriented) <- member_ids
  list(oriented = oriented, medoid = medoid_id, signs = signs)
}

#' Collapse multiple components from one analysis unit
#'
#' @param oriented Named list of sign-oriented loading vectors.
#' @param unit_of Named node ID to analysis-unit map.
#' @param aggregate Within-unit aggregation (`median` or `mean`).
#' @return Named list containing one loading vector per analysis unit.
#' @export
collapse_components_within_unit <- function(
    oriented,
    unit_of,
    aggregate = c("median", "mean")
) {
  aggregate <- match.arg(aggregate)
  ids <- names(oriented)
  if (is.null(ids) || any(!nzchar(ids))) stop("oriented must be named.", call. = FALSE)
  units <- as.character(unit_of[ids])
  if (any(is.na(units) | !nzchar(units))) stop("unit_of missing for oriented members.", call. = FALSE)

  out <- lapply(sort(unique(units)), function(u) {
    mid <- ids[units == u]
    genes <- sort(unique(unlist(lapply(mid, function(id) names(oriented[[id]])))))
    mat <- do.call(cbind, lapply(mid, function(id) {
      v <- rep(NA_real_, length(genes)); names(v) <- genes
      shared <- intersect(names(oriented[[id]]), genes)
      v[shared] <- oriented[[id]][shared]
      v
    }))
    rownames(mat) <- genes
    if (aggregate == "median") {
      apply(mat, 1, stats::median, na.rm = TRUE)
    } else {
      rowMeans(mat, na.rm = TRUE)
    }
  })
  names(out) <- sort(unique(units))
  out
}

#' Build a unit-balanced consensus metagene
#'
#' The workflow is sign orientation -> within-unit collapse -> cross-unit
#' aggregation. Genes must be observed in a minimum fraction of voting units,
#' preventing a gene present in only one cohort/component from defining the
#' cross-cohort consensus.
#'
#' @param loadings Named list of all component loadings.
#' @param member_ids Member nodes of one program/community.
#' @param unit_of Named node-to-analysis-unit map.
#' @param aggregate_within Within-unit aggregation.
#' @param aggregate_across Across-unit aggregation.
#' @param min_gene_unit_fraction Minimum fraction of voting units in which a
#'   gene must have a finite loading.
#' @return Consensus vector plus support/provenance metadata.
#' @export
build_consensus_metagene <- function(
    loadings,
    member_ids,
    unit_of,
    aggregate_within = "median",
    aggregate_across = c("median", "mean"),
    min_gene_unit_fraction = 0.50
) {
  aggregate_across <- match.arg(aggregate_across)
  if (!is.finite(min_gene_unit_fraction) || min_gene_unit_fraction <= 0 || min_gene_unit_fraction > 1) {
    stop("min_gene_unit_fraction must be in (0, 1].", call. = FALSE)
  }
  ori <- orient_component_signs(loadings, member_ids)
  unit_load <- collapse_components_within_unit(
    ori$oriented, unit_of = unit_of, aggregate = aggregate_within
  )
  genes <- sort(unique(unlist(lapply(unit_load, names))))
  matU <- do.call(cbind, lapply(unit_load, function(v) {
    out <- rep(NA_real_, length(genes)); names(out) <- genes
    out[names(v)] <- v
    out
  }))
  rownames(matU) <- genes
  colnames(matU) <- names(unit_load)

  n_support <- rowSums(is.finite(matU))
  min_support <- ceiling(ncol(matU) * min_gene_unit_fraction)
  keep <- n_support >= min_support
  if (!any(keep)) stop("No genes pass the cross-unit support gate.", call. = FALSE)
  matU_keep <- matU[keep, , drop = FALSE]

  cons <- if (aggregate_across == "median") {
    apply(matU_keep, 1, stats::median, na.rm = TRUE)
  } else {
    rowMeans(matU_keep, na.rm = TRUE)
  }
  cons <- cons[is.finite(cons)]
  cons <- .orient_by_max_abs(cons)

  gene_support <- tibble::tibble(
    gene = rownames(matU),
    n_units = as.integer(n_support),
    unit_fraction = n_support / ncol(matU),
    retained = keep
  )

  list(
    consensus = cons,
    gene_support = gene_support,
    medoid = ori$medoid,
    n_member_nodes = length(ori$oriented),
    n_units_voting = length(unit_load),
    units = names(unit_load),
    signs = ori$signs,
    min_gene_unit_fraction = min_gene_unit_fraction
  )
}
