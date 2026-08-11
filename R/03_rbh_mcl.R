# Cross-cohort RBH + MCL -------------------------------------------------------

#' Pairwise absolute Spearman similarity of gene loadings
#'
#' @param loadings Named list of named numeric loading vectors.
#' @param min_shared_genes Minimum shared genes required.
#' @return Square matrix of absolute Spearman correlations.
#' @export
calculate_loading_similarity <- function(loadings, min_shared_genes = 200L) {
  ids <- names(loadings)
  if (is.null(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) {
    stop("loadings must have unique non-empty names.", call. = FALSE)
  }
  n <- length(ids)
  S <- matrix(NA_real_, n, n, dimnames = list(ids, ids))
  if (!n) return(S)
  diag(S) <- 1
  if (n >= 2L) {
    for (i in seq_len(n - 1L)) {
      for (j in seq.int(i + 1L, n)) {
        r <- .abs_spearman(loadings[[ids[i]]], loadings[[ids[j]]], min_shared = min_shared_genes)
        S[i, j] <- r
        S[j, i] <- r
      }
    }
  }
  S
}

.best_hits <- function(x, labels, tolerance = 1e-12) {
  ok <- is.finite(x)
  if (!any(ok)) return(character())
  best <- max(x[ok])
  labels[ok & abs(x - best) <= tolerance]
}

#' Build a threshold-free reciprocal-best-hit graph
#'
#' RBH is evaluated separately for every pair of analysis units. Exact ties can
#' be retained rather than resolved by input order.
#'
#' @param loadings Named list of loading vectors.
#' @param unit_of Named character vector mapping node ID to analysis unit.
#' @param min_shared_genes Minimum genes shared by two loading vectors.
#' @param similarity Optional precomputed absolute-correlation matrix with node
#'   IDs as row/column names.
#' @param tie_policy `"all"` retains all exact/near ties; `"first"` uses a
#'   deterministic first match.
#' @param tie_tolerance Numeric tolerance used for tied best hits.
#' @return `list(edges, graph)`.
#' @export
build_rbh_graph <- function(
    loadings,
    unit_of,
    min_shared_genes = 200L,
    similarity = NULL,
    tie_policy = c("all", "first"),
    tie_tolerance = 1e-12
) {
  tie_policy <- match.arg(tie_policy)
  ids <- names(loadings)
  if (is.null(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) {
    stop("loadings must have unique non-empty names.", call. = FALSE)
  }
  if (is.null(names(unit_of))) stop("unit_of must be named by node ID.", call. = FALSE)
  unit_of <- as.character(unit_of[ids])
  names(unit_of) <- ids
  if (any(is.na(unit_of) | !nzchar(unit_of))) {
    stop("unit_of missing for some loading vectors.", call. = FALSE)
  }

  if (!is.null(similarity)) {
    similarity <- as.matrix(similarity)
    if (is.null(rownames(similarity)) || is.null(colnames(similarity)) ||
        !all(ids %in% rownames(similarity)) || !all(ids %in% colnames(similarity))) {
      stop("similarity must contain all node IDs in row/column names.", call. = FALSE)
    }
  }

  units <- unique(unit_of)
  edge_rows <- list()
  if (length(units) >= 2L) {
    for (ii in seq_len(length(units) - 1L)) {
      for (jj in seq.int(ii + 1L, length(units))) {
        u1 <- units[ii]; u2 <- units[jj]
        n1 <- ids[unit_of == u1]
        n2 <- ids[unit_of == u2]
        S <- if (!is.null(similarity)) {
          similarity[n1, n2, drop = FALSE]
        } else {
          tmp <- matrix(NA_real_, length(n1), length(n2), dimnames = list(n1, n2))
          for (a in n1) for (b in n2) {
            tmp[a, b] <- .abs_spearman(loadings[[a]], loadings[[b]], min_shared = min_shared_genes)
          }
          tmp
        }
        if (!length(S) || all(is.na(S))) next

        best12 <- lapply(seq_len(nrow(S)), function(i) {
          h <- .best_hits(S[i, ], colnames(S), tie_tolerance)
          if (tie_policy == "first" && length(h)) h[1] else h
        })
        names(best12) <- rownames(S)
        best21 <- lapply(seq_len(ncol(S)), function(j) {
          h <- .best_hits(S[, j], rownames(S), tie_tolerance)
          if (tie_policy == "first" && length(h)) h[1] else h
        })
        names(best21) <- colnames(S)

        for (a in n1) {
          for (b in best12[[a]]) {
            if (a %in% best21[[b]]) {
              edge_rows[[length(edge_rows) + 1L]] <- tibble::tibble(
                from = a,
                to = b,
                unit_from = u1,
                unit_to = u2,
                abs_spearman = as.numeric(S[a, b]),
                signed_spearman = .signed_spearman(
                  loadings[[a]], loadings[[b]], min_shared = min_shared_genes
                )
              )
            }
          }
        }
      }
    }
  }

  edges <- if (length(edge_rows)) {
    dplyr::bind_rows(edge_rows)
  } else {
    tibble::tibble(
      from = character(), to = character(),
      unit_from = character(), unit_to = character(),
      abs_spearman = double(), signed_spearman = double()
    )
  }
  if (nrow(edges)) {
    edges <- dplyr::distinct(edges, .data$from, .data$to, .keep_all = TRUE)
  }
  if (!nrow(edges)) {
    g <- igraph::make_empty_graph(n = length(ids), directed = FALSE)
    igraph::V(g)$name <- ids
    return(list(edges = edges, graph = g))
  }
  g <- igraph::graph_from_data_frame(
    edges[, c("from", "to", "abs_spearman")],
    vertices = tibble::tibble(name = ids, analysis_unit = unname(unit_of[ids])),
    directed = FALSE
  )
  list(edges = edges, graph = g)
}

#' Cluster an RBH graph with Markov Clustering
#'
#' Connected components are intentionally not used as biological programs.
#'
#' @param edges Edge table with `from`, `to`, `abs_spearman`.
#' @param node_ids Character vector of all node IDs, including isolates.
#' @param inflation MCL inflation parameter (>0).
#' @param expansion MCL expansion parameter (>1).
#' @param allow1 Allow singleton clusters; singletons will fail recurrence later.
#' @return Named integer community membership.
#' @export
cluster_rbh_mcl <- function(
    edges,
    node_ids,
    inflation = 2,
    expansion = 2,
    allow1 = TRUE
) {
  node_ids <- unique(as.character(node_ids))
  if (!length(node_ids)) return(stats::setNames(integer(), character()))
  if (!is.finite(inflation) || inflation <= 0) stop("inflation must be > 0.", call. = FALSE)
  if (!is.finite(expansion) || expansion <= 1) stop("expansion must be > 1.", call. = FALSE)
  if (!all(c("from", "to", "abs_spearman") %in% names(edges))) {
    stop("edges must contain from, to, abs_spearman.", call. = FALSE)
  }

  n <- length(node_ids)
  A <- matrix(0, n, n, dimnames = list(node_ids, node_ids))
  if (nrow(edges)) {
    for (r in seq_len(nrow(edges))) {
      a <- as.character(edges$from[r]); b <- as.character(edges$to[r])
      w <- as.numeric(edges$abs_spearman[r])
      if (!a %in% node_ids || !b %in% node_ids || !is.finite(w) || w < 0) next
      A[a, b] <- max(A[a, b], w)
      A[b, a] <- max(A[b, a], w)
    }
  }
  diag(A) <- 0
  fit <- MCL::mcl(
    A,
    addLoops = TRUE,
    expansion = expansion,
    inflation = inflation,
    allow1 = allow1,
    max.iter = 200,
    ESM = FALSE
  )
  cl <- as.integer(fit$Cluster)
  if (length(cl) != n) stop("MCL returned unexpected membership length.", call. = FALSE)
  names(cl) <- node_ids
  attr(cl, "method") <- "Markov_Clustering_MCL"
  attr(cl, "not_used") <- "igraph_connected_components"
  cl
}

#' Summarize MCL communities and apply unique-unit recurrence
#'
#' @param membership Named community IDs.
#' @param node_meta Data frame containing at least `node_id`, `analysis_unit`;
#'   optional `source_study` and `platform_class` columns are audited.
#' @param edges Optional RBH edge table.
#' @param min_unique_units Absolute unique-unit recurrence gate.
#' @param n_units_total Recurrence denominator.
#' @return List with node assignments, community ledger and recurrent programs.
#' @export
summarize_program_recurrence <- function(
    membership,
    node_meta,
    edges = NULL,
    min_unique_units = 6L,
    n_units_total = NULL
) {
  node_meta <- as.data.frame(node_meta)
  if (!all(c("node_id", "analysis_unit") %in% names(node_meta))) {
    stop("node_meta needs node_id and analysis_unit.", call. = FALSE)
  }
  if (anyDuplicated(node_meta$node_id)) stop("node_meta node_id must be unique.", call. = FALSE)
  if (is.null(names(membership))) stop("membership must be named by node_id.", call. = FALSE)
  memb <- as.integer(membership[as.character(node_meta$node_id)])
  if (any(is.na(memb))) stop("membership missing for some node_meta rows.", call. = FALSE)

  tmp <- node_meta
  tmp$cluster_raw <- memb
  cmap <- tmp |>
    dplyr::count(.data$cluster_raw, name = "n_nodes") |>
    dplyr::arrange(dplyr::desc(.data$n_nodes), .data$cluster_raw) |>
    dplyr::mutate(community_id = paste0("C", dplyr::row_number()))
  tmp <- dplyr::left_join(tmp, cmap, by = "cluster_raw")

  n_units_total <- if (is.null(n_units_total)) {
    dplyr::n_distinct(tmp$analysis_unit)
  } else {
    as.integer(n_units_total)
  }

  edge_stats <- NULL
  if (!is.null(edges) && nrow(edges)) {
    edge_stats <- dplyr::bind_rows(lapply(unique(tmp$community_id), function(cid) {
      memb_nodes <- tmp$node_id[tmp$community_id == cid]
      subE <- edges[edges$from %in% memb_nodes & edges$to %in% memb_nodes, , drop = FALSE]
      possible <- length(memb_nodes) * (length(memb_nodes) - 1) / 2
      tibble::tibble(
        community_id = cid,
        n_edges = nrow(subE),
        edge_density = if (possible > 0) nrow(subE) / possible else 0,
        median_abs_rho = if (nrow(subE)) stats::median(subE$abs_spearman, na.rm = TRUE) else NA_real_,
        min_abs_rho = if (nrow(subE)) min(subE$abs_spearman, na.rm = TRUE) else NA_real_,
        max_abs_rho = if (nrow(subE)) max(subE$abs_spearman, na.rm = TRUE) else NA_real_
      )
    }))
  }

  qc <- tmp |>
    dplyr::group_by(.data$community_id) |>
    dplyr::summarise(
      n_nodes = dplyr::n(),
      n_unique_units = dplyr::n_distinct(.data$analysis_unit),
      max_nodes_per_unit = max(as.integer(table(.data$analysis_unit))),
      median_nodes_per_unit = as.numeric(stats::median(as.integer(table(.data$analysis_unit)))),
      duplicate_unit_fraction = mean(as.integer(table(.data$analysis_unit)) > 1L),
      unit_ids = paste(sort(unique(.data$analysis_unit)), collapse = ";"),
      n_unique_source_studies = if ("source_study" %in% names(tmp)) {
        dplyr::n_distinct(.data$source_study[!is.na(.data$source_study) & nzchar(.data$source_study)])
      } else NA_integer_,
      source_study_ids = if ("source_study" %in% names(tmp)) {
        paste(sort(unique(.data$source_study[!is.na(.data$source_study) & nzchar(.data$source_study)])), collapse = ";")
      } else NA_character_,
      n_platform_classes = if ("platform_class" %in% names(tmp)) {
        dplyr::n_distinct(.data$platform_class[!is.na(.data$platform_class) & nzchar(.data$platform_class)])
      } else NA_integer_,
      platform_classes = if ("platform_class" %in% names(tmp)) {
        paste(sort(unique(.data$platform_class[!is.na(.data$platform_class) & nzchar(.data$platform_class)])), collapse = ";")
      } else NA_character_,
      .groups = "drop"
    ) |>
    dplyr::mutate(
      unit_fraction = .data$n_unique_units / n_units_total,
      qc_flag = dplyr::case_when(
        .data$max_nodes_per_unit >= 4L ~ "HIGH_DUPLICATE_UNIT_BURDEN",
        .data$max_nodes_per_unit >= 2L ~ "INSPECT_MULTI_COMPONENT_PER_UNIT",
        TRUE ~ "OK"
      )
    )

  qc$recurrent <- mapply(
    assert_valid_recurrence,
    n_nodes = qc$n_nodes,
    n_unique_units = qc$n_unique_units,
    MoreArgs = list(n_units_total = n_units_total, min_unique_units = min_unique_units)
  )
  if (!is.null(edge_stats)) qc <- dplyr::left_join(qc, edge_stats, by = "community_id")
  if (!"n_edges" %in% names(qc)) qc$n_edges <- 0L
  if (!"edge_density" %in% names(qc)) qc$edge_density <- 0
  if (!"median_abs_rho" %in% names(qc)) qc$median_abs_rho <- NA_real_
  if (!"min_abs_rho" %in% names(qc)) qc$min_abs_rho <- NA_real_
  if (!"max_abs_rho" %in% names(qc)) qc$max_abs_rho <- NA_real_
  qc <- qc |>
    dplyr::arrange(
      dplyr::desc(.data$n_unique_units),
      dplyr::desc(.data$n_nodes),
      dplyr::desc(.data$median_abs_rho),
      .data$community_id
    )

  rp_map <- qc |>
    dplyr::filter(.data$recurrent) |>
    dplyr::mutate(program_id = paste0("RP", dplyr::row_number()))
  tmp <- dplyr::left_join(tmp, rp_map[, c("community_id", "program_id")], by = "community_id")
  qc <- dplyr::left_join(qc, rp_map[, c("community_id", "program_id")], by = "community_id")

  list(
    nodes = tibble::as_tibble(tmp),
    ledger = tibble::as_tibble(qc),
    recurrent_programs = tibble::as_tibble(rp_map),
    definition = "MCL_community_with_unique_analysis_unit_gate",
    not_definition = "connected_components_or_n_nodes_alone"
  )
}
