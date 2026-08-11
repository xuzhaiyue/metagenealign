# Post-hoc biological annotation ----------------------------------------------

.clean_rank_vector <- function(consensus) {
  if (is.null(names(consensus))) stop("consensus must be a named numeric vector.", call. = FALSE)
  x <- as.numeric(consensus)
  names(x) <- names(consensus)
  keep <- is.finite(x) & !is.na(names(x)) & nzchar(names(x))
  x <- x[keep]
  if (anyDuplicated(names(x))) stop("consensus has duplicated gene names.", call. = FALSE)
  if (length(x) < 2L) stop("Too few finite ranked genes.", call. = FALSE)
  sort(x, decreasing = TRUE)
}

#' Preranked fgsea annotation of a consensus metagene
#'
#' Allows zero, one, or many significant gene sets and never forces a label.
#'
#' @param consensus Named numeric loading vector.
#' @param gene_sets Named list of gene sets.
#' @param min_size Minimum overlap size.
#' @param max_size Maximum overlap size; defaults to all ranked genes.
#' @param fdr_cutoff FDR threshold for the `significant` flag.
#' @param ... Reserved; outcome-like argument names are rejected.
#' @return Enrichment tibble, possibly empty.
#' @export
annotate_program_fgsea <- function(
    consensus,
    gene_sets,
    min_size = 5L,
    max_size = Inf,
    fdr_cutoff = 0.05,
    ...
) {
  assert_no_outcome_args(...)
  if (!requireNamespace("fgsea", quietly = TRUE)) {
    stop("Package 'fgsea' is required for annotate_program_fgsea().", call. = FALSE)
  }
  if (!is.list(gene_sets) || is.null(names(gene_sets))) {
    stop("gene_sets must be a named list.", call. = FALSE)
  }
  ranks <- .clean_rank_vector(consensus)
  if (!is.finite(max_size)) max_size <- length(ranks)
  paths <- lapply(gene_sets, function(gs) intersect(unique(as.character(gs)), names(ranks)))
  sizes <- vapply(paths, length, integer(1))
  paths <- paths[sizes >= min_size & sizes <= max_size]
  if (!length(paths)) {
    return(tibble::tibble(
      pathway = character(), NES = double(), pval = double(), padj = double(),
      size = integer(), significant = logical(), leadingEdge = character()
    ))
  }

  fg <- tryCatch(
    fgsea::fgsea(
      pathways = paths,
      stats = ranks,
      minSize = as.integer(min_size),
      maxSize = as.integer(max_size)
    ),
    error = function(e) NULL
  )
  if (is.null(fg) || !nrow(fg)) {
    return(tibble::tibble(
      pathway = character(), NES = double(), pval = double(), padj = double(),
      size = integer(), significant = logical(), leadingEdge = character()
    ))
  }

  tibble::as_tibble(fg) |>
    dplyr::transmute(
      pathway = .data$pathway,
      NES = .data$NES,
      pval = .data$pval,
      padj = .data$padj,
      size = .data$size,
      significant = is.finite(.data$padj) & .data$padj <= fdr_cutoff,
      leadingEdge = vapply(.data$leadingEdge, function(x) paste(x, collapse = ";"), character(1))
    ) |>
    dplyr::arrange(.data$padj, dplyr::desc(abs(.data$NES)))
}

#' Map recurrent programs to a prespecified registry
#'
#' This function is post-hoc annotation only. Separation between the discovery
#' universe and annotation registry must be enforced upstream.
#'
#' @param consensus_list Named program consensus vectors.
#' @param registry_gene_sets Named registry gene sets.
#' @param fdr_cutoff FDR cutoff.
#' @param ... Forwarded to [annotate_program_fgsea()].
#' @return Long program-by-registry enrichment table.
#' @export
map_programs_to_registry <- function(
    consensus_list,
    registry_gene_sets,
    fdr_cutoff = 0.05,
    ...
) {
  assert_no_outcome_args(...)
  if (is.null(names(consensus_list)) || any(!nzchar(names(consensus_list)))) {
    stop("consensus_list must have program IDs as names.", call. = FALSE)
  }
  dplyr::bind_rows(lapply(names(consensus_list), function(pid) {
    annotate_program_fgsea(
      consensus_list[[pid]], registry_gene_sets,
      fdr_cutoff = fdr_cutoff, ...
    ) |>
      dplyr::mutate(program_id = pid, .before = 1)
  }))
}

#' Summarize registry-state coverage by recurrent programs
#'
#' @param mapping_tbl Output of [map_programs_to_registry()].
#' @param program_units Named list or data frame linking programs to units.
#' @param registry_order Optional state order; use this to retain honest zero rows.
#' @param n_units_total Total units.
#' @param min_unique_units Recovery gate on unique units covered.
#' @return Coverage ledger.
#' @export
summarize_anchor_recovery <- function(
    mapping_tbl,
    program_units,
    registry_order = NULL,
    n_units_total,
    min_unique_units = 6L
) {
  n_units_total <- as.integer(n_units_total)
  if (n_units_total < 1L) stop("n_units_total must be >=1.", call. = FALSE)
  if (is.data.frame(program_units)) {
    if (!all(c("program_id", "analysis_unit") %in% names(program_units))) {
      stop("program_units data frame needs program_id and analysis_unit.", call. = FALSE)
    }
    pu <- dplyr::distinct(program_units, .data$program_id, .data$analysis_unit)
  } else {
    pu <- dplyr::bind_rows(lapply(names(program_units), function(pid) {
      tibble::tibble(program_id = pid, analysis_unit = as.character(program_units[[pid]]))
    }))
  }
  states <- if (!is.null(registry_order)) registry_order else sort(unique(mapping_tbl$pathway))
  dplyr::bind_rows(lapply(states, function(st) {
    progs_tbl <- mapping_tbl |>
      dplyr::filter(.data$pathway == st, .data$significant) |>
      dplyr::distinct(.data$program_id)
    progs <- progs_tbl$program_id
    units_tbl <- pu |>
      dplyr::filter(.data$program_id %in% progs) |>
      dplyr::distinct(.data$analysis_unit)
    units_hit <- units_tbl$analysis_unit
    tibble::tibble(
      state_id = st,
      n_mapped_programs = length(progs),
      n_units_covered = length(units_hit),
      unit_fraction = length(units_hit) / n_units_total,
      recovered = length(units_hit) >= as.integer(min_unique_units),
      programs = paste(sort(progs), collapse = ";"),
      units = paste(sort(units_hit), collapse = ";")
    )
  }))
}
