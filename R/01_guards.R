# Internal utilities — outcome firewall + guards

.OUTCOME_TOKENS <- c(
  "pcr", "response", "survival", "os", "pfs", "dfs", "event",
  "time_to_event", "outcome", "treatment_response", "residual_cancer_burden",
  "rcb", "auc", "label_y"
)

.canonical_arg_name <- function(x) {
  x <- tolower(as.character(x))
  x <- gsub("[^a-z0-9]+", "_", x)
  gsub("^_+|_+$", "", x)
}

#' Reject outcome / response arguments
#'
#' `metagenealign` is outcome-agnostic by design. This guard inspects argument
#' names, not the provenance of matrices or gene sets; users must still ensure
#' that discovery inputs were created without outcome-driven selection.
#'
#' @param ... Named arguments to inspect.
#' @param extra_forbidden Optional additional forbidden argument names.
#' @return Invisibly `TRUE` if safe; otherwise errors.
#' @export
assert_no_outcome_args <- function(..., extra_forbidden = NULL) {
  dots <- list(...)
  if (!length(dots)) return(invisible(TRUE))
  nms <- names(dots)
  if (is.null(nms)) nms <- rep("", length(dots))
  canon <- .canonical_arg_name(nms[nzchar(nms)])
  forbidden <- unique(.canonical_arg_name(c(.OUTCOME_TOKENS, extra_forbidden)))
  exact_hit <- canon[canon %in% forbidden]
  token_hit <- canon[vapply(canon, function(nm) {
    parts <- strsplit(nm, "_", fixed = TRUE)[[1]]
    any(parts %in% c("pcr", "response", "survival", "outcome", "auc", "rcb")) ||
      grepl("time_to_event|residual_cancer_burden|treatment_response", nm, fixed = FALSE)
  }, logical(1))]
  hit <- unique(c(exact_hit, token_hit))
  if (length(hit)) {
    stop(
      "metagenealign is outcome-agnostic. Refusing argument name(s): ",
      paste(hit, collapse = ", "),
      ". Discovery must not use pCR/response/survival/AUC.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Keep discovery gene universe separate from annotation gene sets
#'
#' @param discovery_genes Character vector of assay-derived discovery features.
#' @param annotation_gene_sets Named list of post-hoc anchors/pathways.
#' @param max_annotation_fraction Soft warning if annotation genes dominate discovery.
#' @return Invisibly `TRUE`.
#' @export
assert_discovery_annotation_separated <- function(
    discovery_genes,
    annotation_gene_sets = NULL,
    max_annotation_fraction = 0.25
) {
  discovery_genes <- unique(as.character(discovery_genes))
  discovery_genes <- discovery_genes[!is.na(discovery_genes) & nzchar(discovery_genes)]
  if (!length(discovery_genes)) stop("discovery_genes is empty.", call. = FALSE)
  if (is.null(annotation_gene_sets)) return(invisible(TRUE))
  if (!is.list(annotation_gene_sets) || is.null(names(annotation_gene_sets)) ||
      any(!nzchar(names(annotation_gene_sets)))) {
    stop("annotation_gene_sets must be a named list with non-empty names.", call. = FALSE)
  }
  ann <- unique(as.character(unlist(annotation_gene_sets, use.names = FALSE)))
  ann <- ann[!is.na(ann) & nzchar(ann)]
  if (setequal(discovery_genes, ann) && length(ann) > 0L) {
    stop(
      "discovery_genes equals the annotation gene-set union. ",
      "Anchors/pathways are post-hoc only and must not define discovery.",
      call. = FALSE
    )
  }
  frac <- length(intersect(discovery_genes, ann)) / length(discovery_genes)
  if (is.finite(frac) && frac > max_annotation_fraction) {
    warning(
      sprintf(
        "Annotation genes cover %.0f%% of discovery_genes (>%.0f%%). ",
        100 * frac, 100 * max_annotation_fraction
      ),
      "Confirm discovery is assay-derived, not registry-restricted.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Validate recurrence bookkeeping
#'
#' @param n_nodes Node count in a community.
#' @param n_unique_units Unique analysis units represented.
#' @param n_units_total Total units in the recurrence denominator.
#' @param recurrence_prop Fraction of units required (alternative to `min_unique_units`).
#' @param min_unique_units Absolute unique-unit gate.
#' @return Logical scalar indicating recurrence.
#' @export
assert_valid_recurrence <- function(
    n_nodes,
    n_unique_units,
    n_units_total,
    recurrence_prop = NULL,
    min_unique_units = NULL
) {
  vals <- as.integer(c(n_nodes, n_unique_units, n_units_total))
  if (any(is.na(vals)) || any(vals < 0L)) stop("Invalid recurrence counts.", call. = FALSE)
  n_nodes <- vals[1]; n_unique_units <- vals[2]; n_units_total <- vals[3]
  if (n_units_total < 1L) stop("n_units_total must be >= 1.", call. = FALSE)
  if (n_unique_units > n_nodes) {
    stop("Invalid recurrence: n_unique_units > n_nodes.", call. = FALSE)
  }
  if (n_unique_units > n_units_total) {
    stop("Invalid recurrence: n_unique_units > n_units_total.", call. = FALSE)
  }
  if (!is.null(min_unique_units)) {
    min_unique_units <- as.integer(min_unique_units)
    if (length(min_unique_units) != 1L || is.na(min_unique_units) || min_unique_units < 1L ||
        min_unique_units > n_units_total) {
      stop("min_unique_units must be between 1 and n_units_total.", call. = FALSE)
    }
    return(n_unique_units >= min_unique_units)
  }
  if (!is.null(recurrence_prop)) {
    recurrence_prop <- as.numeric(recurrence_prop)
    if (length(recurrence_prop) != 1L || !is.finite(recurrence_prop) ||
        recurrence_prop <= 0 || recurrence_prop > 1) {
      stop("recurrence_prop must be in (0, 1].", call. = FALSE)
    }
    return(n_unique_units >= ceiling(n_units_total * recurrence_prop))
  }
  stop("Provide min_unique_units or recurrence_prop.", call. = FALSE)
}

.validate_expression_matrix <- function(x, unit = "matrix") {
  if (!is.matrix(x)) x <- as.matrix(x)
  if (!is.numeric(x)) stop(unit, " must be numeric.", call. = FALSE)
  if (is.null(rownames(x)) || any(!nzchar(rownames(x)))) {
    stop(unit, " must have non-empty gene rownames.", call. = FALSE)
  }
  if (anyDuplicated(rownames(x))) {
    stop(unit, " has duplicated gene rownames; resolve them before discovery.", call. = FALSE)
  }
  if (ncol(x) < 3L) stop(unit, " has too few samples.", call. = FALSE)
  x
}

.stable_string_seed <- function(x, seed = 1L) {
  z <- utf8ToInt(enc2utf8(as.character(x)))
  h <- if (length(z)) sum(z * seq_along(z)) %% 1000000L else 0L
  as.integer(seed) + as.integer(h)
}

.signed_spearman <- function(a, b, min_shared = 50L) {
  shared <- intersect(names(a), names(b))
  if (length(shared) < min_shared) return(NA_real_)
  suppressWarnings(stats::cor(a[shared], b[shared], method = "spearman", use = "complete.obs"))
}

.abs_spearman <- function(a, b, min_shared = 50L) {
  r <- .signed_spearman(a, b, min_shared = min_shared)
  if (!is.finite(r)) NA_real_ else abs(r)
}

.orient_by_max_abs <- function(v) {
  if (!length(v)) return(v)
  w <- which.max(abs(v))
  if (length(w) && is.finite(v[w]) && v[w] < 0) v <- -v
  v
}
