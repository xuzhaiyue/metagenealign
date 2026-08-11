# Human review gate ------------------------------------------------------------

.REVIEW_REQUIRED_COLS <- c(
  "program_id",
  "n_nodes",
  "n_unique_units",
  "n_source_studies",
  "n_platform_classes",
  "max_nodes_per_unit",
  "duplicate_unit_fraction",
  "median_abs_rbh",
  "mcl_stable",
  "leiden_supported",
  "empirical_fdr_supported",
  "candidate_annotation",
  "secondary_annotation",
  "annotation_confidence",
  "reviewed",
  "reviewed_by",
  "reviewed_at",
  "note"
)

#' Create a human review ledger for recurrent programs
#'
#' `reviewed` is always initialized to `FALSE`.
#'
#' @param ledger Program QC ledger from [summarize_program_recurrence()].
#' @param annotations Optional program annotation table.
#' @return Review tibble.
#' @export
create_program_review <- function(ledger, annotations = NULL) {
  ledger <- tibble::as_tibble(ledger)
  if (!"program_id" %in% names(ledger)) {
    stop("ledger must contain program_id.", call. = FALSE)
  }
  ledger <- ledger |> dplyr::filter(!is.na(.data$program_id))
  defaults <- list(
    n_unique_source_studies = NA_integer_,
    n_platform_classes = NA_integer_,
    duplicate_unit_fraction = NA_real_,
    median_abs_rho = NA_real_,
    qc_flag = NA_character_
  )
  for (nm in names(defaults)) if (!nm %in% names(ledger)) ledger[[nm]] <- defaults[[nm]]

  rev <- ledger |>
    dplyr::transmute(
      program_id = .data$program_id,
      n_nodes = .data$n_nodes,
      n_unique_units = .data$n_unique_units,
      n_source_studies = .data$n_unique_source_studies,
      n_platform_classes = .data$n_platform_classes,
      max_nodes_per_unit = .data$max_nodes_per_unit,
      duplicate_unit_fraction = .data$duplicate_unit_fraction,
      median_abs_rbh = .data$median_abs_rho,
      mcl_stable = NA,
      leiden_supported = NA,
      empirical_fdr_supported = NA,
      candidate_annotation = NA_character_,
      secondary_annotation = NA_character_,
      annotation_confidence = NA_character_,
      reviewed = FALSE,
      reviewed_by = NA_character_,
      reviewed_at = as.POSIXct(NA),
      note = .data$qc_flag
    )

  if (!is.null(annotations) && nrow(annotations)) {
    tops <- annotations |>
      dplyr::filter(.data$significant) |>
      dplyr::group_by(.data$program_id) |>
      dplyr::arrange(.data$padj, dplyr::desc(abs(.data$NES)), .by_group = TRUE) |>
      dplyr::summarise(
        candidate_annotation = .data$pathway[1],
        secondary_annotation = if (dplyr::n() >= 2L) .data$pathway[2] else NA_character_,
        .groups = "drop"
      )
    rev <- rev |>
      dplyr::select(-dplyr::any_of(c("candidate_annotation", "secondary_annotation"))) |>
      dplyr::left_join(tops, by = "program_id")
  }
  rev$reviewed <- FALSE
  rev[, .REVIEW_REQUIRED_COLS]
}

#' Validate a filled program review table
#'
#' @param review_tbl Review tibble/data frame.
#' @param require_reviewed If TRUE, every row must be manually reviewed.
#' @return Invisibly the validated table.
#' @export
validate_program_review <- function(review_tbl, require_reviewed = FALSE) {
  missing <- setdiff(.REVIEW_REQUIRED_COLS, names(review_tbl))
  if (length(missing)) {
    stop("Review table missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  if (any(is.na(review_tbl$program_id) | !nzchar(as.character(review_tbl$program_id)))) {
    stop("program_id must be non-missing.", call. = FALSE)
  }
  if (anyDuplicated(review_tbl$program_id)) stop("program_id must be unique.", call. = FALSE)
  if (any(review_tbl$n_unique_units > review_tbl$n_nodes, na.rm = TRUE)) {
    stop("Invalid rows: n_unique_units > n_nodes.", call. = FALSE)
  }
  if (!is.logical(review_tbl$reviewed)) stop("reviewed must be logical.", call. = FALSE)
  if (require_reviewed && !all(!is.na(review_tbl$reviewed) & review_tbl$reviewed)) {
    stop("Not all programs are reviewed=TRUE.", call. = FALSE)
  }
  done <- which(!is.na(review_tbl$reviewed) & review_tbl$reviewed)
  if (length(done)) {
    bad <- done[
      is.na(review_tbl$reviewed_by[done]) |
        !nzchar(as.character(review_tbl$reviewed_by[done])) |
        is.na(review_tbl$reviewed_at[done])
    ]
    if (length(bad)) stop("reviewed=TRUE rows need reviewed_by and reviewed_at.", call. = FALSE)
  }
  invisible(review_tbl)
}
