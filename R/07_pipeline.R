# High-level pipeline ----------------------------------------------------------

#' End-to-end cross-cohort recurrent program discovery
#'
#' Each cohort/unit is decomposed independently. Cross-cohort alignment uses
#' threshold-free RBH and MCL. Recurrence is defined by unique analysis units;
#' consensus metagenes give one vote per unit. Annotation is post-hoc only.
#'
#' @param matrices Named list of genes x samples numeric matrices.
#' @param discovery_genes Assay-derived eligible discovery gene universe.
#' @param annotation_gene_sets Optional post-hoc registry/pathway list.
#' @param unit_metadata Optional data frame with `analysis_unit` and optional
#'   `source_study`, `platform_class`. Extra columns are ignored.
#' @param k_grid ICA k candidates.
#' @param n_runs Stabilized ICA restarts.
#' @param min_stability Component stability gate.
#' @param min_run_fraction Restart-coverage gate for stabilized components.
#' @param k_selection_metric k selection rule passed to [scan_ica_k()].
#' @param max_genes_per_unit Optional within-unit top-variance cap for ICA
#'   fitting. The full `discovery_genes` vector remains the eligibility universe.
#' @param min_shared_genes Minimum loading overlap for RBH.
#' @param mcl_inflation MCL inflation.
#' @param min_unique_units Unique-unit recurrence gate.
#' @param min_gene_unit_fraction Cross-unit support gate for consensus genes.
#' @param strict_units If TRUE, any unit that cannot be analyzed causes an error.
#' @param recurrence_denominator `"analyzed"` (default) or `"input"`.
#' @param seed Random seed.
#' @param ... Reserved; outcome-like argument names are rejected.
#' @return Structured discovery result.
#' @export
discover_recurrent_programs <- function(
    matrices,
    discovery_genes,
    annotation_gene_sets = NULL,
    unit_metadata = NULL,
    k_grid = c(10L, 15L, 20L, 25L, 30L),
    n_runs = 50L,
    min_stability = 0.55,
    min_run_fraction = 0.80,
    k_selection_metric = c("stability_coverage", "mean_stability"),
    max_genes_per_unit = 3000L,
    min_shared_genes = 200L,
    mcl_inflation = 2,
    min_unique_units = 6L,
    min_gene_unit_fraction = 0.50,
    strict_units = TRUE,
    recurrence_denominator = c("analyzed", "input"),
    seed = 1L,
    ...
) {
  assert_no_outcome_args(...)
  assert_discovery_annotation_separated(discovery_genes, annotation_gene_sets)
  k_selection_metric <- match.arg(k_selection_metric)
  recurrence_denominator <- match.arg(recurrence_denominator)
  if (length(strict_units) != 1L || is.na(strict_units) || !is.logical(strict_units)) {
    stop("strict_units must be one logical value.", call. = FALSE)
  }
  if (!is.null(max_genes_per_unit)) {
    max_genes_per_unit <- as.integer(max_genes_per_unit)
    if (length(max_genes_per_unit) != 1L || is.na(max_genes_per_unit) || max_genes_per_unit < 100L) {
      stop("max_genes_per_unit must be NULL or an integer >=100.", call. = FALSE)
    }
  }

  if (!is.list(matrices) || is.null(names(matrices)) || any(!nzchar(names(matrices))) ||
      anyDuplicated(names(matrices))) {
    stop("matrices must be a named list with unique analysis-unit names.", call. = FALSE)
  }
  discovery_genes <- unique(as.character(discovery_genes))
  discovery_genes <- discovery_genes[!is.na(discovery_genes) & nzchar(discovery_genes)]

  meta_units <- tibble::tibble(analysis_unit = names(matrices))
  if (!is.null(unit_metadata)) {
    unit_metadata <- as.data.frame(unit_metadata)
    if (!"analysis_unit" %in% names(unit_metadata)) stop("unit_metadata needs analysis_unit.", call. = FALSE)
    if (anyDuplicated(unit_metadata$analysis_unit)) stop("unit_metadata analysis_unit must be unique.", call. = FALSE)
    keep_meta <- intersect(c("analysis_unit", "source_study", "platform_class"), names(unit_metadata))
    meta_units <- dplyr::left_join(meta_units, unit_metadata[, keep_meta, drop = FALSE], by = "analysis_unit")
  }

  loadings <- list()
  node_meta <- list()
  k_scan <- list()
  unit_audit <- list()
  analyzed_units <- character()

  for (unit in names(matrices)) {
    mat <- .validate_expression_matrix(matrices[[unit]], paste0("matrix '", unit, "'"))
    use <- intersect(discovery_genes, rownames(mat))
    if (length(use) < 100L) {
      msg <- paste0("Unit ", unit, " has too few eligible discovery genes: ", length(use))
      if (strict_units) stop(msg, call. = FALSE) else {
        warning(msg, call. = FALSE)
        unit_audit[[unit]] <- tibble::tibble(
          analysis_unit = unit, analyzed = FALSE, reason = "too_few_discovery_genes",
          n_eligible_genes = length(use), n_fit_genes = NA_integer_, n_samples = ncol(mat)
        )
        next
      }
    }

    sub <- mat[use, , drop = FALSE]
    sds <- matrixStats::rowSds(sub, na.rm = TRUE)
    sub <- sub[is.finite(sds) & sds > 0, , drop = FALSE]
    n_eligible_variable <- nrow(sub)
    if (!is.null(max_genes_per_unit) && is.finite(max_genes_per_unit) &&
        nrow(sub) > as.integer(max_genes_per_unit)) {
      v <- matrixStats::rowVars(sub, na.rm = TRUE)
      ord <- order(v, decreasing = TRUE, na.last = NA)
      sub <- sub[ord[seq_len(min(length(ord), as.integer(max_genes_per_unit)))], , drop = FALSE]
    }

    unit_seed <- .stable_string_seed(unit, seed)
    scanned <- tryCatch(
      scan_ica_k(
        sub,
        k_grid = k_grid,
        n_runs = n_runs,
        seed = unit_seed,
        min_stability = min_stability,
        min_run_fraction = min_run_fraction,
        selection_metric = k_selection_metric
      ),
      error = function(e) e
    )
    if (inherits(scanned, "error")) {
      if (strict_units) stop("Unit ", unit, " failed ICA scan: ", conditionMessage(scanned), call. = FALSE)
      warning("Skipping ", unit, ": ", conditionMessage(scanned), call. = FALSE)
      unit_audit[[unit]] <- tibble::tibble(
        analysis_unit = unit, analyzed = FALSE, reason = "ica_scan_failed",
        n_eligible_genes = n_eligible_variable, n_fit_genes = nrow(sub), n_samples = ncol(sub)
      )
      next
    }

    fit <- scanned$best
    analyzed_units <- c(analyzed_units, unit)
    k_scan[[unit]] <- scanned$scan |>
      dplyr::mutate(analysis_unit = unit, .before = 1)
    unit_audit[[unit]] <- tibble::tibble(
      analysis_unit = unit,
      analyzed = TRUE,
      reason = NA_character_,
      n_eligible_genes = n_eligible_variable,
      n_fit_genes = nrow(sub),
      n_samples = ncol(sub),
      k_selected = fit$k_requested,
      n_components_kept = fit$n_kept,
      mean_stability = fit$mean_stability,
      mean_run_fraction = fit$mean_run_fraction,
      n_runs_ok = fit$n_runs_ok
    )

    unit_meta_row <- meta_units[meta_units$analysis_unit == unit, , drop = FALSE]
    for (j in seq_len(ncol(fit$loadings))) {
      nid <- paste0(unit, "::", colnames(fit$loadings)[j])
      v <- fit$loadings[, j]
      names(v) <- rownames(fit$loadings)
      loadings[[nid]] <- v
      cm <- fit$component_meta[fit$component_meta$component == colnames(fit$loadings)[j], , drop = FALSE]
      row <- tibble::tibble(
        node_id = nid,
        analysis_unit = unit,
        component = colnames(fit$loadings)[j],
        k_selected = fit$k_requested,
        stability = cm$stability,
        run_fraction = cm$run_fraction,
        n_samples = ncol(sub),
        n_genes = nrow(sub)
      )
      if ("source_study" %in% names(unit_meta_row)) row$source_study <- as.character(unit_meta_row$source_study)
      if ("platform_class" %in% names(unit_meta_row)) row$platform_class <- as.character(unit_meta_row$platform_class)
      node_meta[[length(node_meta) + 1L]] <- row
    }
  }

  meta <- dplyr::bind_rows(node_meta)
  if (!nrow(meta)) stop("No components produced.", call. = FALSE)
  n_denominator <- if (recurrence_denominator == "analyzed") length(analyzed_units) else length(matrices)
  if (min_unique_units > n_denominator) {
    stop("min_unique_units exceeds the recurrence denominator.", call. = FALSE)
  }

  unit_of <- stats::setNames(meta$analysis_unit, meta$node_id)
  rbh <- build_rbh_graph(loadings, unit_of = unit_of, min_shared_genes = min_shared_genes)
  memb <- cluster_rbh_mcl(rbh$edges, node_ids = meta$node_id, inflation = mcl_inflation)
  prog <- summarize_program_recurrence(
    memb,
    node_meta = meta,
    edges = rbh$edges,
    min_unique_units = min_unique_units,
    n_units_total = n_denominator
  )

  consensus_list <- list()
  consensus_meta <- list()
  consensus_gene_support <- list()
  if (nrow(prog$recurrent_programs)) {
    for (i in seq_len(nrow(prog$recurrent_programs))) {
      pid <- prog$recurrent_programs$program_id[i]
      member_tbl <- prog$nodes |>
        dplyr::filter(.data$program_id == pid)
      members <- member_tbl$node_id
      cm <- build_consensus_metagene(
        loadings,
        members,
        unit_of = unit_of,
        min_gene_unit_fraction = min_gene_unit_fraction
      )
      consensus_list[[pid]] <- cm$consensus
      consensus_gene_support[[pid]] <- cm$gene_support |>
        dplyr::mutate(program_id = pid, .before = 1)
      consensus_meta[[length(consensus_meta) + 1L]] <- tibble::tibble(
        program_id = pid,
        medoid = cm$medoid,
        n_member_nodes = cm$n_member_nodes,
        n_units_voting = cm$n_units_voting,
        n_consensus_genes = length(cm$consensus),
        min_gene_unit_fraction = cm$min_gene_unit_fraction
      )
    }
  }

  annot <- NULL
  recovery <- NULL
  if (!is.null(annotation_gene_sets) && length(consensus_list)) {
    annot <- map_programs_to_registry(consensus_list, annotation_gene_sets)
    pu <- prog$nodes |>
      dplyr::filter(!is.na(.data$program_id)) |>
      dplyr::distinct(.data$program_id, .data$analysis_unit)
    recovery <- summarize_anchor_recovery(
      annot,
      program_units = pu,
      registry_order = names(annotation_gene_sets),
      n_units_total = n_denominator,
      min_unique_units = min_unique_units
    )
  }

  recurrent_ledger <- prog$ledger |> dplyr::filter(!is.na(.data$program_id))
  review <- create_program_review(recurrent_ledger, annot)

  list(
    components = meta,
    unit_audit = dplyr::bind_rows(unit_audit),
    k_scan = dplyr::bind_rows(k_scan),
    loadings = loadings,
    rbh_edges = rbh$edges,
    membership = memb,
    programs = prog,
    consensus = consensus_list,
    consensus_meta = dplyr::bind_rows(consensus_meta),
    consensus_gene_support = dplyr::bind_rows(consensus_gene_support),
    annotations = annot,
    registry_coverage = recovery,
    review = review,
    method = list(
      within_cohort = "stabilized_ICA_with_restart_coverage",
      k_selection = k_selection_metric,
      cross_cohort = "threshold_free_RBH",
      communities = "MCL_not_connected_components",
      recurrence = "unique_analysis_units",
      recurrence_denominator = recurrence_denominator,
      consensus = "sign_oriented_unit_balanced_with_gene_support_gate",
      outcome = "agnostic_argument_firewall"
    )
  )
}
