test_that("outcome-like arguments are rejected robustly", {
  expect_error(assert_no_outcome_args(pCR = 1), "outcome-agnostic")
  expect_error(assert_no_outcome_args(pcr_status = 1), "outcome-agnostic")
  expect_error(assert_no_outcome_args(clinical_response = TRUE), "outcome-agnostic")
  expect_error(assert_no_outcome_args(AUC = 0.8), "outcome-agnostic")
  expect_silent(assert_no_outcome_args(seed = 1, n_runs = 50))
})

test_that("discovery and annotation gene sets stay separated", {
  disc <- paste0("G", 1:100)
  ann <- list(STATE_A = paste0("G", 1:100))
  expect_error(assert_discovery_annotation_separated(disc, ann), "equals the annotation")
  ann2 <- list(STATE_A = paste0("G", 1:10))
  expect_warning(assert_discovery_annotation_separated(disc, ann2, max_annotation_fraction = 0.05))
  expect_silent(assert_discovery_annotation_separated(disc, ann2, max_annotation_fraction = 0.25))
})

test_that("recurrence uses unique units not node count", {
  expect_true(assert_valid_recurrence(9, 9, 12, min_unique_units = 6))
  expect_false(assert_valid_recurrence(20, 3, 12, min_unique_units = 6))
  expect_error(assert_valid_recurrence(5, 8, 12, min_unique_units = 6))
  expect_true(assert_valid_recurrence(10, 6, 12, recurrence_prop = 0.5))
  expect_false(assert_valid_recurrence(10, 5, 12, recurrence_prop = 0.5))
})

test_that("unit-balanced consensus gives one vote per unit", {
  set.seed(1)
  genes <- paste0("G", 1:50)
  mk <- function(s) stats::setNames(rnorm(50, mean = s), genes)
  loadings <- list(U1_a = mk(1), U1_b = mk(1.1), U2_a = mk(1), U3_a = mk(1))
  unit_of <- c(U1_a = "U1", U1_b = "U1", U2_a = "U2", U3_a = "U3")
  cm <- build_consensus_metagene(
    loadings, names(loadings), unit_of,
    min_gene_unit_fraction = 0.5
  )
  expect_equal(cm$n_member_nodes, 4L)
  expect_equal(cm$n_units_voting, 3L)
  expect_true(all(cm$gene_support$n_units <= 3L))
})

test_that("consensus gene support excludes one-unit-only genes", {
  loadings <- list(
    A = c(G1 = 1, G2 = 1, ONLY_A = 5),
    B = c(G1 = 1.1, G2 = 1.2),
    C = c(G1 = 0.9, G2 = 1.1)
  )
  unit_of <- c(A = "U1", B = "U2", C = "U3")
  cm <- build_consensus_metagene(loadings, names(loadings), unit_of, min_gene_unit_fraction = 2/3)
  expect_false("ONLY_A" %in% names(cm$consensus))
})

test_that("MCL attribute records not connected components", {
  nodes <- paste0("N", 1:6)
  edges <- tibble::tibble(
    from = c("N1", "N1", "N2", "N4", "N4", "N5"),
    to   = c("N2", "N3", "N3", "N5", "N6", "N6"),
    abs_spearman = c(0.8, 0.7, 0.75, 0.8, 0.7, 0.75)
  )
  cl <- cluster_rbh_mcl(edges, nodes, inflation = 2)
  expect_equal(attr(cl, "method"), "Markov_Clustering_MCL")
  expect_equal(attr(cl, "not_used"), "igraph_connected_components")
  expect_equal(length(cl), 6L)
})

test_that("RBH optional similarity is honored and ties can be retained", {
  loadings <- list(
    A1 = c(G1 = 1, G2 = 2, G3 = 3),
    B1 = c(G1 = 1, G2 = 2, G3 = 3),
    B2 = c(G1 = 3, G2 = 2, G3 = 1)
  )
  units <- c(A1 = "A", B1 = "B", B2 = "B")
  S <- matrix(c(
    1, 0.9, 0.9,
    0.9, 1, 0.1,
    0.9, 0.1, 1
  ), 3, 3, byrow = TRUE, dimnames = list(names(loadings), names(loadings)))
  r <- build_rbh_graph(loadings, units, similarity = S, tie_policy = "all", min_shared_genes = 2)
  expect_equal(nrow(r$edges), 2L)
})

test_that("review never auto-approves", {
  ledger <- tibble::tibble(
    program_id = "RP1",
    n_nodes = 9L,
    n_unique_units = 9L,
    max_nodes_per_unit = 1L,
    duplicate_unit_fraction = 0,
    n_unique_source_studies = 8L,
    n_platform_classes = 3L,
    median_abs_rho = 0.5,
    qc_flag = "OK"
  )
  rev <- create_program_review(ledger)
  expect_true(all(rev$reviewed == FALSE))
  expect_error(validate_program_review(rev, require_reviewed = TRUE))
  rev$reviewed[1] <- TRUE
  rev$reviewed_by[1] <- "expert"
  rev$reviewed_at[1] <- Sys.time()
  expect_silent(validate_program_review(rev, require_reviewed = TRUE))
})

test_that("summarize_program_recurrence audits source studies and platforms", {
  meta <- tibble::tibble(
    node_id = paste0("U", 1:6, "_1"),
    analysis_unit = paste0("U", 1:6),
    source_study = c("S1", "S1", "S2", "S3", "S4", "S5"),
    platform_class = c("array", "array", "rna", "rna", "nano", "array")
  )
  memb <- stats::setNames(rep(1L, 6), meta$node_id)
  edges <- tibble::tibble(from = character(), to = character(), abs_spearman = double())
  out <- summarize_program_recurrence(memb, meta, edges, min_unique_units = 6L, n_units_total = 6L)
  expect_true(out$ledger$recurrent[1])
  expect_equal(out$ledger$n_unique_source_studies[1], 5L)
  expect_equal(out$ledger$n_platform_classes[1], 3L)
})

test_that("node count alone never creates recurrence", {
  meta <- tibble::tibble(
    node_id = paste0("U", rep(1:3, each = 4), "_", 1:12),
    analysis_unit = paste0("U", rep(1:3, each = 4))
  )
  memb <- stats::setNames(rep(1L, 12), meta$node_id)
  edges <- tibble::tibble(from = character(), to = character(), abs_spearman = double())
  out <- summarize_program_recurrence(memb, meta, edges, min_unique_units = 6L, n_units_total = 12L)
  expect_equal(nrow(out$recurrent_programs), 0L)
  expect_true(all(!out$ledger$recurrent))
  expect_equal(out$not_definition, "connected_components_or_n_nodes_alone")
})

test_that("singleton stabilization clusters are not called perfectly stable", {
  big <- matrix(rnorm(30), nrow = 10, ncol = 3)
  C <- diag(3)
  cluster <- c(1L, 2L, 2L)
  run_of <- c(1L, 1L, 2L)
  z <- metagenealign:::.cluster_stability_summary(big, C, cluster, run_of, n_runs_total = 3L)
  one <- z[[which(vapply(z, function(x) x$cluster_id, integer(1)) == 1L)]]
  expect_equal(one$stability, 0)
  expect_equal(one$run_fraction, 1/3)
})

test_that("stability balances duplicate estimates within one restart", {
  big <- matrix(rnorm(40), nrow = 10, ncol = 4)
  C <- matrix(c(
    1, .99, .8, .8,
    .99, 1, .79, .79,
    .8, .79, 1, .9,
    .8, .79, .9, 1
  ), 4, 4, byrow = TRUE)
  cluster <- rep(1L, 4)
  run_of <- c(1L, 1L, 2L, 3L)
  z <- metagenealign:::.cluster_stability_summary(big, C, cluster, run_of, n_runs_total = 3L)[[1]]
  expect_equal(z$n_estimates_raw, 4L)
  expect_equal(z$n_estimates_balanced, 3L)
  expect_equal(z$n_runs_represented, 3L)
  expect_equal(z$run_fraction, 1)
})
