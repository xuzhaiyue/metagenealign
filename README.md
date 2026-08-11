# metagenealign

> **Discover transcriptional programs independently within cohorts and retain only those that reproducibly re-emerge across cohorts.**

`metagenealign` is an **outcome-agnostic framework / workflow** for multi-cohort transcriptomic latent-program discovery. It combines established methods; it is **not presented as a novel matrix-factorization algorithm**.

```text
          Cohort 1 ── stabilized ICA ──┐
          Cohort 2 ── stabilized ICA ──┤
          Cohort 3 ── stabilized ICA ──┤
                         ...            ├─ threshold-free RBH graph
          Cohort N ── stabilized ICA ──┘
                                        ↓
                                       MCL
                                        ↓
                           recurrent metagene programs
                                        ↓
                    unit-balanced consensus gene loadings
                                        ↓
                    Hallmark / Reactome / registry (post hoc)
```

## Why this package exists

The framework encodes scientific failure modes that are easy to introduce in cross-cohort discovery:

| Guard | Rule |
|---|---|
| Stabilized ICA | a component needs both high loading similarity **and recurrence across ICA restarts** |
| Recurrence | use `n_distinct(analysis_unit)`, never node count alone |
| Communities | use **MCL on RBH**; connected components are not biological programs |
| Consensus | component → within-unit collapse → cross-unit consensus; one vote per unit |
| Gene support | consensus genes must be represented across a prespecified fraction of voting units |
| Discovery vs annotation | assay-derived discovery universe is separate from post-hoc gene-set annotation |
| Outcomes | high-level APIs reject outcome-like argument names |
| Review | `reviewed = FALSE` until a human reviews the program ledger |

> The outcome firewall can inspect software arguments, not scientific provenance. Users remain responsible for ensuring matrices, feature-selection rules and gene sets were not constructed using outcomes.

## Install

```r
# GitHub
remotes::install_github("xuzhaiyue/metagenealign")
```

`fgsea` is optional and is needed only for post-hoc gene-set annotation.

## Minimal workflow

```r
library(metagenealign)

# matrices: named list, one genes x samples matrix per independent analysis unit
# discovery_genes: assay-derived eligibility universe, not the anchor registry
# unit_metadata: optional source-study/platform audit information

res <- discover_recurrent_programs(
  matrices = cohort_matrices,
  discovery_genes = discovery_genes,
  annotation_gene_sets = registry_sets,  # post hoc only
  unit_metadata = unit_metadata,
  k_grid = c(15, 20, 25, 30),
  n_runs = 50,
  min_run_fraction = 0.80,
  max_genes_per_unit = 3000,
  min_unique_units = 6,
  min_gene_unit_fraction = 0.50
)

res$unit_audit
res$programs$ledger
res$consensus
res$registry_coverage
res$review  # reviewed is always FALSE initially
```

### Important distinction: eligible universe vs ICA fitting features

`discovery_genes` defines the **eligible outcome-blind discovery universe**. If `max_genes_per_unit` is set, ICA is fitted to the top-variable subset selected independently within each unit. The audit records both counts. This distinction should be reported explicitly in manuscripts.

## Layered API

**Within cohort**

- `run_stabilized_ica()`
- `scan_ica_k()`
- `summarize_component_stability()`
- `select_stable_components()`

**Cross cohort**

- `calculate_loading_similarity()`
- `build_rbh_graph()`
- `cluster_rbh_mcl()`
- `summarize_program_recurrence()`

**Consensus**

- `orient_component_signs()`
- `collapse_components_within_unit()`
- `build_consensus_metagene()`

**Annotation**

- `annotate_program_fgsea()`
- `map_programs_to_registry()`
- `summarize_anchor_recovery()`

**Review**

- `create_program_review()`
- `validate_program_review()`

## Interpretation

Programs are **recurrent latent transcriptional programs**, not discrete cell states. Failure of a registry state to map to a recurrent program means that it did not emerge as an independently recurrent latent axis under the frozen discovery rules; it does **not** imply that the biology is absent or unmeasurable by a prespecified score.

## Methodological provenance

The cross-dataset reproducibility framing is inspired by:

- Cantini L, Kairov U, de Reyniès A, Barillot E, Radvanyi F, Zinovyev A. *Assessing reproducibility of matrix factorization methods in independent transcriptomes.* **Bioinformatics** 2019;35:4307–4313. doi:10.1093/bioinformatics/btz225.
- Kairov U, Cantini L, Greco A, et al. *Determining the optimal number of independent components for reproducible transcriptomic data analysis.* **BMC Genomics** 2017;18:712. doi:10.1186/s12864-017-4112-9.

`metagenealign` is an engineering and review-gated workflow around established components (FastICA, RBH and MCL), not a claim of a novel underlying algorithm.

## Relationship to companion utilities

```text
scpipeline.utils       → single-cell workflows
spatialpipeline.utils  → spatial workflows
metagenealign          → multi-cohort latent-program discovery
```

## License

MIT
