# metagenealign

`metagenealign` is an R framework for outcome-blind, cross-cohort discovery of recurrent transcriptional programs.

It combines **unit-wise stabilized ICA**, **threshold-free reciprocal-best-hit (RBH) alignment**, **Markov clustering (MCL)**, a **unique-unit recurrence gate**, **unit-balanced consensus metagenes**, and **post-hoc biological annotation**.

> **Framework, not a novel algorithm.** The package assembles established methods into a reproducible workflow with explicit safeguards against outcome leakage, anchor-driven discovery, recurrence miscounting, and one-cohort-overweighting.

## Workflow

```text
Cohort 1 ── stabilized ICA ──┐
Cohort 2 ── stabilized ICA ──┤
Cohort 3 ── stabilized ICA ──┤
                 ...          ├── loading similarity
Cohort N ── stabilized ICA ──┘
                              ↓
                    threshold-free RBH graph
                              ↓
                             MCL
                              ↓
                  recurrent program communities
                              ↓
                 unique-unit recurrence gate
                              ↓
                unit-balanced consensus metagenes
                              ↓
             Hallmark / Reactome / registry annotation
```

## Design principles

- Analyze cohorts independently; do not pool heterogeneous cohorts for discovery.
- Do not use outcome variables during program discovery or alignment.
- Keep discovery features separate from annotation gene sets.
- Define recurrence using **unique analysis units**, never raw node counts.
- Use MCL communities for biological programs; connected components are not program definitions.
- Give each analysis unit one vote when building consensus metagenes.
- Allow zero, one, or multiple post-hoc annotations per recurrent program.
- Manual review fields are never auto-approved.

## Installation

```r
remotes::install_github("xuzhaiyue/metagenealign")
```

## Status

Version `0.1.0` is the first public framework release candidate. GitHub Actions runs R CMD check on the public repository before a release tag is created.
