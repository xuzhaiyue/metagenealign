# metagenealign 0.1.0

* Initial public framework for outcome-agnostic cross-cohort latent program discovery.
* Implements stabilized ICA, threshold-free RBH, MCL communities, unique-unit
  recurrence, unit-balanced consensus metagenes, post-hoc annotation and human review.
* Stabilization requires both loading similarity and representation across ICA restarts;
  singleton stability clusters are not treated as perfectly stable.
* k selection can use a stability-coverage score so a solution retaining only a few
  components is not favored solely by high mean stability.
* RBH supports explicit tie handling and honors optional precomputed similarity matrices.
* Recurrence audits optional source-study and platform-class metadata.
* Consensus metagenes enforce one vote per analysis unit and a cross-unit gene-support gate.
* Explicit guards cover outcome-like arguments, discovery/annotation separation,
  unique-unit recurrence, MCL-not-connected-components, and manual review.
