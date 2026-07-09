# CellChatAccelRcpp News

## CellChatAccelRcpp 0.1.1

- Made `sparse_stream` the default `computeCommunProbAccelRcpp()` probability kernel.
- Added the `sparse_stream` communication-probability kernel for active sender-receiver ligand-receptor pairs.
- Added sparse benchmark summaries, checkpoint grids, memory-stress evidence and publication figure outputs.
- Added 64-bit R vector indexing in C++ probability and aggregation kernels for large output tensors.
- Updated the README, package metadata and citation metadata for the v0.1.1 GitHub release.
- Clarified that the package preserves the CellChat model and accelerates execution without removing the sender-by-receiver group-pair scaling term.

## CellChatAccelRcpp 0.1.0

- Initial archived release for manuscript submission.
- Included the R/Rcpp acceleration package, benchmark scripts, processed benchmark tables, Supplementary Table S1, manuscript sources and figure-generation materials.
- Reported 216 paired original/accelerated CellChat comparisons with an 11.4x median speedup and numerical agreement to floating-point precision.
