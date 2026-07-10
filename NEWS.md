# CellChatAccelRcpp News

## CellChatAccelRcpp 0.1.3

- Confirmed 64-bit R vector indexing in the C++ probability and aggregation kernels for large output tensors.
- Kept `sparse_stream` as the default `computeCommunProbAccelRcpp()` probability kernel.
- Rebuilt the GitHub and Zenodo release archive so the 64-bit sparse-stream implementation is archived as v0.1.3.

## CellChatAccelRcpp 0.1.2

- Added 64-bit R vector indexing in C++ probability and aggregation kernels for large output tensors.
- Rebuilt the GitHub and Zenodo release archive for the long-vector indexing update.
- Kept `sparse_stream` as the default `computeCommunProbAccelRcpp()` probability kernel.

## CellChatAccelRcpp 0.1.1

- Made `sparse_stream` the default `computeCommunProbAccelRcpp()` probability kernel.
- Added the `sparse_stream` communication-probability kernel for active sender-receiver ligand-receptor pairs.
- Added sparse benchmark summaries, checkpoint grids, memory-stress evidence and publication figure outputs.
- Updated the README, package metadata and citation metadata for the v0.1.1 GitHub release.
- Clarified that the package preserves the CellChat model and accelerates execution without removing the sender-by-receiver group-pair scaling term.

## CellChatAccelRcpp 0.1.0

- Initial archived release for manuscript submission.
- Included the R/Rcpp acceleration package, benchmark scripts, processed benchmark tables, Supplementary Table S1, manuscript sources and figure-generation materials.
- Reported 216 paired original/accelerated CellChat comparisons with an 11.4x median speedup and numerical agreement to floating-point precision.
