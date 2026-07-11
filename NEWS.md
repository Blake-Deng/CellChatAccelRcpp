# CellChatAccelRcpp News

## CellChatAccelRcpp 0.1.4 - current public release

- Keeps the 64-bit `sparse_stream` probability kernel as the default public workflow.
- Aligns repository-facing benchmark summaries with the current manuscript figure and tables, including the 216-run paired benchmark, memory summaries and the large Xenium stress test.
- Adds reviewer-facing package hygiene checks, including source build/check instructions and a smoke test for exported functions and 64-bit R sessions.
- Cleans unused C++ helper code to avoid install-time compiler warnings.

## CellChatAccelRcpp 0.1.3

- Confirmed 64-bit R vector indexing in the C++ probability and aggregation kernels for large output tensors.
- Kept `sparse_stream` as the default `computeCommunProbAccelRcpp()` probability kernel.
- Archived the 64-bit `sparse_stream` implementation as the current GitHub and Zenodo software release.

The public repository and manuscript materials are organized around the v0.1.4 `sparse_stream` release. Internal iteration notes are not part of the current manuscript narrative.
