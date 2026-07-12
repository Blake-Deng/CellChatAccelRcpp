# Experiment Design For The Sparse-Stream Manuscript Benchmark

This benchmark supports the Bioinformatics Application Note for the v0.1.5
`sparse_stream` release of CellChatAccelRcpp. The benchmark is organized around
one public claim: the 64-bit `sparse_stream` implementation preserves CellChat
communication outputs while reducing runtime and memory pressure for repeated
large single-cell RNA workflows.

## Dataset Panel

The primary benchmark uses prepared Seurat RDS inputs from two collections:

1. `3CA_data`, used as Curated Cancer Cell Atlas-derived tumor benchmark inputs.
2. `normal_control`, using PRJNA871268-derived Sample5 and Sample8 controls.

Raw input objects are not committed to GitHub because they are large and
dataset-specific. The repository tracks the accession table, manifest,
experiment grid, processed summaries and source plots needed to reproduce the
analysis when the source objects are available.

## Paired Benchmark

Each dataset is evaluated at six target cell scales:

- 1k
- 5k
- 10k
- 25k
- 50k
- all available cells

Each dataset-scale pair is repeated with three random seeds. For the paired
comparison, original CellChat and CellChatAccelRcpp start from the same prepared
CellChat object and use the same `triMean`, `nboot = 100`, non-spatial RNA
workflow and group labels.

The accelerated branch uses `computeCommunProbAccelRcpp()` with the default
64-bit `sparse_stream` probability kernel.

## Numerical Agreement

Numerical agreement is assessed by comparing flattened CellChat communication
probability tensors between original CellChat and `sparse_stream` outputs.

Primary agreement metrics:

- maximum absolute difference in `cellchat@net$prob`
- Pearson correlation of flattened probability tensors
- ligand-receptor output dimensions and non-zero output consistency

The manuscript reports floating-point agreement across the paired benchmark.

## Runtime And Memory

Runtime is summarized by target cell scale as median and interquartile range.
Speedup is defined as original CellChat elapsed time divided by `sparse_stream`
elapsed time for the same dataset, target cell scale and repeat.

The benchmark also records peak resident memory where available. A separate
large Xenium stress case is reported in the manuscript to demonstrate the
practical memory and runtime barrier of original CellChat under high group
resolution and the completed 64-bit `sparse_stream` run.

## Component Ablation

Component ablations quantify how much of the observed speedup comes from the
probability kernel versus downstream pathway or network aggregation. The
manuscript-facing interpretation is that the communication probability kernel is
the dominant acceleration target, while the accelerated pathway and network
steps preserve compatibility with the standard CellChat workflow.

## Reproducibility Checklist

- Keep package versions in `results/environment/R_packages.txt`.
- Keep the exact experiment grid under the benchmark `results/` directory.
- Keep one metrics file per experiment in `results/runs/`.
- Produce summary tables from `code/summarize_results.R`.
- Keep the current manuscript Figure 1 under `paper/figures/`.
- Do not hand-edit processed benchmark tables.
