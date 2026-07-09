# CellChatAccelRcpp: scalable Rcpp acceleration of CellChat communication inference

> Draft type: Bioinformatics Application Note  
> Target length: <= 4 printed pages; approximately 2000 words plus one main figure.  
> Current evidence base: 12 datasets, 6 cell scales, 3 repeats, 864 benchmark tasks, 1296 metric rows, 0 failed runs.

## Abstract

### Summary
CellChat is widely used to infer cell-cell communication from single-cell transcriptomic data, but its communication probability estimation and downstream aggregation steps become computational bottlenecks as datasets increase in size. We present CellChatAccelRcpp, an Rcpp-based acceleration layer that preserves the CellChat object interface while replacing selected hot paths with compiled implementations. CellChatAccelRcpp accelerates communication probability estimation, pathway-level aggregation and network aggregation, and supports checkpointed benchmark execution for large-scale studies. Across 12 single-cell datasets, six cell-scale settings and three random repeats, CellChatAccelRcpp achieved a median 11.4-fold speedup over the original CellChat workflow while preserving numerically equivalent communication probability tensors. In 216 paired baseline/accelerated comparisons, the maximum absolute probability difference was 1.39e-16 and Pearson correlation was 1.0 in all comparisons.

### Availability and implementation
CellChatAccelRcpp is implemented in R and C++ via Rcpp. Source code, reproducible benchmark scripts, conda environment files, processed benchmark tables and figure-generation scripts are available at https://github.com/Blake-Deng/CellChatAccelRcpp. The archived v0.1.0 release is available at https://doi.org/10.5281/zenodo.21186108.

### Contact
zifengd8@gmail.com.

### Supplementary information
Supplementary data and figures are available in the project repository and submission package.

## 1. Introduction
Cell-cell communication analysis is a central step in interpreting single-cell and spatial transcriptomic datasets. CellChat provides a widely adopted statistical framework and visualization ecosystem for inferring ligand-receptor mediated interactions among annotated cell groups. However, modern atlases increasingly contain tens of thousands to millions of cells, and the computational cost of repeated communication probability estimation can limit routine use in large-scale benchmarking, perturbation analysis and exploratory workflows.

CellChatAccelRcpp addresses this computational bottleneck by accelerating selected CellChat hot paths while keeping the familiar CellChat workflow intact. Rather than introducing a new biological model, the package aims to provide a drop-in acceleration layer for the same probability tensor and downstream network summaries used by CellChat. This design makes it suitable for users who require larger-scale analyses but need outputs that remain directly comparable with existing CellChat results.

## 2. Materials and methods

### Implementation
CellChatAccelRcpp implements compiled Rcpp routines for triMean-based group expression aggregation, ligand-receptor communication probability computation, pathway-level communication aggregation and network aggregation. The accelerated functions operate on standard CellChat objects and return CellChat-compatible outputs. The package currently targets RNA-based CellChat workflows and prioritizes numerical equivalence to the original CellChat implementation.

### Benchmark design
We benchmarked CellChatAccelRcpp against the original CellChat workflow using 12 real single-cell datasets stored in the benchmark workspace. Each dataset was evaluated at six target cell scales: 1k, 5k, 10k, 25k, 50k and all available cells. Each dataset-scale pair was repeated with three random seeds. The main paired benchmark ran both original CellChat and CellChatAccelRcpp from the same prepared CellChat object, enabling direct runtime and numerical equivalence comparisons. Additional component ablations disabled the accelerated probability kernel, pathway aggregation or network aggregation to quantify the contribution of each component.

All experiments were run in a conda-managed R environment on a Linux server with 256 CPU threads and approximately 1 TiB RAM. The benchmark runner used checkpointing at the prepared CellChat object and computed-result levels, enabling resumable execution and protecting against partial failures.

### Metrics
Runtime was measured as elapsed wall-clock time for each engine-specific inference step. Numerical equivalence was assessed by comparing flattened CellChat communication probability tensors between original and accelerated outputs using maximum absolute difference and Pearson correlation. Component ablation was summarized as runtime relative to the full accelerated path.

## 3. Results

CellChatAccelRcpp completed all 864 benchmark tasks without failed metric files. In 216 paired original/accelerated comparisons, the accelerated workflow achieved a median 11.4-fold speedup overall. Median speedup was strongest at smaller scales, reaching 37.4-fold at 1k cells, and remained substantial at larger settings, with approximately 6.0-fold median speedup at the all-cell scale. Median elapsed time decreased from 426.6 seconds for the original CellChat workflow to 36.0 seconds for CellChatAccelRcpp across all paired benchmarks.

The accelerated results were numerically matched to original CellChat outputs. Across all 216 paired comparisons, the maximum absolute probability difference was 1.39e-16 and the minimum Pearson correlation was 1.0. These results indicate that the acceleration preserves the communication probability tensor up to floating-point precision.

Ablation analysis showed that the compiled communication probability kernel contributed most of the runtime reduction. Disabling the accelerated kernel substantially increased runtime relative to the full accelerated path, whereas disabling pathway or network aggregation produced smaller slowdowns. This suggests that probability estimation is the dominant optimization target for large CellChat workflows.

## 4. Discussion
CellChatAccelRcpp provides a practical acceleration layer for CellChat analyses on large single-cell datasets. By retaining CellChat-compatible outputs and matching original probability tensors numerically, the package allows users to scale existing CellChat workflows without changing downstream interpretation. The benchmark results support the use of CellChatAccelRcpp for large-scale comparative studies, repeated subsampling experiments and method-development workflows that require many CellChat runs.

Current limitations include support focused on RNA-based CellChat workflows and selected CellChat hot paths. Future work will extend acceleration coverage, broaden compatibility across CellChat versions and evaluate performance in spatial transcriptomic workflows.

## Figure legend

**Figure 1. CellChatAccelRcpp enables scalable and numerically matched CellChat inference.**  
(A) Schematic of the CellChatAccelRcpp acceleration layer and checkpointed benchmark workflow. (B) Runtime comparison between original CellChat and CellChatAccelRcpp across six target cell scales. Points represent individual dataset-repeat runs; lines and intervals show medians and interquartile ranges. (C) Runtime speedup distribution, calculated as original CellChat runtime divided by accelerated runtime. (D) Numerical equivalence measured by maximum absolute probability difference between original and accelerated communication probability tensors. The dashed red line marks 1e-12. (E) Component ablation showing runtime relative to the full accelerated path when the probability kernel, pathway aggregation or network aggregation is disabled. (F) Dataset-level heatmap of median speedup across cell scales.

## Key claims to support with final tables

- 864/864 benchmark tasks completed.
- 0 error metric rows.
- 216 paired original/accelerated comparisons.
- Median speedup: 11.4x overall.
- Median speedup by scale: 37.4x at 1k, 15.3x at 5k, 11.0x at 10k, 8.0x at 25k, 6.0x at 50k and 6.0x at all cells.
- Maximum absolute probability difference: 1.39e-16.
- Minimum Pearson correlation: 1.0.

## Files generated for this draft

- Main figure: `benchmarks/cellchat_acceleration_2026/results/figures/Fig01_runtime_compression.pdf`
- Main figure PNG: `benchmarks/cellchat_acceleration_2026/results/figures/Fig01_runtime_compression.png`
- Supplementary ablation figure: `benchmarks/cellchat_acceleration_2026/results/figures/Fig05_component_ablation.pdf`
- Supplementary accuracy figure: `benchmarks/cellchat_acceleration_2026/results/figures/Fig04_numerical_equivalence.pdf`
- Runtime summary: `benchmarks/cellchat_acceleration_2026/results/tables/publication_runtime_speedup_summary.tsv`
- Ablation summary: `benchmarks/cellchat_acceleration_2026/results/tables/publication_ablation_slowdown.tsv`
- Accuracy summary: `benchmarks/cellchat_acceleration_2026/results/tables/publication_accuracy_summary.tsv`

## TODO before submission

1. Add license, README, installation instructions and minimal test dataset.
2. Archive a release on Zenodo and add DOI.
3. Decide author order, affiliations, funding and acknowledgements.
4. Confirm whether the benchmark datasets can be redistributed or only referenced by accession/path.
5. Add references in journal style.
6. Include a short disclosure if AI-assisted writing or code generation is used.
