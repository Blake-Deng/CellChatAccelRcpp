# CellChatAccelRcpp

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21186108.svg)](https://doi.org/10.5281/zenodo.21186108)
[![License: GPL-3](https://img.shields.io/badge/license-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html)
[![R](https://img.shields.io/badge/R-%3E%3D%204.1.0-blue.svg)](DESCRIPTION)

CellChatAccelRcpp is an R/Rcpp acceleration layer for large-scale CellChat RNA workflows. It keeps the standard CellChat object interface and replaces selected computational bottlenecks with compiled routines for communication probability estimation, pathway aggregation, network aggregation and group-level expression summaries.

The package is intended for users who need to run many CellChat analyses, larger cell sets or high-resolution groupings while preserving outputs that remain directly comparable with the original CellChat workflow.

Current release: `v0.1.1`

## What Is Accelerated

CellChatAccelRcpp provides accelerated replacements for the main single-dataset RNA workflow steps:

- `computeAveExprAccelRcpp()`: group-level `triMean` expression summaries
- `computeCommunProbAccelRcpp()`: ligand-receptor communication probability inference
- `computeCommunProbPathwayAccelRcpp()`: pathway-level aggregation
- `aggregateNetAccelRcpp()`: network aggregation

The `computeCommunProbAccelRcpp()` implementation includes dense and sparse kernels. The current recommended large-dataset path is:

```r
computeCommunProbAccelRcpp(
  cellchat,
  algorithm = "sparse_stream",
  nboot = 100,
  seed.use = 1L
)
```

The accelerated code preserves the CellChat probability model. It reduces interpreter overhead, avoids redundant work and streams sparse computations, but it does not remove the intrinsic sender-by-receiver group-pair term in the CellChat probability tensor.

## Benchmark Summary

In a paired benchmark across 12 real single-cell datasets, six target cell scales and three repeats, CellChatAccelRcpp completed 864 benchmark jobs without failed metric files.

| metric | result |
| --- | ---: |
| paired original/accelerated comparisons | 216 |
| overall median speedup | 11.4x |
| median original CellChat runtime | 426.6 s |
| median CellChatAccelRcpp runtime | 36.0 s |
| maximum absolute probability difference | 1.39e-16 |
| minimum probability Pearson correlation | 1.000 |

Median speedup by target cell scale:

| cells | median speedup |
| ---: | ---: |
| 1k | 37.4x |
| 5k | 15.3x |
| 10k | 11.0x |
| 25k | 8.0x |
| 50k | 6.0x |
| all available cells | 6.0x |

![Runtime compression](benchmarks/cellchat_acceleration_2026/results/figures/Fig01_runtime_compression.png)

Full benchmark scripts, publication figures and summary tables are in [`benchmarks/cellchat_acceleration_2026`](benchmarks/cellchat_acceleration_2026). Application Note material and LaTeX sources are in [`paper`](paper).

## Installation

Install CellChat and then install this package from GitHub:

```r
install.packages("remotes")
remotes::install_github("jinworks/CellChat")
remotes::install_github("Blake-Deng/CellChatAccelRcpp")
```

For reproducible benchmark work, use the conda environment file:

```bash
mamba env create -f benchmarks/cellchat_acceleration_2026/environment.yml
mamba activate cellchat-accelrcpp
```

## Single Object Usage

```r
library(CellChat)
library(CellChatAccelRcpp)

cellchat <- CellChat::subsetData(cellchat)
cellchat <- CellChat::identifyOverExpressedGenes(cellchat)
cellchat <- CellChat::identifyOverExpressedInteractions(cellchat)

cellchat <- computeCommunProbAccelRcpp(
  cellchat,
  algorithm = "sparse_stream",
  nboot = 100,
  seed.use = 1L
)
cellchat <- CellChat::filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathwayAccelRcpp(cellchat)
cellchat <- aggregateNetAccelRcpp(cellchat)
```

## Batch Usage

Run CellChatAccelRcpp on all Seurat `.rds` files in a directory:

```bash
Rscript scripts/run_cellchat_accel_batch.R \
  --input_dir /path/to/seurat_rds \
  --output_dir /path/to/cellchat_accel_results \
  --group_col openscpca_celltype_annotation \
  --pattern '\\.rds$' \
  --nboot 100 \
  --min_cells 10 \
  --species human
```

Before applying the package to a new dataset type or CellChat version, run an equivalence check:

```bash
Rscript scripts/check_equivalence_one.R \
  /path/to/one_sample_processed_seurat.rds \
  openscpca_celltype_annotation \
  5
```

## Supported Scope

CellChatAccelRcpp currently targets:

- single-dataset CellChat objects
- scRNA-seq and snRNA-seq RNA workflows
- `type = "triMean"`
- `population.size = FALSE`
- non-spatial CellChat workflows
- CellChat v1-style object/API

The package has not yet been validated for spatial distance constraints, `population.size = TRUE`, alternative mean functions, merged CellChat objects or full CellChat v2 workflows. Run the equivalence script before using new settings at scale.

## Repository Layout

```text
R/                      R interface for accelerated CellChat steps
src/                    Rcpp implementations and registration
scripts/                batch and equivalence-check scripts
benchmarks/             benchmark design, scripts, summaries and figures
paper/                  Bioinformatics Application Note draft material
NEWS.md                 release notes
```

## Release Notes

See [`NEWS.md`](NEWS.md) for versioned changes. The `v0.1.1` release adds the sparse streamed probability kernel, expanded benchmark outputs and updated publication material.

## Citation

Please cite the archived software release:

```text
Deng Z. CellChatAccelRcpp: scalable Rcpp acceleration of CellChat inference for large single-cell communication analyses.
DOI: https://doi.org/10.5281/zenodo.21186108
GitHub: https://github.com/Blake-Deng/CellChatAccelRcpp
```
