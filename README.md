# CellChatFastCpp

Fast Rcpp helpers for a common CellChat single-sample RNA workflow.

This package accelerates the CellChat steps that were slow in the tested Seurat RDS workflow:

- `computeCommunProbFastCpp()`
- `computeCommunProbPathwayFastCpp()`
- `aggregateNetFastCpp()`
- `computeAveExprFastCpp()`

The main speedup is in `computeCommunProbFastCpp()`.

## Benchmark On Two Local RDS Files

Tested files:

- `SCPCL000123_processed_seurat.rds`
- `SCPCL000125_processed_seurat.rds`

Original CellChat:

- `computeCommunProb`: `137.199 sec` total

Fast C++:

- `computeCommunProbFastCpp`: `11.055 sec` total
- Speedup: `12.41x`

Per sample:

| sample | CellChat computeCommunProb | Fast C++ |
| --- | ---: | ---: |
| SCPCL000123 | 50.705 sec | 3.600 sec |
| SCPCL000125 | 86.494 sec | 7.455 sec |

The fast C++ outputs were checked against original CellChat / previously validated C++ results:

- `prob` max absolute difference: floating point noise only
- `pval` max absolute difference: `0`
- pathway aggregation: equal within tolerance
- aggregate network count/weight: equal within tolerance

## Scope

This is not a full CellChat replacement. It is a focused accelerator for this setting:

- single-dataset CellChat object
- sc/snRNA-seq RNA data
- `type = "triMean"`
- `population.size = FALSE`
- non-spatial workflow
- CellChat v1-style object/API

It does not currently implement:

- spatial CellChat distance constraints
- `population.size = TRUE`
- `truncatedMean`, `thresholdedMean`, or `median`
- merged CellChat objects
- full CellChat v2 validation

If you change any of these assumptions, run the equivalence script before using the result at scale.

## Install

Install CellChat first. CellChat has moved over time; use the version that works in your environment.

```r
install.packages("remotes")
remotes::install_github("jinworks/CellChat")
```

Then install this package from the local clone:

```r
remotes::install_local("/path/to/CellChatFastCpp")
```

After you upload it to GitHub:

```r
remotes::install_github("YOUR_GITHUB_USERNAME/CellChatFastCpp")
```

## Single Object Usage

```r
library(Seurat)
library(CellChat)
library(CellChatFastCpp)

obj <- readRDS("sample_processed_seurat.rds")
DefaultAssay(obj) <- "RNA"

group_col <- "openscpca_celltype_annotation"
group <- as.character(obj@meta.data[[group_col]])
keep <- !is.na(group) & nzchar(group) & group != "openscpca-excluded"
obj <- subset(obj, cells = colnames(obj)[keep])

group <- as.character(obj@meta.data[[group_col]])
valid <- names(table(group))[table(group) >= 10]
obj <- subset(obj, cells = colnames(obj)[group %in% valid])

cellchat <- createCellChat(object = obj, group.by = group_col, assay = "RNA")
cellchat@DB <- CellChatDB.human

cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

cellchat <- computeCommunProbFastCpp(cellchat, nboot = 100, seed.use = 1L)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathwayFastCpp(cellchat)
cellchat <- aggregateNetFastCpp(cellchat)

saveRDS(cellchat, "sample_cellchat_fast.rds", compress = FALSE)
```

## Batch Usage

Run all `.rds` files in one directory:

```bash
Rscript scripts/run_cellchat_fast_batch.R \
  --input_dir /path/to/seurat_rds \
  --output_dir /path/to/cellchat_fast_results \
  --group_col openscpca_celltype_annotation \
  --pattern '\\.rds$' \
  --nboot 100 \
  --min_cells 10 \
  --species human
```

Recursive search:

```bash
Rscript scripts/run_cellchat_fast_batch.R \
  --input_dir /path/to/project \
  --output_dir /path/to/cellchat_fast_results \
  --group_col openscpca_celltype_annotation \
  --recursive true \
  --pattern '_processed_seurat\\.rds$'
```

Important options:

| option | default | meaning |
| --- | --- | --- |
| `--input_dir` | required | directory containing Seurat `.rds` files |
| `--output_dir` | `cellchat_fast_results` | result directory |
| `--group_col` | `openscpca_celltype_annotation` | Seurat metadata column for cell groups |
| `--pattern` | `\\.rds$` | regex used by `list.files()` |
| `--recursive` | `false` | recursively find input files |
| `--assay` | `RNA` | Seurat assay for CellChat |
| `--species` | `human` | `human` or `mouse` |
| `--min_cells` | `10` | remove groups with fewer cells |
| `--nboot` | `100` | bootstrap permutations |
| `--seed` | `1` | random seed |
| `--exclude_label` | `openscpca-excluded` | annotation label to remove |
| `--skip_existing` | `true` | skip files with existing output RDS |

## Batch Outputs

For each input sample:

- `{sample_id}_cellchat_fast.rds`
- `{sample_id}_LR_detail.csv`
- `{sample_id}_step_timings.csv`

Batch-level summaries:

- `batch_summary.csv`
- `all_step_timings.csv`
- `step_total_by_step.csv`

Use `step_total_by_step.csv` to see which step is still slow across a large run.

## Equivalence Check

Before a new dataset type or CellChat version, run:

```bash
Rscript scripts/check_equivalence_one.R \
  /path/to/one_sample_processed_seurat.rds \
  openscpca_celltype_annotation \
  5
```

This compares:

- original `CellChat::computeCommunProb()`
- `CellChatFastCpp::computeCommunProbFastCpp()`
- pathway aggregation
- network aggregation

Use a small `nboot` first because original CellChat may be slow.

## What Was Not Worth Rewriting

In the tested data:

- `identifyOverExpressedGenes` used `presto::wilcoxauc` and was already fast.
- `identifyOverExpressedInteractions` took about `0.257 sec` total.
- `filterCommunication` was about `0.001 sec` total.
- `subsetCommunication` was about `0.039 sec` total.
- `readRDS`, Seurat filtering, and `saveRDS` are mostly I/O/object-copying costs.

The meaningful target was `computeCommunProb`, especially repeated bootstrap aggregation and ligand/receptor expression calculation.

## GitHub Upload Notes

Do upload:

- `DESCRIPTION`
- `NAMESPACE`
- `R/`
- `src/`
- `scripts/`
- `README.md`
- `.gitignore`

Do not upload:

- local `.rds` files
- `Rlib/`
- large benchmark `results/`
- `.DS_Store`

The `.gitignore` already excludes these files.
