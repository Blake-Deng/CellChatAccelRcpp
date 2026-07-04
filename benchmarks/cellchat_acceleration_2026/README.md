# Large-Scale CellChatAccelRcpp Benchmark

This directory contains the publication-oriented benchmark used to evaluate CellChatAccelRcpp against the original CellChat workflow.

The benchmark was designed for a Bioinformatics Application Note: paired original and accelerated runs were executed from the same prepared CellChat objects, followed by numerical-equivalence checks and component ablations.

This is the main paper benchmark. The separate [`../SCPCP000004`](../SCPCP000004) folder is an earlier cohort-specific validation benchmark and is kept for transparency, but the abstract, main figure and manuscript tables are based on this `cellchat_acceleration_2026` benchmark.

## Design

| component | setting |
| --- | --- |
| datasets | 12 real single-cell datasets from `normal_control` and `3CA_data` |
| cell scales | 1k, 5k, 10k, 25k, 50k, all available cells |
| repeats | 3 random seeds per dataset-scale pair |
| paired comparisons | original CellChat vs CellChatAccelRcpp |
| component ablations | probability kernel, pathway aggregation, network aggregation |
| completed jobs | 864 / 864 |
| failed metric files | 0 |

Raw input objects are not committed to GitHub because they are large and dataset-specific. The tracked files contain the manifest, experiment grid, environment record, summary tables, publication figures and scripts needed to reproduce the analysis when the source data are available.

## Main Results

| metric | result |
| --- | ---: |
| paired original/accelerated comparisons | 216 |
| overall median speedup | 11.4x |
| median original CellChat runtime | 426.6 s |
| median CellChatAccelRcpp runtime | 36.0 s |
| maximum absolute probability difference | 1.39e-16 |
| minimum probability Pearson correlation | 1.000 |

| target cells | median speedup |
| ---: | ---: |
| 1k | 37.4x |
| 5k | 15.3x |
| 10k | 11.0x |
| 25k | 8.0x |
| 50k | 6.0x |
| all | 6.0x |

## Publication Figures

The final figure files are in [`results/figures`](results/figures). PDF files are intended for manuscript submission; PNG files are included for GitHub preview.

| file | content |
| --- | --- |
| `Fig01_runtime_compression` | paired runtime compression by cell scale |
| `Fig02_speedup_distribution` | distribution of speedups across runs |
| `Fig03_paired_runtime_scatter` | paired original vs accelerated runtime |
| `Fig04_numerical_equivalence` | probability tensor numerical agreement |
| `Fig05_component_ablation` | contribution of each accelerated component |
| `Fig06_dataset_speedup_heatmap` | dataset-level speedup landscape |

Alternate Nature-style PNG versions generated during figure refinement are stored in [`results/figures_nature`](results/figures_nature).

## Tables

Summary tables are in [`results/tables`](results/tables).

| table | purpose |
| --- | --- |
| `publication_runtime_speedup_summary_v2.tsv` | final cell-scale runtime and speedup summary used for the manuscript |
| `publication_accuracy_summary_v2.tsv` | final numerical-equivalence summary used for the manuscript |
| `publication_ablation_summary_v2.tsv` | final component-ablation summary used for the manuscript |
| `publication_runtime_speedup_summary.tsv` | earlier publication runtime summary retained for provenance |
| `publication_accuracy_summary.tsv` | earlier publication accuracy summary retained for provenance |
| `publication_ablation_slowdown.tsv` | per-run component ablation results |
| `runtime_summary.tsv` | runtime summaries by dataset and scale |
| `accuracy_summary.tsv` | paired probability agreement details |
| `all_metrics.tsv` | complete benchmark metric table |

The manuscript uses the `_v2` publication summary tables where available. The larger detail tables remain in the repository so the summarized values can be audited.

## Reproducing The Benchmark

Create the environment:

```bash
mamba env create -f environment.yml
mamba activate cellchat-acceleration
```

Install the package from the repository root:

```bash
R CMD INSTALL .
```

Build the data manifest and experiment grid:

```bash
python code/01_build_manifest.py
python code/02_make_experiment_grid.py
```

Run the grid and summarize:

```bash
bash code/launch_full_benchmark.sh
Rscript code/summarize_results.R
ROOT="$(pwd)" Rscript code/04_make_nature_style_figures.R
```

The scripts default to the original server paths used during development. Set `ROOT`, `--manifest`, `--out` and input paths when reproducing the analysis on another machine.
