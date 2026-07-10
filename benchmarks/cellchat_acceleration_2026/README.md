# Large-Scale CellChatAccelRcpp Benchmark

This directory contains the publication-oriented benchmark used to evaluate CellChatAccelRcpp against the original CellChat workflow.

The benchmark was designed for a Bioinformatics Application Note: paired original and accelerated runs were executed from the same prepared CellChat objects, followed by numerical-equivalence checks and component ablations.

## Design

| component | setting |
| --- | --- |
| datasets | 12 real single-cell datasets from `normal_control` and `3CA_data` |
| cell scales | 1k, 5k, 10k, 25k, 50k, all available cells |
| repeats | 3 random seeds per dataset-scale pair |
| paired comparisons | original CellChat vs CellChatAccelRcpp |
| accelerated probability kernel | 64-bit `sparse_stream` |
| component ablations | probability kernel, pathway aggregation, network aggregation |
| completed jobs | 864 / 864 |
| failed metric files | 0 |

Raw input objects are not committed to GitHub because they are large and dataset-specific. The tracked files contain the manifest, experiment grid, environment record, summary tables, source plots and scripts needed to reproduce the analysis when the source data are available.

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

## Source Plots

The benchmark source plots are in [`results/figures`](results/figures). These files document the component plots used to assemble the current manuscript Figure 1 stored under [`../../paper/figures`](../../paper/figures).

The source plots cover runtime compression, speedup distributions, paired runtime scatter, numerical equivalence, component ablation and dataset-level speedup summaries.

## Tables

Summary tables are in [`benchmarks/cellchat_acceleration_2026/results/tables/`](results/tables/).

| table | purpose |
| --- | --- |
| `publication_runtime_speedup_summary.tsv` | cell-scale runtime and speedup summary |
| `publication_accuracy_summary.tsv` | numerical-equivalence summary |
| `publication_ablation_slowdown.tsv` | per-run component ablation results |
| `runtime_summary.tsv` | runtime summaries by dataset and scale |
| `accuracy_summary.tsv` | paired probability agreement details |
| `all_metrics.tsv` | complete benchmark metric table |

## Reproducing The Benchmark

Create the environment:

```bash
cd benchmarks/cellchat_acceleration_2026
mamba env create -f environment.yml
mamba activate cellchat-accelrcpp
```

Install the package from the benchmark directory:

```bash
R CMD INSTALL ../..
```

Build the data manifest and experiment grid:

```bash
python code/01_build_manifest.py --data-root /path/to/source_data
python code/02_make_experiment_grid.py
```

Run the grid and summarize:

```bash
DATA_ROOT=/path/to/source_data bash code/launch_full_benchmark.sh
Rscript code/summarize_results.R
ROOT="$(pwd)" Rscript code/04_make_nature_style_figures.R
```

The scripts write generated files under `benchmarks/cellchat_acceleration_2026/results/` by default. Set `DATA_ROOT`, `ROOT`, `--manifest` and `--out` when reproducing the analysis on another machine.
