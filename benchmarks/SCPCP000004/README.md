# SCPCP000004 Supporting Validation

This directory is retained as supporting validation only. The current
manuscript and GitHub front page start from the v0.1.3 64-bit `sparse_stream`
workflow in `benchmarks/cellchat_acceleration_2026`.

This supporting benchmark compares the original CellChat R implementation with
the CellChatAccelRcpp workflow on the SCPCP000004 Seurat RDS collection.

The raw `.rds` CellChat objects and expression matrices are intentionally not
included in this repository. This directory contains only runnable scripts and
CSV summaries.

## Dataset

- Project: `SCPCP000004`
- Input files: 42 processed Seurat RDS files
- Input path used on the test server:
  `/home/dt2024/share/project/scpca_cellchat/scpca/SCPCP000004/rds`
- Grouping column: `openscpca_celltype_annotation`
- Cell filtering:
  - keep `scpca_filter == "Keep"` when available
  - remove empty labels and `openscpca-excluded`
  - remove groups with fewer than 10 cells
- CellChat parameters:
  - human database
  - RNA assay
  - `type = "triMean"`
  - `population.size = FALSE`
  - `nboot = 100`
  - `seed.use = 1`

## Run Summary

| method | successful samples | failed samples | successful elapsed time |
| --- | ---: | ---: | ---: |
| Original CellChat R | 40 | 2 | 15,202.23 sec |
| CellChatAccelRcpp | 40 | 2 | 1,805.11 sec |

Overall speedup on successful samples:

```text
15,202.23 / 1,805.11 = 8.42x
```

The two failed samples were identical between the two methods:

| sample | reason |
| --- | --- |
| `SCPCL000124` | fewer than 2 annotation groups after filtering |
| `SCPCL001058` | no significant ligand-receptor pairs after `identifyOverExpressedInteractions()` |

## Bottleneck

The main acceleration target remains `computeCommunProb()`.

| step | Original CellChat R | AccelRcpp |
| --- | ---: | ---: |
| communication probability inference | 14,526.05 sec | 957.60 sec |

Step-level speedup:

```text
14,526.05 / 957.60 = 15.17x
```

## Output Validation

The successful sample outputs were compared at the saved CellChat object level.

| validation | result |
| --- | ---: |
| samples compared | 40 |
| identical cell groups | 40 / 40 |
| identical LR pair counts | 40 / 40 |
| LR pair overlap | 100% |
| median probability matrix correlation | 1 |
| median network weight correlation | 1 |

This means the accelerated workflow reproduced the original CellChat R outputs
for the tested successful samples while reducing total runtime from about 4.22
hours to about 30.1 minutes.

## Files

Scripts:

- `scripts/run_SCPCP000004_official_cellchat_R.R`: runs the original CellChat R
  workflow on the RDS directory using the same filtering and annotation rules.
- `scripts/compare_SCPCP000004_official_R_vs_accelrcpp.R`: compares original R
  summaries against AccelRcpp summaries.

Result summaries:

- `results/official_R_run_summary.csv`
- `results/accelrcpp_run_summary.csv`
- `results/status_counts.csv`
- `results/sample_elapsed_comparison.csv`
- `results/computeCommunProb_step_comparison.csv`
- `results/step_total_comparison_long.csv`
- `results/object_level_output_validation.csv`

Full per-step timing tables:

- `results/official_R_all_step_timings.csv`
- `results/accelrcpp_all_step_timings.csv`
- `results/official_R_step_total_by_step.csv`
- `results/accelrcpp_step_total_by_step.csv`

## Reproduce

Run the original R baseline:

```bash
SCPCP_INPUT_DIR=/path/to/SCPCP000004/rds \
SCPCP_RESULT_DIR=/path/to/results/SCPCP000004_official_cellchat_R \
SCPCP_WORKERS=12 \
SCPCP_NBOOT=100 \
Rscript benchmarks/SCPCP000004/scripts/run_SCPCP000004_official_cellchat_R.R
```

Run the CellChatAccelRcpp workflow with the package batch script. The current
package default is `algorithm = "sparse_stream"`:

```bash
Rscript scripts/run_cellchat_accel_batch.R \
  --input_dir /path/to/SCPCP000004/rds \
  --output_dir /path/to/results/SCPCP000004 \
  --group_col openscpca_celltype_annotation \
  --pattern '_processed_seurat\\.rds$' \
  --nboot 100 \
  --min_cells 10 \
  --species human
```

Then compare summaries:

```bash
SCPCP_ACCELRCPP_DIR=/path/to/results/SCPCP000004 \
SCPCP_OFFICIAL_DIR=/path/to/results/SCPCP000004_official_cellchat_R \
SCPCP_COMPARE_DIR=/path/to/results/SCPCP000004_official_R_vs_accelrcpp_compare \
Rscript benchmarks/SCPCP000004/scripts/compare_SCPCP000004_official_R_vs_accelrcpp.R
```
