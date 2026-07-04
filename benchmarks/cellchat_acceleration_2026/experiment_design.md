# Experiment design for CellChat acceleration manuscript

## Manuscript target

Target article type: Bioinformatics software/original paper.

Core claim to support:

> The accelerated implementation reproduces CellChat communication probability
> results while reducing runtime and enabling larger real-world single-cell
> cohorts with checkpointable execution.

## Dataset panel

Use only real datasets as the primary evidence.

1. `3CA_data`: transferred complete under `/home/dzf/share/3CA_data`.
2. `normal_control`: transferred complete under `/home/dzf/share/normal_control`.
3. `Pediatric`: currently partial under `/home/dzf/share/Pediatric`; use after transfer resumes.

The benchmark should prefer ready-to-use `.rds` Seurat objects. `.h5ad`, `.h5`,
and `.mtx` files are included in the manifest as conversion candidates.

## E1. Numerical agreement

Purpose: prove the accelerated engine preserves the CellChat result.

Design:

- Use 6-12 representative Seurat RDS datasets.
- Downsample to 1k, 5k, 10k, 25k, 50k, and all available cells.
- Run baseline and accelerated engines with identical random seeds and grouping.
- Metrics:
  - maximum absolute difference in `cellchat@net$prob`
  - Pearson correlation of flattened probability tensors
  - number of ligand-receptor interactions
  - pathway-level interaction count agreement

Expected figure/table:

- Figure 1: probability scatter or Bland-Altman plot.
- Table 1: agreement statistics by dataset and scale.

## E2. Runtime and scalability

Purpose: show the algorithm supports large-scale workloads.

Design:

- Same dataset panel as E1.
- Three repeats per scale.
- Record elapsed time, success/failure, and optional peak memory.
- Report median and IQR.

Metrics:

- speedup = baseline median elapsed / accelerated median elapsed
- maximum completed cell count per engine
- throughput = cells processed per second

Expected figure/table:

- Figure 2: runtime vs number of cells, log scale.
- Figure 3: speedup vs number of cells.
- Table 2: largest completed workload per engine.

## E3. Cohort throughput

Purpose: show practical batch throughput.

Design:

- Run all eligible RDS datasets in `3CA_data` and `normal_control`.
- Use checkpoint/resume mode.
- Summarize total wall-clock time, completed samples, failed samples, and errors.

Expected figure/table:

- Figure 4: per-sample runtime distribution.
- Supplementary Table: per-sample status and ligand-receptor summary.

## E4. Checkpoint/resume reliability

Purpose: support "breakpoint experiment" claims.

Design:

- Select at least one medium and one large sample.
- Start accelerated run with checkpointing enabled.
- Interrupt after a fixed step or after N minutes.
- Restart with `--resume`.
- Compare resumed output to uninterrupted accelerated output.

Metrics:

- resumed result status
- skipped/completed checkpoint count
- max absolute probability difference vs uninterrupted run
- total time lost after interruption

Expected figure/table:

- Figure 5: resume timeline.
- Table 3: resumed vs uninterrupted agreement.

## E5. Ablation experiments

Purpose: isolate which design choices contribute speed.

Ablation modes used in the grid:

- `full`: all acceleration enabled.
- `no_accel_kernel`: disable accelerated probability kernel.
- `no_sparse_prefilter`: disable sparse/pre-filter optimization if supported.
- `no_parallel`: single-thread or serial execution if supported.

The benchmark script sets `CELLCHAT_ACCEL_ABLATION` and
`options(CellChatAccelRcpp.ablation = ...)`. The accelerated package should
read one of these controls or expose a matching function.

Expected figure/table:

- Figure 6: ablation runtime bars.
- Table 4: contribution of each component to speedup.

## Statistical reporting

- Use three repeats for runtime experiments.
- Report median and IQR.
- Use paired comparisons within the same dataset and scale.
- Do not use p-values as the main evidence; effect sizes and successful scale
  completion are more important for a software paper.

## Reproducibility checklist

- Pin package versions in `results/environment/R_packages.txt`.
- Keep the exact experiment grid in `results/experiment_grid.csv`.
- Keep one metrics file per experiment in `results/runs/`.
- Produce final tables from `code/summarize_results.R`.
- Do not hand-edit benchmark tables.

