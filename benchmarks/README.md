# Benchmarks

This directory contains the primary v0.1.3 `sparse_stream` benchmark track and one supporting validation track. The GitHub and manuscript narrative starts from the `sparse_stream` benchmark.

## Benchmark Tracks

| directory | role | status | main result |
| --- | --- | --- | --- |
| [`cellchat_acceleration_2026`](cellchat_acceleration_2026) | Main v0.1.3 `sparse_stream` benchmark for the Bioinformatics Application Note | Complete | 12 real datasets, 216 paired original/accelerated comparisons, 648 ablation runs, median 11.4x speedup |
| [`SCPCP000004`](SCPCP000004) | Supporting dataset-specific validation on the OpenScPCA SCPCP000004 cohort | Complete | 40 successful sample-level comparisons, 8.42x end-to-end speedup, matched CellChat outputs |

## What To Cite In The Manuscript

The manuscript benchmark uses [`cellchat_acceleration_2026`](cellchat_acceleration_2026). This is the paper-facing `sparse_stream` benchmark with the final tables, source plot scripts, environment records and Supplementary Table S1 links.

[`SCPCP000004`](SCPCP000004) is retained as supporting evidence that the accelerated workflow also reproduces a separate cohort-style CellChat analysis. It is useful for repository transparency, but it is not the primary benchmark summarized in the Application Note abstract.

## Main Publication Benchmark Contents

Key files under [`cellchat_acceleration_2026`](cellchat_acceleration_2026):

- `README.md`: benchmark-specific summary and reproduction instructions.
- `experiment_design.md`: design rationale for agreement, runtime, checkpointing and ablation experiments.
- `environment.yml`: conda environment used for benchmark work.
- `data_manifest.tsv` and `dataset_candidates.tsv`: dataset inventory and selected RDS inputs.
- `code/`: manifest construction, grid generation, paired benchmark runner, summarization and figure scripts.
- `results/tables/`: processed benchmark tables used by the manuscript.
- `results/figures/`: PDF and PNG source plots used to assemble the current manuscript figure.
- `results/figures_nature/`: alternate PNG source plots produced during figure refinement.
- `results/environment/`: package and system environment records.

Raw `.rds`, `.h5ad`, `.h5`, `.mtx` and other large source data files are intentionally not stored in this GitHub repository.
