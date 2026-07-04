# CellChatAccelRcpp enables scalable cell-cell communication inference for large single-cell atlases

## Abstract

### Motivation

Cell-cell communication inference is widely used to interpret single-cell
transcriptomic data, but the computational cost of CellChat can become a
practical bottleneck when analyses are expanded from individual samples to
large cohorts and atlases. This limits routine use of repeated benchmarking,
parameter exploration, and checkpointed large-scale workflows.

### Results

We present CellChatAccelRcpp, an accelerated implementation designed to preserve
the CellChat analysis interface while reducing the runtime of communication
probability computation. Across real single-cell datasets from 3CA,
normal-control, and pediatric atlas resources, CellChatAccelRcpp reproduced the
baseline CellChat communication probability outputs with [TO_FILL: agreement
statistics] while reducing median runtime by [TO_FILL: speedup] and enabling
successful completion of datasets up to [TO_FILL: maximum cell count]. A
checkpoint/resume workflow allowed interrupted large-scale runs to continue
without restarting completed steps, producing results consistent with
uninterrupted runs.

### Availability and implementation

The software is implemented in R and Rcpp. Source code, benchmark scripts,
experiment grids, and processed benchmark outputs will be made available at
[TO_FILL: public repository URL]. The exact datasets or access instructions
used for benchmarking will be provided in the Data Availability statement.

### Contact

[TO_FILL]

## 1 Introduction

Single-cell transcriptomics has made it possible to study intercellular
communication across diverse tissues, disease states, and developmental
contexts. CellChat is a widely used framework for inferring ligand-receptor
communication patterns from single-cell expression matrices. However, as
studies move from single datasets to multi-sample atlases, computational cost
can constrain practical workflows. Large cell counts, repeated downsampling,
batch analyses, and parameter sweeps all multiply this cost.

Here we introduce CellChatAccelRcpp, an acceleration layer for CellChat-focused
communication analysis. The method targets the computationally intensive
probability calculation stage while maintaining compatibility with standard
CellChat objects. We evaluate the method with real datasets and emphasize four
questions required for a robust software manuscript: numerical agreement,
runtime scalability, large-cohort throughput, and checkpoint/restart behavior.

## 2 Materials and methods

### 2.1 Datasets

Benchmarks use real single-cell datasets stored under `/home/dzf/share`,
including 3CA-derived datasets, normal-control datasets, and pediatric atlas
data. The primary benchmark uses Seurat RDS objects when available, because
these can be analyzed directly by CellChat. H5AD, H5, and MTX files are
recorded in the data manifest as conversion candidates.

Final manuscript text should replace this paragraph with accession identifiers,
sample counts, cell counts, tissue/disease descriptions, and preprocessing
steps.

### 2.2 Baseline CellChat workflow

For each dataset, the baseline workflow constructs a CellChat object from a
Seurat object, assigns cell groups from existing annotation columns when
available, and otherwise derives clusters using Seurat normalization, variable
feature selection, scaling, PCA, nearest-neighbor graph construction, and
clustering. The workflow then runs the standard CellChat sequence:
`subsetData`, `identifyOverExpressedGenes`,
`identifyOverExpressedInteractions`, `computeCommunProb`,
`filterCommunication`, `computeCommunProbPathway`, and `aggregateNet`.

### 2.3 Accelerated workflow

The accelerated workflow uses the same input objects, group labels, random
seeds, and downstream CellChat steps. The baseline `computeCommunProb` call is
replaced by the CellChatAccelRcpp accelerated probability computation. The
benchmark wrapper records the accelerated function name and package version in
the environment report.

### 2.4 Scalability benchmark

Each dataset is evaluated at increasing cell counts: 1,000, 5,000, 10,000,
25,000, 50,000, and all available cells. Each run is repeated three times.
Runtime is summarized as median and interquartile range. Speedup is computed as
the median baseline runtime divided by the median accelerated runtime for the
same dataset and scale.

### 2.5 Numerical agreement

Agreement is assessed by comparing flattened CellChat communication probability
tensors between baseline and accelerated runs. Metrics include maximum absolute
difference, Pearson correlation, number of inferred ligand-receptor
interactions, and pathway-level network summaries. Agreement is evaluated
within the same dataset, scale, and random seed.

### 2.6 Checkpoint/resume experiment

The checkpoint experiment interrupts a large accelerated run after a defined
time point or checkpoint boundary, then restarts the same experiment with
resume enabled. The resumed result is compared with an uninterrupted run using
the same agreement metrics. The manuscript should report whether any completed
step was recomputed, how much time was lost after interruption, and whether the
final result matched the uninterrupted run.

### 2.7 Ablation analysis

Ablation experiments disable individual acceleration components where supported:
accelerated kernel, sparse pre-filtering, and parallel execution. Each ablation
is run on the same dataset-scale pairs as the full accelerated method. The
effect of each component is reported as runtime change relative to the full
accelerated method and relative to baseline CellChat.

## 3 Results

### 3.1 Dataset inventory and benchmark coverage

The benchmark manifest identified [TO_FILL] candidate files across
`3CA_data`, `normal_control`, and `Pediatric`, including [TO_FILL] directly
usable RDS objects. The final benchmark grid included [TO_FILL] datasets,
[TO_FILL] cell-count scales, and [TO_FILL] repeated runs.

### 3.2 CellChatAccelRcpp preserves CellChat communication results

Across matched baseline and accelerated runs, communication probability tensors
showed [TO_FILL] agreement. The maximum absolute difference was [TO_FILL], and
the Pearson correlation was [TO_FILL]. Ligand-receptor interaction counts and
pathway-level summaries were [TO_FILL].

### 3.3 Acceleration improves runtime and scale

CellChatAccelRcpp reduced median runtime by [TO_FILL] relative to baseline
CellChat. The largest completed baseline workload contained [TO_FILL] cells,
whereas the accelerated workflow completed [TO_FILL] cells. Runtime gains were
most pronounced at [TO_FILL] scale.

### 3.4 Checkpointing supports robust large-scale experiments

After forced interruption, the checkpointed accelerated workflow resumed from
completed outputs and finished successfully in [TO_FILL]. The resumed output
matched the uninterrupted accelerated run with [TO_FILL] agreement.

### 3.5 Ablation identifies the main acceleration contributors

Disabling [TO_FILL] produced the largest runtime regression, indicating that
[TO_FILL] contributes most to the observed speedup. The full accelerated method
remained the most efficient configuration across [TO_FILL] datasets.

## 4 Discussion

CellChatAccelRcpp addresses a practical limitation in large-scale cell-cell
communication studies: the cost of repeated CellChat inference across many
samples and cell-count scales. The proposed benchmark emphasizes real datasets,
agreement with the established baseline, and operational behavior under
interruption. These properties are important for atlas-scale workflows where
analyses may run for many hours and need to be resumed reliably.

Limitations should be reported transparently. If acceleration is limited to
specific CellChat stages, the manuscript should avoid claiming full-pipeline
acceleration. If memory use remains high for specific datasets, report the
failure modes and explain when users should downsample or split analyses.

## 5 Conclusion

CellChatAccelRcpp provides a scalable path for CellChat-style communication
analysis in large single-cell cohorts while preserving baseline output
agreement. The checkpointable benchmark framework supports reproducible
evaluation and can be reused as the software evolves.

## Availability of data and materials

[TO_FILL: public repository, version tag, Zenodo/OSF archive if applicable,
dataset accession IDs or access instructions.]

## Funding

[TO_FILL]

## Conflict of interest

[TO_FILL]

## References

[TO_FILL]
