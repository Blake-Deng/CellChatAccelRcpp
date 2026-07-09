#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default = NA_character_) {
  hit <- grep(paste0("^", name, "="), args, value = TRUE)
  if (length(hit) == 0) return(default)
  sub(paste0("^", name, "="), "", hit[[1]])
}

all_args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", all_args, value = TRUE)
script_file <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[[1]]) else NA_character_
default_root <- if (!is.na(script_file)) normalizePath(file.path(dirname(script_file), ".."), mustWork = FALSE) else getwd()

root <- get_arg("--root", Sys.getenv("CELLCHAT_BENCH_ROOT", default_root))
old_run_dir <- get_arg("--old-run-dir", "/home/dzf/cellchat_acceleration/results/runs")
sparse_run_dir <- get_arg("--sparse-run-dir", file.path(root, "results", "runs"))
table_dir <- get_arg("--table-dir", file.path(root, "results", "tables"))
compare_outputs <- tolower(get_arg("--compare-outputs", "TRUE")) %in% c("1", "true", "yes", "y")

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

read_metrics <- function(files) {
  rows <- lapply(files, function(path) {
    x <- utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
    x$metrics_file <- basename(path)
    x
  })
  cols <- unique(unlist(lapply(rows, names)))
  filled <- lapply(rows, function(x) {
    missing <- setdiff(cols, names(x))
    for (nm in missing) x[[nm]] <- NA
    x[cols]
  })
  do.call(rbind, filled)
}

pick_one <- function(x) {
  if (length(x) == 0) return(NA)
  x[[1]]
}

safe_num <- function(x) suppressWarnings(as.numeric(x))

old_files <- list.files(old_run_dir, pattern = "\\.metrics\\.tsv$", full.names = TRUE)
sparse_files <- list.files(sparse_run_dir, pattern = "__algorithm-sparse_exact\\.metrics\\.tsv$", full.names = TRUE)
if (length(old_files) == 0) stop("No old metrics files found in ", old_run_dir)
if (length(sparse_files) == 0) stop("No sparse_exact metrics files found in ", sparse_run_dir)

old_metrics <- read_metrics(old_files)
sparse_metrics <- read_metrics(sparse_files)
sparse_metrics$base_experiment_id <- sub("__algorithm-sparse_exact$", "", sparse_metrics$experiment_id)
old_metrics$base_experiment_id <- old_metrics$experiment_id

sparse_ids <- unique(sparse_metrics$base_experiment_id)
old_metrics <- old_metrics[old_metrics$base_experiment_id %in% sparse_ids, , drop = FALSE]

comparison <- do.call(rbind, lapply(sparse_ids, function(id) {
  old <- old_metrics[old_metrics$base_experiment_id == id, , drop = FALSE]
  sparse <- sparse_metrics[sparse_metrics$base_experiment_id == id, , drop = FALSE]
  baseline <- old[old$engine == "baseline", , drop = FALSE]
  dense <- old[old$engine == "accelerated", , drop = FALSE]
  old_cmp <- old[old$engine == "comparison", , drop = FALSE]

  baseline_sec <- safe_num(pick_one(baseline$elapsed_sec))
  dense_sec <- safe_num(pick_one(dense$elapsed_sec))
  sparse_sec <- safe_num(pick_one(sparse$elapsed_sec))

  data.frame(
    experiment_id = id,
    dataset_id = pick_one(sparse$dataset_id),
    n_cells_target = pick_one(sparse$n_cells_target),
    seed = pick_one(sparse$seed),
    baseline_status = pick_one(baseline$status),
    dense_status = pick_one(dense$status),
    sparse_status = pick_one(sparse$status),
    baseline_elapsed_sec = baseline_sec,
    dense_rcpp_elapsed_sec = dense_sec,
    sparse_exact_elapsed_sec = sparse_sec,
    baseline_vs_dense_speedup = baseline_sec / dense_sec,
    baseline_vs_sparse_exact_speedup = baseline_sec / sparse_sec,
    dense_vs_sparse_exact_speedup = dense_sec / sparse_sec,
    n_groups = safe_num(pick_one(sparse$n_groups)),
    total_lr_pairs = safe_num(pick_one(sparse$total_lr_pairs)),
    active_pairs = safe_num(pick_one(sparse$active_pairs)),
    skipped_pairs = safe_num(pick_one(sparse$skipped_pairs)),
    total_pairs_sparse = safe_num(pick_one(sparse$total_pairs_sparse)),
    active_fraction = safe_num(pick_one(sparse$active_fraction)),
    n_genes_kernel = safe_num(pick_one(sparse$n_genes_kernel)),
    old_baseline_dense_max_abs_prob_diff = safe_num(pick_one(old_cmp$max_abs_prob_diff)),
    old_baseline_dense_pearson_prob = safe_num(pick_one(old_cmp$pearson_prob)),
    stringsAsFactors = FALSE
  )
}))

comparison <- comparison[order(comparison$dataset_id, comparison$n_cells_target, comparison$seed), , drop = FALSE]
runtime_path <- file.path(table_dir, "sparse_exact_runtime_comparison.tsv")
utils::write.table(comparison, runtime_path, sep = "\t", quote = FALSE, row.names = FALSE)

summary_row <- data.frame(
  n_experiments = nrow(comparison),
  n_sparse_ok = sum(comparison$sparse_status == "ok", na.rm = TRUE),
  median_baseline_elapsed_sec = median(comparison$baseline_elapsed_sec, na.rm = TRUE),
  median_dense_rcpp_elapsed_sec = median(comparison$dense_rcpp_elapsed_sec, na.rm = TRUE),
  median_sparse_exact_elapsed_sec = median(comparison$sparse_exact_elapsed_sec, na.rm = TRUE),
  median_baseline_vs_dense_speedup = median(comparison$baseline_vs_dense_speedup, na.rm = TRUE),
  median_baseline_vs_sparse_exact_speedup = median(comparison$baseline_vs_sparse_exact_speedup, na.rm = TRUE),
  median_dense_vs_sparse_exact_speedup = median(comparison$dense_vs_sparse_exact_speedup, na.rm = TRUE),
  min_dense_vs_sparse_exact_speedup = min(comparison$dense_vs_sparse_exact_speedup, na.rm = TRUE),
  max_dense_vs_sparse_exact_speedup = max(comparison$dense_vs_sparse_exact_speedup, na.rm = TRUE),
  median_active_fraction = median(comparison$active_fraction, na.rm = TRUE),
  max_old_baseline_dense_max_abs_prob_diff = max(comparison$old_baseline_dense_max_abs_prob_diff, na.rm = TRUE),
  min_old_baseline_dense_pearson_prob = min(comparison$old_baseline_dense_pearson_prob, na.rm = TRUE),
  stringsAsFactors = FALSE
)

equiv_path <- file.path(table_dir, "sparse_exact_output_equivalence.tsv")
if (compare_outputs) {
  old_checkpoint_dir <- file.path(old_run_dir, "checkpoints")
  sparse_checkpoint_dir <- file.path(sparse_run_dir, "checkpoints")

  compare_one <- function(id) {
    dense_path <- file.path(old_checkpoint_dir, id, "accelerated.computed.rds")
    sparse_path <- file.path(sparse_checkpoint_dir, paste0(id, "__algorithm-sparse_exact"), "accelerated.computed.rds")
    if (!file.exists(dense_path) || !file.exists(sparse_path)) {
      return(data.frame(
        experiment_id = id,
        status = "missing_output",
        prob_max_abs_diff = NA_real_,
        pval_max_abs_diff = NA_real_,
        prob_pearson = NA_real_,
        pval_pearson = NA_real_,
        stringsAsFactors = FALSE
      ))
    }

    dense_obj <- readRDS(dense_path)
    sparse_obj <- readRDS(sparse_path)
    dense_prob <- as.numeric(dense_obj@net$prob)
    sparse_prob <- as.numeric(sparse_obj@net$prob)
    dense_pval <- as.numeric(dense_obj@net$pval)
    sparse_pval <- as.numeric(sparse_obj@net$pval)

    prob_cor <- if (stats::sd(dense_prob) == 0 && stats::sd(sparse_prob) == 0) {
      if (identical(dense_prob, sparse_prob)) 1 else NA_real_
    } else {
      suppressWarnings(stats::cor(dense_prob, sparse_prob))
    }
    pval_cor <- if (stats::sd(dense_pval) == 0 && stats::sd(sparse_pval) == 0) {
      if (identical(dense_pval, sparse_pval)) 1 else NA_real_
    } else {
      suppressWarnings(stats::cor(dense_pval, sparse_pval))
    }

    data.frame(
      experiment_id = id,
      status = "ok",
      prob_max_abs_diff = max(abs(dense_prob - sparse_prob), na.rm = TRUE),
      pval_max_abs_diff = max(abs(dense_pval - sparse_pval), na.rm = TRUE),
      prob_pearson = prob_cor,
      pval_pearson = pval_cor,
      stringsAsFactors = FALSE
    )
  }

  equiv <- do.call(rbind, lapply(comparison$experiment_id, compare_one))
  utils::write.table(equiv, equiv_path, sep = "\t", quote = FALSE, row.names = FALSE)
  summary_row$sparse_dense_equivalence_n <- nrow(equiv)
  summary_row$sparse_dense_equivalence_ok <- sum(equiv$status == "ok", na.rm = TRUE)
  summary_row$sparse_dense_max_abs_prob_diff <- max(equiv$prob_max_abs_diff, na.rm = TRUE)
  summary_row$sparse_dense_max_abs_pval_diff <- max(equiv$pval_max_abs_diff, na.rm = TRUE)
  summary_row$sparse_dense_min_prob_pearson <- min(equiv$prob_pearson, na.rm = TRUE)
  summary_row$sparse_dense_min_pval_pearson <- min(equiv$pval_pearson, na.rm = TRUE)
} else {
  summary_row$sparse_dense_equivalence_n <- NA_integer_
  summary_row$sparse_dense_equivalence_ok <- NA_integer_
  summary_row$sparse_dense_max_abs_prob_diff <- NA_real_
  summary_row$sparse_dense_max_abs_pval_diff <- NA_real_
  summary_row$sparse_dense_min_prob_pearson <- NA_real_
  summary_row$sparse_dense_min_pval_pearson <- NA_real_
}

summary_path <- file.path(table_dir, "sparse_exact_summary.tsv")
utils::write.table(summary_row, summary_path, sep = "\t", quote = FALSE, row.names = FALSE)

cat("runtime_comparison=", runtime_path, "\n", sep = "")
cat("summary=", summary_path, "\n", sep = "")
if (compare_outputs) cat("output_equivalence=", equiv_path, "\n", sep = "")
