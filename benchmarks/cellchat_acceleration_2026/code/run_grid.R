#!/usr/bin/env Rscript

parse_args <- function(argv) {
  args <- list()
  i <- 1
  while (i <= length(argv)) {
    key <- argv[[i]]
    if (!startsWith(key, "--")) stop("Unexpected argument: ", key)
    key <- sub("^--", "", key)
    args[[key]] <- argv[[i + 1]]
    i <- i + 2
  }
  args
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
grid_path <- args[["grid"]]
task_id <- as.integer(args[["task-id"]])
root <- args[["root"]]
all_args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", all_args, value = TRUE)
script_file <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[[1]]) else NA_character_
default_root <- if (!is.na(script_file)) normalizePath(file.path(dirname(script_file), ".."), mustWork = FALSE) else getwd()
if (is.null(root)) root <- default_root
if (is.na(task_id) || task_id < 1) stop("--task-id must be a positive integer")

grid <- utils::read.csv(grid_path, stringsAsFactors = FALSE)
if (task_id > nrow(grid)) stop("task_id exceeds grid size")
row <- grid[task_id, ]

rscript <- file.path(R.home("bin"), "Rscript")
bench <- file.path(root, "code", "cellchat_benchmark.R")
out_dir <- file.path(root, "results", "runs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cmd_args <- c(
  bench,
  "--out-dir", out_dir,
  "--experiment-id", row$experiment_id,
  "--dataset-id", row$dataset_id,
  "--engine", row$engine,
  "--ablation", row$ablation,
  "--n-cells", as.character(row$n_cells),
  "--seed", as.character(row$seed),
  "--label-col", row$label_col,
  "--resume", "true"
)
if ("input_path" %in% names(row) && nzchar(row$input_path)) {
  cmd_args <- c(cmd_args, "--input", row$input_path)
}
if ("accel_algorithm" %in% names(row) && nzchar(row$accel_algorithm)) {
  cmd_args <- c(cmd_args, "--accel-algorithm", row$accel_algorithm)
}
if ("prepared_cellchat" %in% names(row) && nzchar(row$prepared_cellchat)) {
  cmd_args <- c(cmd_args, "--prepared-cellchat", row$prepared_cellchat)
}

status <- system2(rscript, cmd_args)
quit(save = "no", status = status)

