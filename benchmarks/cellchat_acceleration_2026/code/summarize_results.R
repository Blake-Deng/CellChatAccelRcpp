#!/usr/bin/env Rscript

root <- Sys.getenv("CELLCHAT_BENCH_ROOT", "/home/dzf/cellchat_acceleration")
run_dir <- file.path(root, "results", "runs")
table_dir <- file.path(root, "results", "tables")
figure_dir <- file.path(root, "results", "figures")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

files <- list.files(run_dir, pattern = "\\.metrics\\.tsv$", full.names = TRUE)
if (length(files) == 0) stop("No metrics files found in ", run_dir)

read_one <- function(path) {
  x <- utils::read.delim(path, stringsAsFactors = FALSE)
  x$metrics_file <- basename(path)
  x
}
metrics <- do.call(rbind, lapply(files, read_one))
utils::write.table(metrics, file.path(table_dir, "all_metrics.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

ok <- subset(metrics, status == "ok" & engine %in% c("baseline", "accelerated"))
if (nrow(ok) > 0) {
  runtime <- aggregate(
    elapsed_sec ~ dataset_id + n_cells_target + engine + ablation,
    data = ok,
    FUN = function(x) c(median = median(x), iqr = IQR(x), n = length(x))
  )
  runtime <- do.call(data.frame, runtime)
  names(runtime) <- sub("elapsed_sec\\.", "", names(runtime))
  utils::write.table(runtime, file.path(table_dir, "runtime_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
}

cmp <- subset(metrics, engine == "comparison")
if (nrow(cmp) > 0) {
  utils::write.table(cmp, file.path(table_dir, "accuracy_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
}

if (requireNamespace("ggplot2", quietly = TRUE) && nrow(ok) > 0) {
  library(ggplot2)
  p <- ggplot(ok, aes(x = n_cells_target, y = elapsed_sec, color = engine)) +
    geom_point(alpha = 0.7) +
    stat_summary(fun = median, geom = "line", aes(group = engine)) +
    facet_wrap(~ dataset_id, scales = "free_y") +
    theme_bw(base_size = 10) +
    labs(x = "Cells sampled", y = "Elapsed time (s)", color = "Engine")
  ggsave(file.path(figure_dir, "runtime_by_scale.pdf"), p, width = 10, height = 7)
}

cat("summary_tables=", table_dir, "\n")
cat("figures=", figure_dir, "\n")

