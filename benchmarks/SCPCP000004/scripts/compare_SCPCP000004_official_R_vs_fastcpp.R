suppressPackageStartupMessages({
  user_lib <- "/home/dt2024/dt2024020307/R/x86_64-pc-linux-gnu-library/4.2"
  if (dir.exists(user_lib)) .libPaths(c(user_lib, .libPaths()))
})

options(stringsAsFactors = FALSE)

fast_dir <- Sys.getenv(
  "SCPCP_FASTCPP_DIR",
  "/home/dt2024/dt2024020307/cellchat/results/SCPCP000004"
)
official_dir <- Sys.getenv(
  "SCPCP_OFFICIAL_DIR",
  "/home/dt2024/dt2024020307/cellchat/results/SCPCP000004_official_cellchat_R"
)
out_dir <- Sys.getenv(
  "SCPCP_COMPARE_DIR",
  "/home/dt2024/dt2024020307/cellchat/results/SCPCP000004_official_R_vs_fastcpp_compare"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_csv_safe <- function(path) {
  if (!file.exists(path)) return(data.frame())
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

fast_summary <- read_csv_safe(file.path(fast_dir, "run_summary.csv"))
official_summary <- read_csv_safe(file.path(official_dir, "run_summary.csv"))
fast_timing <- read_csv_safe(file.path(fast_dir, "all_step_timings.csv"))
official_timing <- read_csv_safe(file.path(official_dir, "all_step_timings.csv"))

if (!nrow(fast_summary)) stop("Missing or empty FastCpp summary: ", file.path(fast_dir, "run_summary.csv"))
if (!nrow(official_summary)) stop("Missing or empty official summary: ", file.path(official_dir, "run_summary.csv"))

fast_summary$method <- "FastCpp"
official_summary$method <- "Official_R"
summary_all <- rbind(
  fast_summary[, intersect(names(fast_summary), names(official_summary)), drop = FALSE],
  official_summary[, intersect(names(fast_summary), names(official_summary)), drop = FALSE]
)
utils::write.csv(summary_all, file.path(out_dir, "combined_run_summary.csv"), row.names = FALSE)

status_tab <- aggregate(sample_id ~ method + status, summary_all, length)
names(status_tab)[names(status_tab) == "sample_id"] <- "n_samples"
utils::write.csv(status_tab, file.path(out_dir, "status_counts.csv"), row.names = FALSE)

fast_ok <- fast_summary[fast_summary$status == "success", , drop = FALSE]
official_ok <- official_summary[official_summary$status == "success", , drop = FALSE]
sample_cmp <- merge(
  fast_ok[, c("sample_id", "elapsed_sec", "n_cells_used", "n_groups", "n_lr", "n_pathways")],
  official_ok[, c("sample_id", "elapsed_sec", "n_cells_used", "n_groups", "n_lr", "n_pathways")],
  by = "sample_id",
  suffixes = c("_fastcpp", "_official_R")
)
if (nrow(sample_cmp)) {
  sample_cmp$speedup_total <- sample_cmp$elapsed_sec_official_R / sample_cmp$elapsed_sec_fastcpp
  sample_cmp <- sample_cmp[order(sample_cmp$speedup_total, decreasing = TRUE), , drop = FALSE]
}
utils::write.csv(sample_cmp, file.path(out_dir, "sample_elapsed_comparison.csv"), row.names = FALSE)

step_summary <- function(df, method) {
  if (!nrow(df)) return(data.frame())
  out <- aggregate(elapsed_sec ~ step, df, sum)
  names(out)[names(out) == "elapsed_sec"] <- "total_elapsed_sec"
  out$method <- method
  out[, c("method", "step", "total_elapsed_sec")]
}

step_all <- rbind(step_summary(fast_timing, "FastCpp"), step_summary(official_timing, "Official_R"))
utils::write.csv(step_all, file.path(out_dir, "step_total_comparison_long.csv"), row.names = FALSE)

fast_prob_steps <- fast_timing[grepl("computeCommunProb", fast_timing$step), , drop = FALSE]
official_prob_steps <- official_timing[grepl("computeCommunProb", official_timing$step), , drop = FALSE]
prob_cmp <- merge(
  fast_prob_steps[, c("sample_id", "step", "elapsed_sec")],
  official_prob_steps[, c("sample_id", "step", "elapsed_sec")],
  by = "sample_id",
  suffixes = c("_fastcpp", "_official_R")
)
if (nrow(prob_cmp)) {
  prob_cmp$speedup_computeCommunProb <- prob_cmp$elapsed_sec_official_R / prob_cmp$elapsed_sec_fastcpp
  prob_cmp <- prob_cmp[order(prob_cmp$speedup_computeCommunProb, decreasing = TRUE), , drop = FALSE]
}
utils::write.csv(prob_cmp, file.path(out_dir, "computeCommunProb_step_comparison.csv"), row.names = FALSE)

cat("Comparison written to:", out_dir, "\n")
cat("FastCpp success:", sum(fast_summary$status == "success"), "/", nrow(fast_summary), "\n")
cat("Official R success:", sum(official_summary$status == "success"), "/", nrow(official_summary), "\n")
if (nrow(sample_cmp)) {
  cat("Median total speedup official/FastCpp:", stats::median(sample_cmp$speedup_total, na.rm = TRUE), "\n")
}
if (nrow(prob_cmp)) {
  cat("Median computeCommunProb speedup official/FastCpp:", stats::median(prob_cmp$speedup_computeCommunProb, na.rm = TRUE), "\n")
}
