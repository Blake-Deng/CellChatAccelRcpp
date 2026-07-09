#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args) >= 1) args[[1]] else "benchmarks/cellchat_acceleration_2026/results/tables/sparse_exact_microbenchmark.tsv"

need <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package not available: ", pkg, call. = FALSE)
  }
}
need("CellChatAccelRcpp")

make_case <- function(G = 300L, K = 32L, nLR = 600L, nboot = 40L,
                      active_prob = 0.20, seed = 1L) {
  set.seed(seed)
  avg <- matrix(0, nrow = G, ncol = K)
  active <- matrix(stats::runif(G * K) < active_prob, nrow = G, ncol = K)
  avg[active] <- stats::rgamma(sum(active), shape = 2, rate = 8)

  avgBoot <- array(0, dim = c(G, K, nboot))
  for (b in seq_len(nboot)) {
    boot_active <- matrix(stats::runif(G * K) < active_prob, nrow = G, ncol = K)
    vals <- matrix(0, nrow = G, ncol = K)
    vals[boot_active] <- stats::rgamma(sum(boot_active), shape = 2, rate = 8)
    avgBoot[, , b] <- vals
  }

  ligandIdx <- matrix(sample.int(G, nLR, replace = TRUE), ncol = 1)
  receptorIdx <- matrix(sample.int(G, nLR, replace = TRUE), ncol = 1)
  emptyIdx <- matrix(0L, nrow = nLR, ncol = 1)
  list(
    avg = avg,
    avgBoot = avgBoot,
    ligandIdx = ligandIdx,
    receptorIdx = receptorIdx,
    coAIdx = emptyIdx,
    coIIdx = emptyIdx,
    agonistIdx = emptyIdx,
    antagonistIdx = emptyIdx,
    hasAgonist = rep(FALSE, nLR),
    hasAntagonist = rep(FALSE, nLR),
    Kh = 0.5,
    n_power = 1
  )
}

run_kernel <- function(fun, x) {
  do.call(fun, x)
}

time_kernel <- function(fun, x, reps = 3L) {
  elapsed <- numeric(reps)
  result <- NULL
  for (i in seq_len(reps)) {
    gc(FALSE)
    elapsed[[i]] <- system.time({ result <- run_kernel(fun, x) })[["elapsed"]]
  }
  list(result = result, elapsed = elapsed, median_elapsed = stats::median(elapsed))
}

cases <- expand.grid(
  K = c(16L, 32L, 48L),
  active_prob = c(0.15, 0.25),
  stringsAsFactors = FALSE
)
cases$nLR <- 600L
cases$nboot <- 40L
cases$G <- 300L
cases$seed <- seq_len(nrow(cases)) + 20260706L

rows <- vector("list", nrow(cases))
for (i in seq_len(nrow(cases))) {
  cfg <- cases[i, ]
  message(sprintf("case %d/%d: K=%d active_prob=%.2f", i, nrow(cases), cfg$K, cfg$active_prob))
  x <- make_case(G = cfg$G, K = cfg$K, nLR = cfg$nLR, nboot = cfg$nboot,
                 active_prob = cfg$active_prob, seed = cfg$seed)

  dense <- time_kernel(CellChatAccelRcpp:::cellchat_prob_from_avg_cpp, x)
  sparse <- time_kernel(CellChatAccelRcpp:::cellchat_prob_from_avg_sparse_cpp, x)

  prob_dense <- as.numeric(dense$result$prob)
  prob_sparse <- as.numeric(sparse$result$prob)
  pval_dense <- as.numeric(dense$result$pval)
  pval_sparse <- as.numeric(sparse$result$pval)

  rows[[i]] <- data.frame(
    G = cfg$G,
    K = cfg$K,
    nLR = cfg$nLR,
    nboot = cfg$nboot,
    active_prob = cfg$active_prob,
    dense_median_sec = dense$median_elapsed,
    sparse_median_sec = sparse$median_elapsed,
    sparse_speedup_vs_dense = dense$median_elapsed / sparse$median_elapsed,
    active_pairs = sparse$result$active_pairs,
    skipped_pairs = sparse$result$skipped_pairs,
    total_pairs = sparse$result$total_pairs,
    active_fraction = sparse$result$active_fraction,
    max_abs_prob_diff = max(abs(prob_dense - prob_sparse), na.rm = TRUE),
    max_abs_pval_diff = max(abs(pval_dense - pval_sparse), na.rm = TRUE),
    pearson_prob = suppressWarnings(stats::cor(prob_dense, prob_sparse, use = "pairwise.complete.obs")),
    pearson_pval = suppressWarnings(stats::cor(pval_dense, pval_sparse, use = "pairwise.complete.obs")),
    dense_elapsed_all = paste(sprintf("%.6f", dense$elapsed), collapse = ";"),
    sparse_elapsed_all = paste(sprintf("%.6f", sparse$elapsed), collapse = ";"),
    stringsAsFactors = FALSE
  )
}

res <- do.call(rbind, rows)
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
utils::write.table(res, out, sep = "\t", quote = FALSE, row.names = FALSE)
print(res[, c("K", "active_prob", "dense_median_sec", "sparse_median_sec",
              "sparse_speedup_vs_dense", "active_fraction",
              "max_abs_prob_diff", "max_abs_pval_diff")])
cat("wrote=", out, "\n", sep = "")
