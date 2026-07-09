#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args) >= 1) args[[1]] else "benchmarks/cellchat_acceleration_2026/results/tables/sparse_exact_real_check.tsv"
nboot <- as.integer(Sys.getenv("SPARSE_EXACT_NBOOT", "20"))
seed <- as.integer(Sys.getenv("SPARSE_EXACT_SEED", "20260706"))

need <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Package not available: ", pkg, call. = FALSE)
}
need("CellChatAccelRcpp")

`%||%` <- function(a, b) if (is.null(a)) b else a

prepared_paths <- c(
  namecheck_3ca_10x_n200 = "/home/dzf/cellchat_acceleration/results/namecheck_20260704/checkpoints/namecheck_3ca_10x_accel_n200/prepared_cellchat.rds",
  prjna490728_5k = "/home/dzf/cellchat_acceleration/results/runs/checkpoints/3CA_data__PRJNA490728__cells-5000__rep-2__engine-both__ablation-full/prepared_cellchat.rds"
)
prepared_paths <- prepared_paths[file.exists(prepared_paths)]
if (!length(prepared_paths)) stop("No prepared CellChat objects found.", call. = FALSE)

compare_vec <- function(a, b) {
  x <- as.numeric(a)
  y <- as.numeric(b)
  list(
    max_abs = max(abs(x - y), na.rm = TRUE),
    cor = suppressWarnings(stats::cor(x, y, use = "pairwise.complete.obs"))
  )
}

time_call <- function(expr) {
  gc(FALSE)
  result <- NULL
  elapsed <- system.time({ result <- force(expr) })[["elapsed"]]
  list(result = result, elapsed = elapsed)
}

rows <- list()
for (nm in names(prepared_paths)) {
  message("reading ", nm)
  obj <- readRDS(prepared_paths[[nm]])
  dims <- dim(obj@data.signaling)
  cells <- if (length(dims) >= 2) dims[[2]] else length(obj@idents)
  genes <- if (length(dims) >= 1) dims[[1]] else NA_integer_
  groups <- length(levels(obj@idents))
  lr <- nrow(obj@LR$LRsig)

  message("full dense ", nm)
  dense <- time_call(CellChatAccelRcpp::computeCommunProbAccelRcpp(
    obj, nboot = nboot, seed.use = seed, algorithm = "dense"
  ))
  message("full sparse ", nm)
  sparse <- time_call(CellChatAccelRcpp::computeCommunProbAccelRcpp(
    obj, nboot = nboot, seed.use = seed, algorithm = "sparse_exact"
  ))
  message("full subset sparse ", nm)
  subset_sparse <- time_call(CellChatAccelRcpp::computeCommunProbAccelRcpp(
    obj, nboot = nboot, seed.use = seed, algorithm = "sparse_exact_subset_boot"
  ))
  message("full ondemand simple ", nm)
  ondemand_sparse <- time_call(CellChatAccelRcpp::computeCommunProbAccelRcpp(
    obj, nboot = nboot, seed.use = seed, algorithm = "sparse_exact_ondemand_simple"
  ))

  add_full_row <- function(result, elapsed, level) {
    prob_cmp <- compare_vec(dense$result@net$prob, result@net$prob)
    pval_cmp <- compare_vec(dense$result@net$pval, result@net$pval)
    stats <- result@options$accelrcpp
    data.frame(
      object = nm,
      level = level,
      cells = cells,
      genes = genes,
      groups = groups,
      lr = lr,
      nboot = nboot,
      dense_sec = dense$elapsed,
      sparse_sec = elapsed,
      sparse_speedup_vs_dense = dense$elapsed / elapsed,
      active_pairs = stats$active_pairs %||% NA_real_,
      skipped_pairs = stats$skipped_pairs %||% NA_real_,
      total_pairs = stats$total_pairs %||% NA_real_,
      active_fraction = stats$active_fraction %||% NA_real_,
      n_genes_kernel = stats$n_genes_kernel %||% NA_real_,
      n_simple_ondemand_lr = stats$n_simple_ondemand_lr %||% NA_real_,
      n_other_lr = stats$n_other_lr %||% NA_real_,
      boot_tri_mean_evals = stats$boot_tri_mean_evals %||% NA_real_,
      max_abs_prob_diff = prob_cmp$max_abs,
      max_abs_pval_diff = pval_cmp$max_abs,
      pearson_prob = prob_cmp$cor,
      pearson_pval = pval_cmp$cor,
      stringsAsFactors = FALSE
    )
  }

  rows[[length(rows) + 1L]] <- add_full_row(sparse$result, sparse$elapsed, "full_sparse_exact")
  rows[[length(rows) + 1L]] <- add_full_row(subset_sparse$result, subset_sparse$elapsed, "full_sparse_exact_subset_boot")
  rows[[length(rows) + 1L]] <- add_full_row(ondemand_sparse$result, ondemand_sparse$elapsed, "full_sparse_exact_ondemand_simple")


  message("kernel setup ", nm)
  kernel_row <- tryCatch({
    data <- as.matrix(obj@data.signaling)
    max_data <- max(data)
    data.use <- data / max_data
    group <- obj@idents
    K <- nlevels(group)
    group_int <- as.integer(group)
    data.use.avg <- CellChatAccelRcpp:::group_tri_mean_cpp(data.use, group_int, K)
    colnames(data.use.avg) <- levels(group)
    rownames(data.use.avg) <- rownames(data.use)

    set.seed(seed)
    permutation <- replicate(nboot, sample.int(ncol(data.use), size = ncol(data.use)))
    group_boot <- matrix(group_int[permutation], nrow = ncol(data.use), ncol = nboot)
    avg_boot <- CellChatAccelRcpp:::group_tri_mean_boot_cpp(data.use, group_boot, K)$avg_boot

    pairLRsig <- obj@LR$LRsig
    complex_input <- obj@DB$complex
    cofactor_input <- obj@DB$cofactor
    data_genes <- rownames(data.use)
    cofactor_cols <- grepl("cofactor", colnames(cofactor_input))
    ligand_idx <- CellChatAccelRcpp:::make_index_matrix(as.character(pairLRsig$ligand), data_genes,
                                                        complex_input = complex_input, keep_duplicates = TRUE)
    receptor_idx <- CellChatAccelRcpp:::make_index_matrix(as.character(pairLRsig$receptor), data_genes,
                                                          complex_input = complex_input, keep_duplicates = TRUE)
    coA_idx <- CellChatAccelRcpp:::make_index_matrix(as.character(pairLRsig$co_A_receptor), data_genes,
                                                     cofactor_input = cofactor_input, cofactor_cols = cofactor_cols)
    coI_idx <- CellChatAccelRcpp:::make_index_matrix(as.character(pairLRsig$co_I_receptor), data_genes,
                                                     cofactor_input = cofactor_input, cofactor_cols = cofactor_cols)
    agonist_names <- as.character(pairLRsig$agonist)
    antagonist_names <- as.character(pairLRsig$antagonist)
    hasAgonist <- !is.na(agonist_names) & nzchar(agonist_names)
    hasAntagonist <- !is.na(antagonist_names) & nzchar(antagonist_names)
    agonist_idx <- CellChatAccelRcpp:::make_index_matrix(agonist_names, data_genes,
                                                        cofactor_input = cofactor_input, cofactor_cols = cofactor_cols)
    antagonist_idx <- CellChatAccelRcpp:::make_index_matrix(antagonist_names, data_genes,
                                                           cofactor_input = cofactor_input, cofactor_cols = cofactor_cols)
    kernel_args <- list(
      avg = data.use.avg,
      avgBoot = avg_boot,
      ligandIdx = ligand_idx,
      receptorIdx = receptor_idx,
      coAIdx = coA_idx,
      coIIdx = coI_idx,
      agonistIdx = agonist_idx,
      antagonistIdx = antagonist_idx,
      hasAgonist = hasAgonist,
      hasAntagonist = hasAntagonist,
      Kh = 0.5,
      n_power = 1
    )

    message("kernel dense ", nm)
    dense_kernel <- time_call(do.call(CellChatAccelRcpp:::cellchat_prob_from_avg_cpp, kernel_args))
    message("kernel sparse ", nm)
    sparse_kernel <- time_call(do.call(CellChatAccelRcpp:::cellchat_prob_from_avg_sparse_cpp, kernel_args))
    kernel_prob_cmp <- compare_vec(dense_kernel$result$prob, sparse_kernel$result$prob)
    kernel_pval_cmp <- compare_vec(dense_kernel$result$pval, sparse_kernel$result$pval)

    data.frame(
      object = nm,
      level = "kernel_only_after_common_precompute",
      cells = cells,
      genes = genes,
      groups = groups,
      lr = lr,
      nboot = nboot,
      dense_sec = dense_kernel$elapsed,
      sparse_sec = sparse_kernel$elapsed,
      sparse_speedup_vs_dense = dense_kernel$elapsed / sparse_kernel$elapsed,
      active_pairs = sparse_kernel$result$active_pairs %||% NA_real_,
      skipped_pairs = sparse_kernel$result$skipped_pairs %||% NA_real_,
      total_pairs = sparse_kernel$result$total_pairs %||% NA_real_,
      active_fraction = sparse_kernel$result$active_fraction %||% NA_real_,
      n_genes_kernel = NA_real_,
      n_simple_ondemand_lr = NA_real_,
      n_other_lr = NA_real_,
      boot_tri_mean_evals = NA_real_,
      max_abs_prob_diff = kernel_prob_cmp$max_abs,
      max_abs_pval_diff = kernel_pval_cmp$max_abs,
      pearson_prob = kernel_prob_cmp$cor,
      pearson_pval = kernel_pval_cmp$cor,
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    warning("kernel-only check failed for ", nm, ": ", conditionMessage(e))
    NULL
  })
  if (!is.null(kernel_row)) rows[[length(rows) + 1L]] <- kernel_row
}

res <- do.call(rbind, rows)
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
utils::write.table(res, out, sep = "\t", quote = FALSE, row.names = FALSE)
print(res)
cat("wrote=", out, "\n", sep = "")
