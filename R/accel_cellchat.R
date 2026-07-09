computeCommunProbAccelRcpp <- function(object, type = "triMean", raw.use = TRUE, population.size = FALSE,
                                     nboot = 100, seed.use = 1L, Kh = 0.5, n = 1,
                                     distance.use = NULL,
                                     algorithm = c("sparse_stream", "dense", "sparse_exact", "sparse_exact_subset_boot", "sparse_exact_ondemand_simple")) {
  algorithm <- match.arg(algorithm)
  check_supported_object(object, type, population.size, distance.use)

  ptm <- Sys.time()
  data <- if (raw.use) {
    as.matrix(object@data.signaling)
  } else {
    if (!"data.smooth" %in% methods::slotNames(object)) {
      stop("object@data.smooth is missing. Run CellChat::projectData() first or use raw.use = TRUE.", call. = FALSE)
    }
    as.matrix(object@data.smooth)
  }

  if (!length(data)) stop("No signaling data found in object@data.signaling.", call. = FALSE)
  max_data <- max(data)
  if (!is.finite(max_data) || max_data <= 0) {
    stop("Signaling expression matrix has no positive finite values.", call. = FALSE)
  }

  pairLRsig <- object@LR$LRsig
  if (is.null(pairLRsig) || !nrow(pairLRsig)) {
    stop("object@LR$LRsig is empty. Run identifyOverExpressedInteractions() first.", call. = FALSE)
  }

  complex_input <- object@DB$complex
  cofactor_input <- object@DB$cofactor
  if (is.null(complex_input) || is.null(cofactor_input)) {
    stop("object@DB must include CellChat complex and cofactor tables.", call. = FALSE)
  }

  group <- object@idents
  if (!is.factor(group)) {
    stop("object@idents must be a factor for single-dataset CellChat objects.", call. = FALSE)
  }
  K <- nlevels(group)
  if (K != length(unique(group))) {
    stop("Unused factor levels found in object@idents. Drop unused levels before running.", call. = FALSE)
  }

  data.use <- data / max_data
  nC <- ncol(data.use)
  group_int <- as.integer(group)

  data_genes <- rownames(data.use)
  cofactor_cols <- grepl("cofactor", colnames(cofactor_input))
  ligand_idx <- make_index_matrix(as.character(pairLRsig$ligand), data_genes,
                                  complex_input = complex_input, keep_duplicates = TRUE)
  receptor_idx <- make_index_matrix(as.character(pairLRsig$receptor), data_genes,
                                    complex_input = complex_input, keep_duplicates = TRUE)
  coA_idx <- make_index_matrix(as.character(pairLRsig$co_A_receptor), data_genes,
                               cofactor_input = cofactor_input, cofactor_cols = cofactor_cols)
  coI_idx <- make_index_matrix(as.character(pairLRsig$co_I_receptor), data_genes,
                               cofactor_input = cofactor_input, cofactor_cols = cofactor_cols)
  agonist_names <- as.character(pairLRsig$agonist)
  antagonist_names <- as.character(pairLRsig$antagonist)
  hasAgonist <- !is.na(agonist_names) & nzchar(agonist_names)
  hasAntagonist <- !is.na(antagonist_names) & nzchar(antagonist_names)
  agonist_idx <- make_index_matrix(agonist_names, data_genes,
                                   cofactor_input = cofactor_input, cofactor_cols = cofactor_cols)
  antagonist_idx <- make_index_matrix(antagonist_names, data_genes,
                                      cofactor_input = cofactor_input, cofactor_cols = cofactor_cols)

  subset_boot <- algorithm %in% c("sparse_exact_subset_boot", "sparse_exact_ondemand_simple")
  ondemand_simple <- identical(algorithm, "sparse_exact_ondemand_simple")
  sparse_stream <- identical(algorithm, "sparse_stream")
  n_genes_full <- nrow(data.use)
  n_genes_kernel <- n_genes_full

  set.seed(seed.use)
  permutation <- replicate(nboot, sample.int(nC, size = nC))
  group_boot <- matrix(group_int[permutation], nrow = nC, ncol = nboot)

  row_count <- function(idx) rowSums(idx > 0L)
  row_empty <- function(idx) rowSums(idx > 0L) == 0L
  direct_simple <- !hasAgonist & !hasAntagonist &
    row_count(ligand_idx) == 1L & row_count(receptor_idx) == 1L &
    row_empty(coA_idx) & row_empty(coI_idx)

  run_sparse_subset <- function(lr_use) {
    idx_list <- list(
      ligand = ligand_idx[lr_use, , drop = FALSE],
      receptor = receptor_idx[lr_use, , drop = FALSE],
      coA = coA_idx[lr_use, , drop = FALSE],
      coI = coI_idx[lr_use, , drop = FALSE],
      agonist = agonist_idx[lr_use, , drop = FALSE],
      antagonist = antagonist_idx[lr_use, , drop = FALSE]
    )
    if (subset_boot) {
      used_genes <- sort(unique(unlist(idx_list, use.names = FALSE)))
      used_genes <- used_genes[used_genes > 0L]
      if (!length(used_genes)) stop("No LR/cofactor genes are present in the signaling matrix.", call. = FALSE)
      old_to_new <- integer(n_genes_full)
      old_to_new[used_genes] <- seq_along(used_genes)
      remap_idx <- function(idx) {
        out <- idx
        positive <- out > 0L
        out[positive] <- old_to_new[out[positive]]
        out
      }
      data.kernel <- data.use[used_genes, , drop = FALSE]
      idx_list <- lapply(idx_list, remap_idx)
      n_kernel <- length(used_genes)
    } else {
      data.kernel <- data.use
      n_kernel <- n_genes_full
    }

    data.use.avg <- group_tri_mean_cpp(data.kernel, group_int, K)
    colnames(data.use.avg) <- levels(group)
    rownames(data.use.avg) <- rownames(data.kernel)
    avg_boot <- group_tri_mean_boot_cpp(data.kernel, group_boot, K)$avg_boot

    prob_fun <- if (algorithm %in% c("sparse_exact", "sparse_exact_subset_boot", "sparse_exact_ondemand_simple")) {
      cellchat_prob_from_avg_sparse_cpp
    } else {
      cellchat_prob_from_avg_cpp
    }
    res_part <- prob_fun(
      avg = data.use.avg,
      avgBoot = avg_boot,
      ligandIdx = idx_list$ligand,
      receptorIdx = idx_list$receptor,
      coAIdx = idx_list$coA,
      coIIdx = idx_list$coI,
      agonistIdx = idx_list$agonist,
      antagonistIdx = idx_list$antagonist,
      hasAgonist = hasAgonist[lr_use],
      hasAntagonist = hasAntagonist[lr_use],
      Kh = Kh,
      n_power = n
    )
    res_part$n_genes_kernel <- n_kernel
    res_part
  }

  if (sparse_stream) {
    data.use.avg <- group_tri_mean_cpp(data.use, group_int, K)
    colnames(data.use.avg) <- levels(group)
    rownames(data.use.avg) <- rownames(data.use)
    res <- cellchat_prob_sparse_stream_cpp(
      data = data.use,
      group = group_int,
      groupBoot = group_boot,
      avg = data.use.avg,
      ligandIdx = ligand_idx,
      receptorIdx = receptor_idx,
      coAIdx = coA_idx,
      coIIdx = coI_idx,
      agonistIdx = agonist_idx,
      antagonistIdx = antagonist_idx,
      hasAgonist = hasAgonist,
      hasAntagonist = hasAntagonist,
      Kh = Kh,
      n_power = n
    )
    res$n_genes_kernel <- n_genes_full
  } else if (ondemand_simple && any(direct_simple)) {
    nLR <- nrow(pairLRsig)
    Prob <- array(0, dim = c(K, K, nLR))
    Pval <- array(1, dim = c(K, K, nLR))
    simple_idx <- which(direct_simple)
    other_idx <- which(!direct_simple)

    ligand_gene <- as.integer(ligand_idx[simple_idx, 1])
    receptor_gene <- as.integer(receptor_idx[simple_idx, 1])
    res_simple <- cellchat_prob_simple_ondemand_cpp(
      data = data.use,
      group = group_int,
      groupBoot = group_boot,
      ligandGene = ligand_gene,
      receptorGene = receptor_gene,
      Kh = Kh,
      n_power = n
    )
    Prob[, , simple_idx] <- res_simple$prob
    Pval[, , simple_idx] <- res_simple$pval

    active_pairs <- res_simple$active_pairs %||% 0
    skipped_pairs <- res_simple$skipped_pairs %||% 0
    total_pairs <- res_simple$total_pairs %||% 0
    kernel_genes <- length(unique(c(ligand_gene, receptor_gene)))

    if (length(other_idx)) {
      res_other <- run_sparse_subset(other_idx)
      Prob[, , other_idx] <- res_other$prob
      Pval[, , other_idx] <- res_other$pval
      active_pairs <- active_pairs + (res_other$active_pairs %||% 0)
      skipped_pairs <- skipped_pairs + (res_other$skipped_pairs %||% 0)
      total_pairs <- total_pairs + (res_other$total_pairs %||% 0)
      kernel_genes <- kernel_genes + (res_other$n_genes_kernel %||% 0)
    }
    res <- list(
      prob = Prob,
      pval = Pval,
      active_pairs = active_pairs,
      skipped_pairs = skipped_pairs,
      total_pairs = total_pairs,
      active_fraction = if (total_pairs > 0) active_pairs / total_pairs else NA_real_,
      n_genes_kernel = kernel_genes,
      n_simple_ondemand_lr = length(simple_idx),
      n_other_lr = length(other_idx),
      boot_tri_mean_evals = res_simple$boot_tri_mean_evals %||% NA_real_
    )
  } else {
    res <- run_sparse_subset(rep(TRUE, nrow(pairLRsig)))
    n_genes_kernel <- res$n_genes_kernel %||% n_genes_full
  }

  Prob <- res$prob
  Pval <- res$pval
  dimnames(Prob) <- list(levels(group), levels(group), rownames(pairLRsig))
  dimnames(Pval) <- dimnames(Prob)
  object@net <- list(prob = Prob, pval = Pval)
  object@options$run.time <- as.numeric(Sys.time() - ptm, units = "secs")
  object@options$parameter <- list(type.mean = type, trim = NULL, raw.use = raw.use,
                                   population.size = population.size, nboot = nboot,
                                   seed.use = seed.use, Kh = Kh, n = n,
                                   distance.use = distance.use,
                                   interaction.range = NULL, ratio = NULL, tol = NULL,
                                   k.min = NULL, contact.dependent = FALSE,
                                   contact.range = NULL, contact.knn.k = NULL,
                                   contact.dependent.forced = FALSE)
  object@options$accelrcpp <- list(
    algorithm = algorithm,
    active_pairs = res$active_pairs %||% NA_real_,
    skipped_pairs = res$skipped_pairs %||% NA_real_,
    total_pairs = res$total_pairs %||% NA_real_,
    active_fraction = res$active_fraction %||% NA_real_,
    n_genes_full = n_genes_full,
    n_genes_kernel = res$n_genes_kernel %||% n_genes_kernel,
    n_simple_ondemand_lr = res$n_simple_ondemand_lr %||% NA_integer_,
    n_other_lr = res$n_other_lr %||% NA_integer_,
    boot_tri_mean_evals = res$boot_tri_mean_evals %||% NA_real_,
    cache_hits = res$cache_hits %||% NA_real_,
    cache_genes = res$cache_genes %||% NA_real_,
    cache_slots = res$cache_slots %||% NA_real_,
    streamed_lr = res$streamed_lr %||% NA_real_,
    max_lr_genes = res$max_lr_genes %||% NA_real_
  )
  object
}

computeCommunProbPathwayAccelRcpp <- function(object = NULL, net = NULL, pairLR.use = NULL, thresh = 0.05) {
  if (is.null(net)) net <- object@net
  if (is.null(pairLR.use)) pairLR.use <- object@LR$LRsig
  pathways <- unique(pairLR.use$pathway_name)
  pathway_id <- as.integer(factor(pairLR.use$pathway_name, levels = pathways))
  res <- pathway_sum_cpp(net$prob, net$pval, pathway_id, length(pathways), thresh)
  prob.pathways <- res$prob_pathway
  dimnames(prob.pathways) <- list(dimnames(net$prob)[[1]], dimnames(net$prob)[[2]], pathways)

  LR <- dimnames(net$prob)[[3]]
  LR.sig <- LR[res$lr_sum != 0]
  pathways.sig <- pathways[res$pathway_sum != 0]
  prob.pathways.sig <- prob.pathways[, , pathways.sig, drop = FALSE]
  if (length(pathways.sig)) {
    idx <- order(apply(prob.pathways.sig, 3, sum), decreasing = TRUE)
    pathways.sig <- pathways.sig[idx]
    prob.pathways.sig <- prob.pathways.sig[, , idx, drop = FALSE]
  }

  if (is.null(object)) {
    list(pathways = pathways.sig, prob = prob.pathways.sig)
  } else {
    object@net$LRs <- LR.sig
    object@netP$pathways <- pathways.sig
    object@netP$prob <- prob.pathways.sig
    object
  }
}

aggregateNetAccelRcpp <- function(object, thresh = 0.05) {
  res <- aggregate_net_cpp(object@net$prob, object@net$pval, thresh)
  dimnames(res$count) <- dimnames(object@net$prob)[1:2]
  dimnames(res$weight) <- dimnames(object@net$prob)[1:2]
  object@net$count <- res$count
  object@net$weight <- res$weight
  object
}

computeAveExprAccelRcpp <- function(object, features = NULL, group.by = NULL,
                                  type = "triMean", slot.name = c("data.signaling", "data"),
                                  data.use = NULL) {
  if (!identical(type, "triMean")) {
    stop("CellChatAccelRcpp currently supports only type = 'triMean'.", call. = FALSE)
  }
  slot.name <- match.arg(slot.name)
  if (is.null(data.use)) data.use <- methods::slot(object, slot.name)
  features.use <- if (is.null(features)) row.names(data.use) else intersect(features, row.names(data.use))
  data.use <- as.matrix(data.use[features.use, , drop = FALSE])
  labels <- if (is.null(group.by)) object@idents else object@meta[[group.by]]
  if (!is.factor(labels)) labels <- factor(labels)
  avg <- group_tri_mean_cpp(data.use, as.integer(labels), nlevels(labels))
  rownames(avg) <- features.use
  colnames(avg) <- levels(labels)
  avg
}
