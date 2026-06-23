computeCommunProbFastCpp <- function(object, type = "triMean", raw.use = TRUE, population.size = FALSE,
                                     nboot = 100, seed.use = 1L, Kh = 0.5, n = 1,
                                     distance.use = NULL) {
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
  data.use.avg <- group_tri_mean_cpp(data.use, group_int, K)
  colnames(data.use.avg) <- levels(group)
  rownames(data.use.avg) <- rownames(data.use)

  set.seed(seed.use)
  permutation <- replicate(nboot, sample.int(nC, size = nC))
  group_boot <- matrix(group_int[permutation], nrow = nC, ncol = nboot)
  avg_boot <- group_tri_mean_boot_cpp(data.use, group_boot, K)$avg_boot

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

  res <- cellchat_prob_from_avg_cpp(
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
    Kh = Kh,
    n_power = n
  )

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
  object
}

computeCommunProbPathwayFastCpp <- function(object = NULL, net = NULL, pairLR.use = NULL, thresh = 0.05) {
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

aggregateNetFastCpp <- function(object, thresh = 0.05) {
  res <- aggregate_net_cpp(object@net$prob, object@net$pval, thresh)
  dimnames(res$count) <- dimnames(object@net$prob)[1:2]
  dimnames(res$weight) <- dimnames(object@net$prob)[1:2]
  object@net$count <- res$count
  object@net$weight <- res$weight
  object
}

computeAveExprFastCpp <- function(object, features = NULL, group.by = NULL,
                                  type = "triMean", slot.name = c("data.signaling", "data"),
                                  data.use = NULL) {
  if (!identical(type, "triMean")) {
    stop("CellChatFastCpp currently supports only type = 'triMean'.", call. = FALSE)
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
