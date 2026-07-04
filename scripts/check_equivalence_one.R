#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(CellChat)
  library(CellChatAccelRcpp)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) {
  stop("Usage: Rscript scripts/check_equivalence_one.R input.rds group_col [nboot]", call. = FALSE)
}

path <- args[[1]]
group_col <- args[[2]]
nboot <- if (length(args) >= 3L) as.integer(args[[3]]) else 5L

filter_obj <- function(obj, group_col, min_cells = 10) {
  group <- as.character(obj@meta.data[[group_col]])
  keep <- !is.na(group) & nzchar(group) & group != "openscpca-excluded"
  obj <- subset(obj, cells = colnames(obj)[keep])
  group <- as.character(obj@meta.data[[group_col]])
  valid <- names(table(group))[table(group) >= min_cells]
  subset(obj, cells = colnames(obj)[group %in% valid])
}

obj <- readRDS(path)
Seurat::DefaultAssay(obj) <- "RNA"
obj <- filter_obj(obj, group_col)
cc <- CellChat::createCellChat(object = obj, group.by = group_col, assay = "RNA")
cc@DB <- CellChatDB.human
cc <- CellChat::subsetData(cc)
cc <- CellChat::identifyOverExpressedGenes(cc)
cc <- CellChat::identifyOverExpressedInteractions(cc)

t0 <- proc.time()
cc_ref <- CellChat::computeCommunProb(cc, type = "triMean", nboot = nboot, seed.use = 1L)
t_ref <- proc.time() - t0

t0 <- proc.time()
cc_accel <- CellChatAccelRcpp::computeCommunProbAccelRcpp(cc, nboot = nboot, seed.use = 1L)
t_accel <- proc.time() - t0

cc_ref <- CellChat::computeCommunProbPathway(cc_ref)
cc_accel <- CellChatAccelRcpp::computeCommunProbPathwayAccelRcpp(cc_accel)
cc_ref <- CellChat::aggregateNet(cc_ref)
cc_accel <- CellChatAccelRcpp::aggregateNetAccelRcpp(cc_accel)

cat("reference elapsed:", unname(t_ref[["elapsed"]]), "\n")
cat("accelerated elapsed:", unname(t_accel[["elapsed"]]), "\n")
cat("prob max abs diff:", max(abs(cc_ref@net$prob - cc_accel@net$prob), na.rm = TRUE), "\n")
cat("pval max abs diff:", max(abs(cc_ref@net$pval - cc_accel@net$pval), na.rm = TRUE), "\n")
cat("pathway prob all.equal:", isTRUE(all.equal(cc_ref@netP$prob, cc_accel@netP$prob, tolerance = 1e-12)), "\n")
cat("aggregate count all.equal:", isTRUE(all.equal(cc_ref@net$count, cc_accel@net$count, tolerance = 1e-12)), "\n")
cat("aggregate weight all.equal:", isTRUE(all.equal(cc_ref@net$weight, cc_accel@net$weight, tolerance = 1e-12)), "\n")
