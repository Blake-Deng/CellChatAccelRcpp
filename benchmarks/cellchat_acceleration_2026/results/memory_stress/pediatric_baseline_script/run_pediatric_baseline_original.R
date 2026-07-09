suppressPackageStartupMessages(library(CellChat))
prepared <- "/home/dzf/cellchat_acceleration/CellChatAccelRcpp_e517ce2/benchmarks/cellchat_acceleration_2026/results/memory_stress/pediatric_sparse_stream_input/pediatric_library_id_cellchat_prepared_legal_lr.rds"
out <- "/home/dzf/cellchat_acceleration/CellChatAccelRcpp_e517ce2/benchmarks/cellchat_acceleration_2026/results/memory_stress/pediatric_library_id_baseline_original_120GB.rds"
tsv <- "/home/dzf/cellchat_acceleration/CellChatAccelRcpp_e517ce2/benchmarks/cellchat_acceleration_2026/results/memory_stress/pediatric_library_id_baseline_original_120GB.tsv"
obj <- readRDS(prepared)
cat("loaded prepared\n")
cat("data.signaling", paste(dim(obj@data.signaling), collapse=" x "), "groups", length(levels(obj@idents)), "LR", nrow(obj@LR$LRsig), "\n")
ptm <- Sys.time()
status <- "ok"
err <- ""
res <- tryCatch({
  obj <- CellChat::computeCommunProb(obj)
  compute_elapsed <- as.numeric(difftime(Sys.time(), ptm, units="secs"))
  cat("CellChat::computeCommunProb done elapsed", compute_elapsed, "\n")
  obj <- CellChat::filterCommunication(obj, min.cells = 10)
  obj <- CellChat::computeCommunProbPathway(obj)
  obj <- CellChat::aggregateNet(obj)
  saveRDS(obj, out)
  cat("saved", out, "\n")
  compute_elapsed
}, error=function(e) {
  status <<- "error"
  err <<- conditionMessage(e)
  cat("ERROR", err, "\n")
  NA_real_
})
row <- data.frame(
  engine="baseline_original",
  algorithm="CellChat::computeCommunProb",
  status=status,
  elapsed_sec=as.numeric(difftime(Sys.time(), ptm, units="secs")),
  compute_elapsed_sec=res,
  cells=ncol(obj@data.signaling),
  genes=nrow(obj@data.signaling),
  groups=length(levels(obj@idents)),
  lr=nrow(obj@LR$LRsig),
  error=err,
  stringsAsFactors=FALSE
)
write.table(row, tsv, sep="\t", quote=FALSE, row.names=FALSE)
print(row)
quit(save="no", status=ifelse(status == "ok", 0L, 2L))
