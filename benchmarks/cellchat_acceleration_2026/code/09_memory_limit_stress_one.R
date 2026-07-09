#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default = NA_character_) {
  hit <- grep(paste0('^', name, '='), args, value = TRUE)
  if (!length(hit)) return(default)
  sub(paste0('^', name, '='), '', hit[[1]])
}
engine <- get_arg('--engine', 'baseline')
prepared <- get_arg('--prepared')
out <- get_arg('--out')
algorithm <- get_arg('--algorithm', 'sparse_stream')
if (!nzchar(prepared) || !file.exists(prepared)) stop('missing --prepared')
if (!nzchar(out)) stop('missing --out')
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)

start <- Sys.time()
status <- 'ok'
err <- ''
obj <- NULL
peak_notes <- Sys.getenv('MEMORY_LIMIT_LABEL', '')

safe_dim <- function(x) {
  d <- tryCatch(dim(x), error=function(e) NULL)
  if (is.null(d)) c(NA_integer_, NA_integer_) else d
}
scalar_or_na <- function(x, na_value) {
  if (is.null(x) || !length(x)) return(na_value)
  x[[1]]
}

tryCatch({
  suppressPackageStartupMessages({
    library(CellChat)
    library(CellChatAccelRcpp)
  })
  base <- readRDS(prepared)
  dims <- safe_dim(base@data.signaling)
  cells <- dims[[2]]
  genes <- dims[[1]]
  groups <- length(levels(base@idents))
  lr <- nrow(base@LR$LRsig)
  if (identical(engine, 'baseline')) {
    obj <- CellChat::computeCommunProb(base)
    obj <- CellChat::filterCommunication(obj, min.cells = 10)
    obj <- CellChat::computeCommunProbPathway(obj)
    obj <- CellChat::aggregateNet(obj)
  } else if (identical(engine, 'sparse_stream')) {
    obj <- CellChatAccelRcpp::computeCommunProbAccelRcpp(base, algorithm = algorithm)
    obj <- CellChat::filterCommunication(obj, min.cells = 10)
    obj <- CellChatAccelRcpp::computeCommunProbPathwayAccelRcpp(obj)
    obj <- CellChatAccelRcpp::aggregateNetAccelRcpp(obj)
  } else {
    stop('unsupported engine: ', engine)
  }
}, error=function(e) {
  status <<- 'error'
  err <<- conditionMessage(e)
})

elapsed <- as.numeric(difftime(Sys.time(), start, units='secs'))
if (exists('base')) {
  dims <- safe_dim(base@data.signaling)
  cells <- dims[[2]]
  genes <- dims[[1]]
  groups <- length(levels(base@idents))
  lr <- nrow(base@LR$LRsig)
} else {
  cells <- genes <- groups <- lr <- NA_integer_
}
accel_algorithm <- if (!is.null(obj)) {
  scalar_or_na(tryCatch(obj@options$accelrcpp$algorithm, error=function(e) NA_character_), NA_character_)
} else NA_character_
active_fraction <- if (!is.null(obj)) {
  scalar_or_na(tryCatch(obj@options$accelrcpp$active_fraction, error=function(e) NA_real_), NA_real_)
} else NA_real_
row <- data.frame(
  engine=engine,
  algorithm=ifelse(engine == 'sparse_stream', algorithm, ''),
  memory_limit=peak_notes,
  status=status,
  elapsed_sec=elapsed,
  cells=cells,
  genes=genes,
  groups=groups,
  lr=lr,
  accel_algorithm=accel_algorithm,
  active_fraction=active_fraction,
  error=err,
  stringsAsFactors=FALSE
)
write.table(row, out, sep='\t', quote=FALSE, row.names=FALSE)
print(row)
if (!is.null(obj) && identical(status, 'ok')) {
  saveRDS(obj, sub('\\.tsv$', '.computed.rds', out))
}
quit(save='no', status=ifelse(status == 'ok', 0L, 2L))
