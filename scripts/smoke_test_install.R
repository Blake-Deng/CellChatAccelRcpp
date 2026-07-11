#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(CellChatAccelRcpp)
})

expected_exports <- c(
  "computeAveExprAccelRcpp",
  "computeCommunProbAccelRcpp",
  "computeCommunProbPathwayAccelRcpp",
  "aggregateNetAccelRcpp"
)

missing_exports <- setdiff(expected_exports, getNamespaceExports("CellChatAccelRcpp"))
if (length(missing_exports) > 0L) {
  stop(
    "Missing expected exported functions: ",
    paste(missing_exports, collapse = ", "),
    call. = FALSE
  )
}

if (.Machine$sizeof.pointer < 8L) {
  stop("CellChatAccelRcpp requires a 64-bit R session for large-output workflows.", call. = FALSE)
}

cat("CellChatAccelRcpp smoke test OK\n")
cat("version:", as.character(utils::packageVersion("CellChatAccelRcpp")), "\n")
cat("R pointer size:", .Machine$sizeof.pointer, "bytes\n")
cat("exports:", paste(expected_exports, collapse = ", "), "\n")
