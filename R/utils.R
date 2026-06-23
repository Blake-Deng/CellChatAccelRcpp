make_index_matrix <- function(items, data_genes, complex_input = NULL,
                              cofactor_input = NULL, cofactor_cols = NULL,
                              keep_duplicates = FALSE) {
  pos <- match(data_genes, data_genes)
  names(pos) <- data_genes
  sets <- vector("list", length(items))

  for (i in seq_along(items)) {
    x <- items[[i]]
    genes <- character()
    if (!is.na(x) && nzchar(x)) {
      if (!is.null(complex_input) && x %in% rownames(complex_input)) {
        subunits <- unlist(complex_input[x, grep("^subunit", colnames(complex_input)), drop = FALSE],
                           use.names = FALSE)
        genes <- subunits[nzchar(subunits)]
      } else if (!is.null(cofactor_input) && x %in% rownames(cofactor_input)) {
        subunits <- unlist(cofactor_input[x, cofactor_cols, drop = FALSE], use.names = FALSE)
        genes <- subunits[nzchar(subunits)]
      } else {
        genes <- x
      }
    }
    idx <- if (keep_duplicates) {
      unname(pos[genes])
    } else {
      unname(pos[intersect(genes, data_genes)])
    }
    sets[[i]] <- idx[!is.na(idx)]
  }

  max_len <- max(1L, max(lengths(sets)))
  out <- matrix(0L, nrow = length(items), ncol = max_len)
  for (i in seq_along(sets)) {
    if (length(sets[[i]])) out[i, seq_along(sets[[i]])] <- as.integer(sets[[i]])
  }
  out
}

check_supported_object <- function(object, type, population.size, distance.use) {
  if (!identical(type, "triMean")) {
    stop("CellChatFastCpp currently supports only type = 'triMean'.", call. = FALSE)
  }
  if (!identical(population.size, FALSE)) {
    stop("CellChatFastCpp currently supports only population.size = FALSE.", call. = FALSE)
  }
  if (!is.null(distance.use) && !identical(distance.use, FALSE)) {
    stop("CellChatFastCpp currently supports only non-spatial RNA workflows with distance.use = NULL/FALSE.", call. = FALSE)
  }
  if (!identical(object@options$datatype, "RNA")) {
    stop("CellChatFastCpp currently supports only object@options$datatype == 'RNA'.", call. = FALSE)
  }
}
