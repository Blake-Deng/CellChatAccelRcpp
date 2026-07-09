\name{computeCommunProbAccelRcpp}
\alias{computeCommunProbAccelRcpp}
\title{Accelerated CellChat Communication Probability}
\usage{
computeCommunProbAccelRcpp(object, type = "triMean", raw.use = TRUE,
  population.size = FALSE, nboot = 100, seed.use = 1L, Kh = 0.5,
  n = 1, distance.use = NULL, algorithm = c("sparse_stream", "dense", "sparse_exact", "sparse_exact_subset_boot", "sparse_exact_ondemand_simple"))
}
\arguments{
\item{object}{A preprocessed single-dataset CellChat object.}
\item{type}{Mean type. Only \code{"triMean"} is currently supported.}
\item{raw.use}{Use \code{object@data.signaling}.}
\item{population.size}{Only \code{FALSE} is currently supported.}
\item{nboot}{Number of bootstrap permutations.}
\item{seed.use}{Random seed.}
\item{Kh}{Hill function parameter.}
\item{n}{Hill function exponent.}
\item{distance.use}{Spatial distance option. Only non-spatial workflows are supported.}
\item{algorithm}{Probability kernel. The default \code{"sparse_stream"} mode computes bootstrap triMean values lazily for genes and groups touched by active sender-receiver pairs, reusing cached values across ligand-receptor pairs. \code{"dense"} evaluates all sender-receiver group pairs. \code{"sparse_exact"} skips bootstrap comparisons for observed-zero group pairs whose final CellChat p-value is necessarily 1. \code{"sparse_exact_subset_boot"} additionally computes bootstrap group averages only for genes referenced by ligand-receptor and cofactor tables. \code{"sparse_exact_ondemand_simple"} is an experimental hybrid that computes bootstrap triMean values on demand for simple single-gene ligand-receptor pairs and falls back to the sparse subset path for complex pairs.}
}
\value{A CellChat object with \code{object@net$prob} and \code{object@net$pval}.}
\description{
Rcpp implementation of the common CellChat RNA \code{computeCommunProb}
workflow for \code{type = "triMean"}. The optional sparse algorithms preserve
dense output tensors while avoiding bootstrap work for observed-zero
sender-receiver pairs. The default \code{"sparse_stream"} mode is an
active-pair implementation that avoids dense bootstrap-array precomputation and
reuses lazy bootstrap triMean values for active ligand-receptor group pairs.
}
