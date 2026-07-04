\name{computeCommunProbAccelRcpp}
\alias{computeCommunProbAccelRcpp}
\title{Accelerated CellChat Communication Probability}
\usage{
computeCommunProbAccelRcpp(object, type = "triMean", raw.use = TRUE,
  population.size = FALSE, nboot = 100, seed.use = 1L, Kh = 0.5,
  n = 1, distance.use = NULL)
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
}
\value{A CellChat object with \code{object@net$prob} and \code{object@net$pval}.}
\description{
Rcpp implementation of the common CellChat RNA \code{computeCommunProb}
workflow for \code{type = "triMean"}.
}
