#' dppca: Differentially Private Principal Component Analysis Visualization
#'
#' Tools for estimating principal component directions and constructing
#' differentially private PCA visualizations. The package provides private
#' principal component directions; clipped, modified-winsorized, and Huber
#' scree and proportion-of-variance-explained summaries; additive and sparse
#' two-dimensional score histograms; group-wise score visualizations; synthetic
#' score sampling from released private histograms; and an interactive Shiny
#' application.
#'
#' [dp_pca()] stores a shared set of private estimates for summaries and
#' scree, score, loading, variable, and biplot displays. Its default settings
#' disable sample-mean centering and omit non-private reference results.
#' Standalone plotting wrappers and the Shiny app also support diagnostic
#' comparisons that include non-private data summaries. See the individual
#' functions for preprocessing assumptions and privacy-budget accounting.
#'
#' @keywords internal
#' @useDynLib dppca, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom Rdpack reprompt
"_PACKAGE"
