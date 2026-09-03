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
#' @keywords internal
#' @useDynLib dppca, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom Rdpack reprompt
"_PACKAGE"
