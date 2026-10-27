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
#' [dp_loading()], [dp_scree()], and [dp_score()] also
#' support print, summary and plot methods for their stored results.
#' For example, call `dp_loading()` with a privacy budget, then plot the fit.
#' The Shiny app also supports diagnostic comparisons that include non-private
#' data summaries.
#'
#' Data frames may include numeric and non-numeric columns. PCA uses the
#' numeric feature columns and warns when non-numeric columns are excluded.
#' Set `group` to a column name in [dp_score()] or [dp_pca()] to obtain grouped
#' score results; the selected group column is excluded from the PCA features,
#' even when it is numeric. Leave `group = NULL` for an ungrouped analysis.
#' Summaries identify the selected group variable and its groups. Numeric
#' features are used without automatic centering or standardization. See the
#' individual functions for preprocessing assumptions and privacy-budget
#' accounting.
#'
#' @keywords internal
#' @useDynLib dppca, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom Rdpack reprompt
"_PACKAGE"

# Common argument checks ----------------------------------------------------

#' Validate a logical scalar
#'
#' @param x Object to validate.
#' @param arg Argument name used in an error message.
#'
#' @return Invisibly returns `TRUE`.
#' @noRd
validate_logical_value <- function(x, arg) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop("`", arg, "` must be `TRUE` or `FALSE`.", call. = FALSE)
  }

  invisible(TRUE)
}

#' Validate a positive integer scalar
#'
#' @param x Object to validate.
#' @param arg Argument name used in an error message.
#'
#' @return Invisibly returns `TRUE`.
#' @noRd
validate_positive_integer <- function(x, arg) {
  if (
    !is.numeric(x) || length(x) != 1L || !is.finite(x) ||
    x < 1 || x != as.integer(x)
  ) {
    stop("`", arg, "` must be a positive integer.", call. = FALSE)
  }

  invisible(TRUE)
}
