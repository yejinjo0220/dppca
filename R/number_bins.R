
#' Recommend the Number of Histogram Bins for 2D DP Visualization
#'
#' @description
#' Computes a recommended number of bins per axis for 2-dimensional
#' differentially private histograms, using heuristics proposed in the
#' literature. The function supports two methods:
#' \itemize{
#'   \item \code{"W"}: Wasserman's rule-of-thumb for 2D histograms,
#'     yielding \eqn{m_{\text{axis}} \asymp n^{1/4}}.
#'   \item \code{"J"}: Jing's log-based rule,
#'     yielding \eqn{m_{\text{axis}} \asymp (n / \log n)^{1/3}}.
#' }
#'
#' @param X A data frame or matrix with \eqn{n} rows representing observations.
#'   Only the number of rows is used; the dimensionality of \code{X} is ignored,
#'   but the heuristics are intended for 2D visualizations.
#' @param method Character string specifying the bin selection method.
#'   Either \code{"W"} (Wasserman) or \code{"J"} (Jing).
#'
#' @details
#' For 2D histograms, if the total number of bins scales as
#' \eqn{m_{\text{total}} \asymp n^{1/2}}, then the number of bins per axis is
#' \eqn{m_{\text{axis}} \asymp n^{1/4}} (Wasserman, 2010).
#'
#' Alternatively, Jing's criterion derives a bandwidth
#' \eqn{h \asymp (\log n / n)^{1/3}}, implying that the number of bins per axis is
#' \eqn{m_{\text{axis}} \asymp (n / \log n)^{1/3}}.
#'
#' Because these formulas often underestimate the optimal visual smoothness,
#' practitioners frequently choose slightly larger values than the recommendation.
#'
#' @return An integer giving the suggested number of bins for each axis.
#'
#' @seealso \code{\link{dp_hist}}, \code{\link{dp_score_plot}}
#'
#' @export
number_bins <- function(X, method = c("W", "J")) {
  X <- as.data.frame(X)
  n <- nrow(X)
  if (n <= 1) stop("Number of rows must be greater than 1.")

  method <- match.arg(method)

  if (method == "W") {
    # Wasserman: d = 2 → m_total ~ n^(1/2), m_axis ~ n^(1/4)
    m_axis <- n^(1/4)

  } else if (method == "J") {
    # Jing: d = 2 → h ~ (log n / n)^(1/3), m_axis ~ (n / log n)^(1/3)
    if (n <= 2) stop("Log-based method requires n > 2.")
    m_axis <- (n / log(n))^(1/3)
  }

  as.integer(round(m_axis))
}
















