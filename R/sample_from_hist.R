
#' Sample Synthetic Points From a DP Histogram
#'
#' Generates synthetic 2D samples by drawing from a differentially private
#' histogram. Each histogram bin is selected according to its (normalized)
#' probability mass, and the returned points are sampled uniformly within
#' the corresponding rectangular bin.
#'
#' @param X A data frame representing a 2D histogram. It must contain the
#'   following columns:
#'   \code{xmin}, \code{xmax}, \code{ymin}, \code{ymax}, \code{prob},
#'   where \code{prob} denotes the noisy bin probability.
#' @param k Integer; the number of synthetic points to generate.
#'
#' @details
#' The function first normalizes the (possibly unnormalized) histogram
#' probabilities to form a valid discrete distribution. Then, \code{k} bins
#' are drawn with replacement according to these probabilities. For each
#' selected bin, a point is sampled uniformly from the rectangle
#' \eqn{[x_{\min}, x_{\max}] \times [y_{\min}, y_{\max}]}.
#'
#' This allows the user to visualize or analyze the distribution implied by a
#' DP histogram without directly exposing individual (private) data points.
#'
#' @return A data frame with \code{k} rows and two columns:
#'   \code{x} and \code{y}, representing the sampled 2D points.
#'
#' @seealso \code{\link{dp_hist}}, \code{\link{dp_score_plot}}
#'
#' @export
sample_from_hist <- function(X, k) {
  # hist_df: columns must include xmin, xmax, ymin, ymax, prob
  stopifnot(all(c("xmin", "xmax", "ymin", "ymax", "prob") %in% names(X)))

  m <- nrow(X)
  if (m == 0) stop("hist_df must have at least one row.")

  q_hat <- X$prob
  if (any(q_hat < 0)) {
    stop("Histogram probabilities must be nonnegative.")
  }

  q_hat <- q_hat / sum(q_hat)

  selected_bins <- sample.int(m, size = k, replace = TRUE, prob = q_hat)

  xmin <- X$xmin[selected_bins]
  xmax <- X$xmax[selected_bins]
  ymin <- X$ymin[selected_bins]
  ymax <- X$ymax[selected_bins]

  # Uniform sampling for each bins
  x <- xmin + stats::runif(k) * (xmax - xmin)
  y <- ymin + stats::runif(k) * (ymax - ymin)

  out <- data.frame(x = x, y = y)
  rownames(out) <- NULL
  out
}


