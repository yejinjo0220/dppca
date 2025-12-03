#' Differentially Private Plotting Frame for 2D Projections
#'
#' @description
#' Constructs a differentially private plotting frame (x- and y-limits) for
#' 2-dimensional data projected onto principal component axes. The frame is
#' created using differentially private estimates of the minimum and maximum
#' values based on smooth sensitivity.
#'
#' @param X A numeric matrix with exactly two columns representing the 2D
#'   projected data (e.g., PC1 and PC2 scores).
#' @param eps_q Privacy budget \eqn{\varepsilon_q} allocated for estimating the
#'   minimum and maximum (split equally between lower and upper quantiles).
#' @param delta_q Privacy parameter \eqn{\delta_q} for the smooth-sensitivity
#'   quantile estimator (also split equally).
#' @param inflate A non-negative numeric value indicating how much wider the
#'   DP frame should be relative to the DP min--max interval. A value of
#'   \code{inflate = 0.10} enlarges the interval by 10\%.
#' @param q Optional symmetric quantile level in \eqn{(0,1)}. If provided, the
#'   function uses the \eqn{q}- and \eqn{1-q}-quantiles instead of the extreme
#'   order statistics \eqn{1/n} and \eqn{(n-1)/n}. Useful for robustness when
#'   the data contain extreme outliers.
#' @param clip_to_data_range Logical; if \code{TRUE}, privatized quantiles are
#'   clipped to remain within the empirical data range.
#'
#' @details
#' To ensure differential privacy, the function first flattens the 2D input
#' matrix \code{X} into a vector and estimates privatized lower and upper
#' quantiles via \code{dp_quantile_ss()}, using the smooth sensitivity method
#' of Nissim, Raskhodnikova, and Smith (2007). The resulting DP min and max are
#' symmetrized to obtain a center point and half-length, which is optionally
#' expanded by the multiplicative factor \code{1 + inflate}.
#'
#' The returned \code{xlim} and \code{ylim} limits are identical, ensuring that
#' the plotting frame remains square—important for PCA score plots.
#'
#' @return A list with two components:
#' \itemize{
#'   \item \code{xlim}: A numeric vector of length 2 giving the DP x-axis limits.
#'   \item \code{ylim}: A numeric vector of length 2 giving the DP y-axis limits.
#' }
#'
#' @seealso \code{\link{dp_quantile_ss}}, \code{\link{dp_score_plot}},
#'   \code{\link{dp_score_plot_group}}
#'
#' @export
dp_frame <- function(X, eps_q, delta_q, inflate = 0.10,
                     q = NULL,
                     clip_to_data_range = TRUE) {

  stopifnot(is.matrix(X), ncol(X) == 2)
  z <- as.numeric(X)
  n <- length(z)
  stopifnot(n >= 2)

  if (is.null(q)) {
    q_min <- 1 / n
    q_max <- (n - 1) / n
  } else {
    stopifnot(is.numeric(q), length(q) == 1, q > 0, q < 1)

    q_max <- q
    q_min <- 1 - q

    if (q_min > q_max) {# If q < 0.5
      q_max_2 <- q_min
      q_min <- q_max
      q_max <- q_max_2
    }
  }
  eps_each   <- eps_q / 2
  delta_each <- delta_q / 2

  dp_min <- dp_quantile_ss(z, q = q_min, epsilon = eps_each, delta = delta_each)
  dp_max <- dp_quantile_ss(z, q = q_max, epsilon = eps_each, delta = delta_each)

  if (dp_min > dp_max) {
    dp_max_2 <- dp_min
    dp_min <- dp_max
    dp_max <- dp_max_2
  }

  center   <- (dp_min + dp_max) / 2
  half_len <- (dp_max - dp_min) / 2

  inflated_half_len <- (1 + inflate) * half_len

  frame_min <- center - inflated_half_len
  frame_max <- center + inflated_half_len

  list(xlim = c(frame_min, frame_max),
       ylim = c(frame_min, frame_max))
}

