#' Differentially Private Quantile via Smooth Sensitivity
#'
#' @description
#' Computes a differentially private estimate of a quantile using the
#' smooth sensitivity framework of Nissim, Raskhodnikova, and Smith (2007).
#' This function is primarily used internally by \code{dp_frame()} to obtain
#' DP estimates of the minimum and maximum of the data.
#'
#' @param x Numeric vector of observations.
#' @param q Target quantile in \eqn{(0,1)}.
#' @param epsilon Privacy budget \eqn{\varepsilon}.
#' @param delta Privacy parameter \eqn{\delta}.
#' @param clip_to_data_range Logical; if \code{TRUE}, the returned noisy quantile
#'   is clipped to lie within \code{range(x)}.
#'
#' @details
#' The smooth sensitivity at scale \eqn{k} is computed by examining all
#' \eqn{(k+1)}-neighborhoods of the empirical quantile order statistic.
#' The global sensitivity decays as \eqn{\exp(-\beta k)}, where
#' \eqn{\beta = \varepsilon / (2 \log(1/\delta))}.
#' The resulting smooth sensitivity value determines the Laplace noise scale.
#'
#' If \code{clip_to_data_range = TRUE}, the privatized quantile estimate is
#' restricted to the empirical data range for numerical stability.
#'
#' @return A numeric value representing the differentially private estimate
#'   of the specified quantile.
#'
#' @keywords internal
#'
dp_quantile_ss <- function(x, q, epsilon, delta,
                           clip_to_data_range = TRUE) {

  stopifnot(length(x) >= 1, is.finite(epsilon), epsilon > 0,
            is.finite(delta), delta > 0, q > 0, q < 1)
  xs <- sort(x)
  n  <- length(xs)
  r  <- max(1, min(n, ceiling(q * n)))

  beta <- epsilon / (2 * log(1 / delta))

  ss_beta <- 0
  for (k in 0:(n - 1)) {
    t_vals <- 0:(k + 1)
    left   <- r + t_vals - (k + 1)
    right  <- r + t_vals
    valid  <- (left >= 1) & (right <= n)
    if (any(valid)) {
      A_k <- max(xs[right[valid]] - xs[left[valid]])
      val <- exp(-beta * k) * A_k
      if (val > ss_beta) {
        ss_beta <- val
      }
    }
  }

  b <- (2 * ss_beta) / epsilon
  q_noise <- xs[r] + VGAM::rlaplace(1, location = 0, scale = b)

  if (clip_to_data_range) {
    q_noise <- max(min(q_noise, max(xs)), min(xs))
  }

  q_noise
}
