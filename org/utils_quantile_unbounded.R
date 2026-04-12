#' Unbounded private upper quantile via log-binning
#'
#' Internal helper implementing the upper-tail unbounded quantile estimator
#' used by the PMWM scree estimator.
#'
#' @param data Numeric vector.
#' @param l Finite lower anchor for the transformed support.
#' @param b Log-binning base; must be greater than 1.
#' @param q Quantile level in (0,1).
#' @param eps_1,delta_1 Privacy parameters for the noisy threshold.
#' @param eps_2,delta_2 Privacy parameters for the noisy cumulative scan.
#' @param max_extra_bins Integer buffer added beyond the largest occupied bin.
#'
#' @return A numeric scalar quantile estimate.
#' @keywords internal
unbounded_quantile_upper <- function(
    data, l, b, q,
    eps_1, delta_1,
    eps_2, delta_2,
    max_extra_bins = 1000
) {
  data <- as.numeric(data)
  n <- length(data)

  if (n < 1) stop("data must have length >= 1.")
  if (!is.finite(l)) stop("l must be finite.")
  if (!is.finite(b) || b <= 1) stop("b must be > 1.")
  if (!is.finite(q) || q <= 0 || q >= 1) stop("q must be in (0,1).")

  if (!is.finite(eps_1) || eps_1 <= 0) stop("eps_1 must be > 0.")
  if (!is.finite(delta_1) || delta_1 <= 0 || delta_1 >= 1) {
    stop("delta_1 must be in (0,1).")
  }
  if (!is.finite(eps_2) || eps_2 <= 0) stop("eps_2 must be > 0.")
  if (!is.finite(delta_2) || delta_2 <= 0 || delta_2 >= 1) {
    stop("delta_2 must be in (0,1).")
  }

  max_extra_bins <- as.integer(max_extra_bins)
  if (!is.finite(max_extra_bins) || max_extra_bins < 0) {
    stop("max_extra_bins must be a nonnegative integer.")
  }

  # Gaussian noise for sensitivity-1 counts
  sd1 <- sqrt(2 * log(1.25 / delta_1)) / eps_1
  sd2 <- sqrt(2 * log(1.25 / delta_2)) / eps_2

  # log-binning indices i = floor(log_b(max(x - l + 1, b)))
  z <- pmax(data - l + 1, b)
  idx <- floor(log(z, base = b))
  idx <- as.integer(idx)

  counts <- table(idx)

  get_count <- function(i) {
    key <- as.character(i)
    if (key %in% names(counts)) as.integer(counts[[key]]) else 0L
  }

  t_noisy <- q * n + stats::rnorm(1, mean = 0, sd = sd1)

  cur <- 0L
  i <- 0L
  max_bin <- if (length(idx) > 0) max(idx) else 0L
  i_max <- max_bin + max_extra_bins

  repeat {
    cur <- cur + get_count(i)
    i <- i + 1L

    if (cur + stats::rnorm(1, mean = 0, sd = sd2) > t_noisy) {
      break
    }

    if (i > i_max) {
      break
    }
  }

  res <- tryCatch({
    val <- b^i + l - 1
    if (!is.finite(val)) .Machine$double.xmax else val
  }, error = function(e) {
    .Machine$double.xmax
  })

  res
}

#' Unbounded private quantile via log-binning
#'
#' Internal helper implementing a lower/upper quantile estimator by reducing
#' lower quantiles to upper quantiles on negated data.
#'
#' @param data Numeric vector.
#' @param l,u Finite bounds used to truncate the returned estimate.
#' @param b Log-binning base; must be greater than 1.
#' @param q Quantile level in (0,1).
#' @param eps_1,delta_1 Privacy parameters for the noisy threshold.
#' @param eps_2,delta_2 Privacy parameters for the noisy cumulative scan.
#' @param max_extra_bins Integer buffer added beyond the largest occupied bin.
#'
#' @return A numeric scalar quantile estimate truncated to \code{[l,u]}.
#' @keywords internal
unbounded_quantile <- function(
    data, l, u, b, q,
    eps_1, delta_1,
    eps_2, delta_2,
    max_extra_bins = 1000
) {
  data <- as.numeric(data)

  if (!is.finite(l) || !is.finite(u) || l > u) {
    stop("Need finite l <= u.")
  }
  if (!is.finite(q) || q <= 0 || q >= 1) {
    stop("q must be in (0,1).")
  }

  if (q < 0.5) {
    est <- -unbounded_quantile_upper(
      data = -data,
      l = -u,
      b = b,
      q = 1 - q,
      eps_1 = eps_1,
      delta_1 = delta_1,
      eps_2 = eps_2,
      delta_2 = delta_2,
      max_extra_bins = max_extra_bins
    )
    return(max(est, l))
  }

  est <- unbounded_quantile_upper(
    data = data,
    l = l,
    b = b,
    q = q,
    eps_1 = eps_1,
    delta_1 = delta_1,
    eps_2 = eps_2,
    delta_2 = delta_2,
    max_extra_bins = max_extra_bins
  )
  min(est, u)
}
