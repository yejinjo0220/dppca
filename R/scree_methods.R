# ============================================================
# scree_methods.R
# Internal method implementations and helpers for scree estimation
# ============================================================

#' Validate common scree inputs
#'
#' Checks common inputs used by scree-estimation routines. This helper validates
#' that `X` can be treated as a data matrix, that `k` is a valid number of
#' leading components, and that `eps` and `delta` are valid privacy
#' parameters.
#'
#' @param X Matrix-like data object with observations in rows and variables in
#'   columns.
#' @param k Number of leading principal components. Must be an integer between
#'   `1` and the number of columns in `X`.
#' @param eps Positive number defining the total `epsilon` privacy
#'   parameter for scree estimation.
#' @param delta Number in `(0, 1)` defining the total `delta` privacy
#'   parameter for scree estimation.
#'
#' @return Invisibly returns `TRUE` if all checks pass.
#' @noRd
validate_scree_inputs <- function(X, k, eps, delta) {
  X <- as.matrix(X)
  n <- nrow(X); d <- ncol(X)
  if (n < 2) stop("Need n >= 2.")
  if (d < 1) stop("Need ncol(X) >= 1.")
  if (!is.numeric(k) || length(k) != 1 || !is.finite(k) || k < 1 || k > d) stop("k must be an integer in {1, ..., ncol(X)}.")
  if (!is.numeric(eps) || length(eps) != 1 || !is.finite(eps) || eps <= 0) stop("eps must be a single positive number.")
  if (!is.numeric(delta) || length(delta) != 1 || !is.finite(delta) || delta <= 0 || delta >= 1) stop("delta must be a single number in (0, 1).")
  invisible(TRUE)
}

#' Post-process a scree sequence
#'
#' Truncates scree estimates at zero and, if needed, applies isotonic regression
#' to enforce a nonincreasing sequence. This is a valid post-processing step for
#' differentially private outputs because it does not access the original data.
#'
#' @param x Numeric vector of raw scree estimates.
#'
#' @return Numeric vector of nonnegative, nonincreasing scree estimates.
#' @noRd
scree_post_processing <- function(x) {
  y <- pmax(as.numeric(x), 0)
  fit <- stats::isoreg(seq_along(y), -y)
  pmax(-fit$yf, 0)
}

#' Convert scree values to proportions of variance explained
#'
#' Converts scree estimates to proportions of variance explained (PVE). If the
#' total scree is not positive and finite, a zero vector is returned.
#'
#' @param scree Numeric vector of scree estimates.
#'
#' @return Numeric vector with the same length as `scree`.
#' @noRd
scree_to_pve <- function(scree) {
  scree <- as.numeric(scree)
  s <- sum(scree)
  if (!is.finite(s) || s <= 0) return(rep(0, length(scree)))
  scree / s
}

#' Winsorize numeric values
#'
#' Clamps numeric values to the closed interval `[lo, hi]`.
#'
#' @param x Numeric vector.
#' @param lo Finite lower bound.
#' @param hi Finite upper bound satisfying `lo <= hi`.
#'
#' @return Numeric vector with all entries truncated to `[lo, hi]`.
#' @noRd
winsorization <- function(x, lo, hi) {
  if (!is.finite(lo) || !is.finite(hi) || lo > hi) stop("Need finite lo <= hi.")
  pmin(pmax(x, lo), hi)
}


#' Estimate private scree values with clipped means
#'
#' Internal implementation of the clipped-mean scree estimator. The method first
#' preprocesses the data, computes non-private or private principal component
#' directions, projects the data onto those directions, and then estimates the
#' variance of each score vector using a clipped mean with Gaussian noise.
#'
#' For component `ell`, the squared centered scores are clipped at `C_clip` and
#' the noisy clipped mean is rescaled by `n / (n - 1)` to match the usual sample
#' variance convention. If `mono = TRUE`, the final scree vector is
#' post-processed to be nonnegative and nonincreasing.
#'
#' @param X Numeric data matrix with observations in rows.
#' @param k Number of leading principal components.
#' @param eps Positive number defining the total `epsilon` privacy
#'   parameter for the scree routine.
#' @param delta Number in `(0, 1)` defining the total `delta` privacy
#'   parameter for the scree routine.
#' @param center A logical value indicating whether to center the columns of `X`
#'   before computing principal component directions.
#' @param standardize A logical value indicating whether to scale the columns of
#'   `X` by their sample standard deviations after optional centering.
#' @param C_clip Positive clipping threshold applied to squared centered scores.
#' @param g_dppca A logical value indicating whether to compute private
#'   principal component directions.
#' @param cpp.option A logical value passed to `dp_pc_dir()` when private
#'   directions are computed.
#' @param mono A logical value indicating whether to enforce a nonnegative and
#'   nonincreasing scree sequence by post-processing.
#'
#' @return A list with components `scree`, `scree_np`, and `pve`.
#' @noRd
dp_scree_clipped <- function(X, k, eps, delta,
                             center = TRUE, standardize = TRUE,
                             C_clip,
                             g_dppca = FALSE, cpp.option = FALSE,
                             mono = TRUE) {
  validate_scree_inputs(
    X = X,
    k = k,
    eps = eps,
    delta = delta
  )

  X <- as.matrix(X)
  n <- nrow(X)
  k <- as.integer(k)

  if (!is.numeric(C_clip) || length(C_clip) != 1 ||
      !is.finite(C_clip) || C_clip <= 0) {
    stop("C_clip must be a single positive number.")
  }

  if (isTRUE(g_dppca)) {
    eps_dir <- eps / 2
    delta_dir <- delta / 2
    eps_scree <- eps / 2
    delta_scree <- delta / 2
  } else {
    eps_dir <- NULL
    delta_dir <- NULL
    eps_scree <- eps
    delta_scree <- delta
  }

  X_proc <- prep_matrix_for_pca(
    X = X,
    center = center,
    standardize = standardize
  )

  V_used <- dp_pc_dir(
    X = X,
    k = k,
    g_dppca = g_dppca,
    eps_dir = eps_dir,
    delta_dir = delta_dir,
    center = center,
    standardize = standardize,
    cpp.option = cpp.option
  )

  Y <- X_proc %*% V_used

  eps_ell <- eps_scree / k
  delta_ell <- delta_scree / k

  scree_np <- numeric(k)
  scree <- numeric(k)

  for (ell in seq_len(k)) {
    y <- Y[, ell]
    ybar <- mean(y)
    w <- (y - ybar)^2

    w_clip <- pmin(w, C_clip)
    mu_hat <- mean(w_clip)

    scree_np[ell] <- (n / (n - 1)) * mu_hat

    Delta_ell <- C_clip / (n - 1)
    sd_noise <- Delta_ell * sqrt(2 * log(1.25 / delta_ell)) / eps_ell

    scree[ell] <- max(
      scree_np[ell] + stats::rnorm(1, mean = 0, sd = sd_noise),
      0
    )
  }

  if (isTRUE(mono)) {
    scree <- scree_post_processing(scree)
  }

  list(
    scree = scree,
    scree_np = scree_np,
    pve = scree_to_pve(scree)
  )
}


#' Select a private dyadic scale by a noisy histogram
#'
#' Internal helper for the Huber scree estimator. Nonnegative values are assigned
#' to dyadic bins indexed by powers of two, the bin counts are perturbed with
#' Laplace noise, and the dyadic scale corresponding to the largest noisy count
#' is returned.
#'
#' @param u Nonnegative numeric vector whose scale is summarized.
#' @param eps_m2 Positive number defining the `epsilon` privacy parameter for
#'   the noisy histogram.
#' @param k_min_m2 Integer lower bound for the dyadic bin index.
#' @param k_max_m2 Integer upper bound for the dyadic bin index.
#'
#' @return Positive numeric scalar giving the selected dyadic scale.
#' @noRd
dp_hist_m2 <- function(u, eps_m2, k_min_m2, k_max_m2) {
  u <- as.numeric(u)

  if (length(u) < 1) stop("u must have length >= 1.")
  if (!is.finite(eps_m2) || eps_m2 <= 0) stop("eps_m2 must be > 0.")
  if (missing(k_min_m2) || missing(k_max_m2)) {
    stop("k_min_m2 and k_max_m2 must be supplied.")
  }
  if (!is.numeric(k_min_m2) || length(k_min_m2) != 1 ||
      !is.finite(k_min_m2)) {
    stop("k_min_m2 must be a single finite number.")
  }
  if (!is.numeric(k_max_m2) || length(k_max_m2) != 1 ||
      !is.finite(k_max_m2)) {
    stop("k_max_m2 must be a single finite number.")
  }

  k_min_m2 <- as.integer(k_min_m2)
  k_max_m2 <- as.integer(k_max_m2)

  if (k_min_m2 > k_max_m2) stop("Need k_min_m2 <= k_max_m2.")

  bins <- k_min_m2:k_max_m2
  counts <- integer(length(bins))
  names(counts) <- as.character(bins)

  for (val in u) {
    if (!is.finite(val) || val < 0) next

    kk <- if (val == 0) {
      k_min_m2
    } else {
      as.integer(winsorization(floor(log(val, base = 2)), k_min_m2, k_max_m2))
    }

    counts[[as.character(kk)]] <- counts[[as.character(kk)]] + 1L
  }

  noisy <- as.numeric(counts) + VGAM::rlaplace(length(counts), scale = 1 / eps_m2)
  noisy <- pmax(noisy, 0)

  2^as.integer(names(counts)[which.max(noisy)])
}

#' Estimate a private scalar scale proxy
#'
#' Internal helper for the Huber scree estimator. The input is paired into
#' adjacent differences, converted to squared paired differences, summarized by
#' block medians, and passed to `dp_hist_m2()`.
#'
#' @param w Numeric vector, typically squared centered projected scores for one
#'   principal component.
#' @param eps_m2 Positive number defining the `epsilon` privacy parameter for
#'   the scale-proxy step.
#' @param k_min_m2 Integer lower bound for the dyadic bin index.
#' @param k_max_m2 Integer upper bound for the dyadic bin index.
#' @param M Optional number of blocks for the block-median step. If `NULL`, a
#'   default based on `sqrt(n / 2)` is used.
#'
#' @return Positive numeric scalar giving a private scale proxy.
#' @noRd
dp_m2 <- function(w, eps_m2, k_min_m2, k_max_m2, M = NULL) {
  w <- as.numeric(w)
  n <- length(w)

  if (n < 4) stop("Need length(w) >= 4.")
  if (!is.finite(eps_m2) || eps_m2 <= 0) stop("eps_m2 must be > 0.")
  if (missing(k_min_m2) || missing(k_max_m2)) {
    stop("k_min_m2 and k_max_m2 must be supplied.")
  }

  m_pairs <- floor(n / 2)
  w1 <- w[seq(1, 2 * m_pairs, by = 2)]
  w2 <- w[seq(2, 2 * m_pairs, by = 2)]

  s <- (w1 - w2)^2 / 2

  if (is.null(M)) M <- floor(sqrt(n / 2))
  M <- as.integer(M)

  if (M < 1) stop("M must be >= 1.")
  if (M > length(s)) M <- length(s)

  block_size <- floor(length(s) / M)
  if (block_size < 1) {
    M <- length(s)
    block_size <- 1
  }

  u <- numeric(M)
  for (b in seq_len(M)) {
    start <- (b - 1) * block_size + 1
    end <- if (b == M) length(s) else b * block_size
    u[b] <- stats::median(s[start:end])
  }

  dp_hist_m2(u = u, eps_m2 = eps_m2, k_min_m2 = k_min_m2, k_max_m2 = k_max_m2
  )
}

#' Convert a scale proxy to a Huber threshold
#'
#' Converts a private scale proxy into a robustification threshold for the scalar
#' Huber mean step.
#'
#' @param m2_hat Nonnegative private scale proxy.
#' @param eps_tau Positive number defining the `epsilon` privacy parameter for
#'   the Huber noisy-gradient-descent step for one component.
#' @param n_tau Effective sample size used in the threshold calculation.
#'
#' @return Positive numeric scalar Huber threshold.
#' @noRd
tau_from_m2 <- function(m2_hat, eps_tau, n_tau) {
  m2_hat <- max(as.numeric(m2_hat), 0)

  if (!is.finite(eps_tau) || eps_tau <= 0) stop("eps_tau must be > 0.")
  if (!is.numeric(n_tau) || length(n_tau) != 1 || n_tau < 2) {
    stop("n_tau must be >= 2.")
  }

  ln <- log(n_tau)
  denom <- sqrt((1 + ln) * ln)

  if (!is.finite(denom) || denom <= 0) denom <- 1

  sqrt(m2_hat) * sqrt((eps_tau * n_tau) / denom)
}

#' Estimate a private scalar mean by Huber noisy gradient descent
#'
#' Internal implementation of a scalar Huber-type private mean estimator. At
#' each iteration, residuals are clipped at `tau`, the average clipped residual
#' is used as a gradient step, and Gaussian noise is added for privacy.
#'
#' @param w Numeric vector representing a one-dimensional sample.
#' @param eps_gd Positive number defining the `epsilon` privacy parameter for
#'   noisy gradient descent.
#' @param delta_gd Number in `(0, 1)` defining the `delta` privacy parameter for
#'   noisy gradient descent.
#' @param tau Positive Huber threshold.
#' @param T Positive number of gradient-descent iterations.
#' @param mu0 Initial value for the iterative procedure.
#' @param eta0 Positive step size.
#'
#' @return Numeric scalar giving the final private estimate.
#' @noRd
dp_huber_noisy_gd <- function(w, eps_gd, delta_gd, tau, T, mu0 = 0, eta0 = 1) {
  w <- as.numeric(w)
  n <- length(w)

  if (n < 2) stop("Need length(w) >= 2.")
  if (!is.finite(eps_gd) || eps_gd <= 0) stop("eps_gd must be > 0.")
  if (!is.finite(delta_gd) || delta_gd <= 0 || delta_gd >= 1) {
    stop("delta_gd must be in (0, 1).")
  }
  if (!is.finite(tau) || tau <= 0) stop("tau must be > 0.")
  if (!is.finite(eta0) || eta0 <= 0) stop("eta0 must be > 0.")
  if (!is.numeric(T) || T < 1) stop("T must be >= 1.")

  T <- as.integer(T)
  mu <- as.numeric(mu0)

  eps_step <- eps_gd / T
  del_step <- delta_gd / T

  for (t in 0:(T - 1L)) {
    r <- w - mu
    psi <- winsorization(r, -tau, tau)
    g <- mean(psi)

    Delta_step <- (2 * eta0 * tau) / n
    sd_noise <- Delta_step * sqrt(2 * log(1.25 / del_step)) / eps_step

    mu <- mu + eta0 * g + stats::rnorm(1, mean = 0, sd = sd_noise)
  }

  mu
}


#' Estimate private scree values with Huber-type private means
#'
#' Internal implementation of the Huber scree estimator. The method preprocesses
#' the data, computes non-private or private principal component directions, and
#' estimates the variance of each score vector using a private Huber-type scalar
#' mean estimator with noisy gradient descent.
#'
#' For each component, a private scale proxy is first obtained with `dp_m2()`,
#' converted to a Huber threshold with `tau_from_m2()`, and then used in
#' `dp_huber_noisy_gd()`. The result is rescaled by `n / (n - 1)`. If
#' `mono = TRUE`, the final scree vector is post-processed to be nonnegative and
#' nonincreasing.
#'
#' @param X Numeric data matrix with observations in rows.
#' @param k Number of leading principal components.
#' @param eps Positive number defining the total `epsilon` privacy
#'   parameter for the scree routine.
#' @param delta Number in `(0, 1)` defining the total `delta` privacy
#'   parameter for the scree routine.
#' @param g_dppca A logical value indicating whether to compute private
#'   principal component directions.
#' @param cpp.option A logical value passed to `dp_pc_dir()` when private
#'   directions are computed.
#' @param center A logical value indicating whether to center the columns of `X`
#'   before computing principal component directions.
#' @param standardize A logical value indicating whether to scale the columns of
#'   `X` by their sample standard deviations after optional centering.
#' @param mu0 Initial value for noisy gradient descent.
#' @param eta0 Positive step size for noisy gradient descent.
#' @param T Optional number of gradient-descent iterations. If `NULL`, a default
#'   based on `ceiling(log(n))` is used.
#' @param M Optional number of blocks used in `dp_m2()`. If `NULL`, a default
#'   based on `floor(sqrt(n / 2))` is used.
#' @param k_min_m2 Integer lower bound for dyadic histogram bins used in
#'   `dp_hist_m2()`. This value must be supplied by the user.
#' @param k_max_m2 Integer upper bound for dyadic histogram bins used in
#'   `dp_hist_m2()`. This value must be supplied by the user.
#' @param m2_frac Fraction of the scree privacy parameter allocated to the
#'   private scale-proxy step. This value must be supplied by the user.
#' @param mono A logical value indicating whether to enforce a nonnegative and
#'   nonincreasing scree sequence by post-processing.
#'
#' @return A list with components `scree`, `scree_np`, and `pve`.
#' @noRd
dp_scree_huber <- function(X, k, eps, delta,
                           k_min_m2, k_max_m2, m2_frac,
                           g_dppca = FALSE, cpp.option = FALSE,
                           center = TRUE, standardize = FALSE,
                           mu0 = 0, eta0 = 1, T = NULL, M = NULL,
                           mono = TRUE) {
  validate_scree_inputs(X = X, k = k, eps = eps, delta = delta)

  X <- as.matrix(X)
  n <- nrow(X)
  k <- as.integer(k)

  if (missing(k_min_m2) || missing(k_max_m2) || missing(m2_frac)) {
    stop("k_min_m2, k_max_m2, and m2_frac must be supplied.")
  }
  if (!is.finite(eta0) || eta0 <= 0) {
    stop("eta0 must be > 0.")
  }
  if (!is.numeric(k_min_m2) || length(k_min_m2) != 1 ||
      !is.finite(k_min_m2)) {
    stop("k_min_m2 must be a single finite number.")
  }
  if (!is.numeric(k_max_m2) || length(k_max_m2) != 1 ||
      !is.finite(k_max_m2)) {
    stop("k_max_m2 must be a single finite number.")
  }
  if (as.integer(k_min_m2) > as.integer(k_max_m2)) {
    stop("Need k_min_m2 <= k_max_m2.")
  }
  if (!is.numeric(m2_frac) || length(m2_frac) != 1 ||
      !is.finite(m2_frac) || m2_frac <= 0 || m2_frac >= 1) {
    stop("m2_frac must be a single number in (0, 1).")
  }

  if (is.null(T)) T <- ceiling(log(n))
  T <- max(1L, as.integer(T))

  if (is.null(M)) M <- floor(sqrt(n) / 2)
  M <- max(1L, as.integer(M))

  if (isTRUE(g_dppca)) {
    eps_dir <- eps / 2
    delta_dir <- delta / 2
    eps_scree <- eps / 2
    delta_scree <- delta / 2
  } else {
    eps_dir <- NULL
    delta_dir <- NULL
    eps_scree <- eps
    delta_scree <- delta
  }

  eps_m2 <- eps_scree * m2_frac
  delta_m2 <- delta_scree * m2_frac

  eps_gd <- eps_scree * (1 - m2_frac)
  delta_gd <- delta_scree * (1 - m2_frac)

  eps_m2_ell <- eps_m2 / k
  delta_m2_ell <- delta_m2 / k

  eps_gd_ell <- eps_gd / k
  delta_gd_ell <- delta_gd / k

  X_proc <- prep_matrix_for_pca(X = X, center = center, standardize = standardize)

  V_used <- dp_pc_dir(
    X = X,
    k = k,
    g_dppca = g_dppca,
    eps_dir = eps_dir,
    delta_dir = delta_dir,
    center = center,
    standardize = standardize,
    cpp.option = cpp.option
  )

  Y <- X_proc %*% V_used

  scree_np <- numeric(k)
  scree <- numeric(k)

  for (ell in seq_len(k)) {
    y <- Y[, ell]
    ybar <- mean(y)
    w <- (y - ybar)^2

    m2_hat <- dp_m2(
      w = w,
      eps_m2 = eps_m2_ell,
      k_min_m2 = k_min_m2,
      k_max_m2 = k_max_m2,
      M = M
    )

    tau_ell <- tau_from_m2(
      m2_hat = m2_hat,
      eps_tau = eps_gd_ell,
      n_tau = n
    )

    if (!is.finite(tau_ell) || tau_ell <= 0) {
      tau_ell <- 1
    }

    muT <- dp_huber_noisy_gd(
      w = w,
      eps_gd = eps_gd_ell,
      delta_gd = delta_gd_ell,
      tau = tau_ell,
      T = T,
      mu0 = mu0,
      eta0 = eta0
    )

    scree_np[ell] <- (n / (n - 1)) * mean(w)
    scree[ell] <- (n / (n - 1)) * max(muT, 0)
  }

  if (isTRUE(mono)) {
    scree <- scree_post_processing(scree)
  }

  list(
    scree = scree,
    scree_np = scree_np,
    pve = scree_to_pve(scree)
  )
}


#' Estimate a private upper quantile by log-binning
#'
#' Internal implementation of an upper-tail private quantile estimator based on
#' a logarithmic grid, a noisy threshold, and noisy cumulative counts. It is used
#' by the PMWM scree estimator to construct private winsorization bounds.
#'
#' @param data Numeric vector.
#' @param l Finite lower anchor for the search grid.
#' @param beta Log-binning base for the geometric grid. Must be greater than
#'   `1`.
#' @param q Quantile level in `(0, 1)`.
#' @param eps_1 Positive number defining the `epsilon` privacy parameter for the
#'   noisy threshold.
#' @param delta_1 Number in `(0, 1)` defining the `delta` privacy parameter for
#'   the noisy threshold.
#' @param eps_2 Positive number defining the `epsilon` privacy parameter for the
#'   noisy cumulative scan.
#' @param delta_2 Number in `(0, 1)` defining the `delta` privacy parameter for
#'   the noisy cumulative scan.
#' @param max_extra_bins Nonnegative number of additional bins to search beyond
#'   the largest occupied bin.
#'
#' @return Numeric scalar giving a private upper-tail quantile estimate.
#' @noRd
unbounded_quantile_upper <- function(data, l, beta, q,
                                     eps_1, delta_1,
                                     eps_2, delta_2,
                                     max_extra_bins = 1000) {
  data <- as.numeric(data)
  n <- length(data)

  if (n < 1) stop("data must have length >= 1.")
  if (!is.finite(l)) stop("l must be finite.")
  if (!is.finite(beta) || beta <= 1) stop("beta must be > 1.")
  if (!is.finite(q) || q <= 0 || q >= 1) stop("q must be in (0, 1).")
  if (!is.finite(eps_1) || eps_1 <= 0) stop("eps_1 must be > 0.")
  if (!is.finite(delta_1) || delta_1 <= 0 || delta_1 >= 1) {
    stop("delta_1 must be in (0, 1).")
  }
  if (!is.finite(eps_2) || eps_2 <= 0) stop("eps_2 must be > 0.")
  if (!is.finite(delta_2) || delta_2 <= 0 || delta_2 >= 1) {
    stop("delta_2 must be in (0, 1).")
  }

  max_extra_bins <- as.integer(max_extra_bins)
  if (!is.finite(max_extra_bins) || max_extra_bins < 0) {
    stop("max_extra_bins must be a nonnegative integer.")
  }

  sd1 <- sqrt(2 * log(1.25 / delta_1)) / eps_1
  sd2 <- sqrt(2 * log(1.25 / delta_2)) / eps_2

  z <- pmax(data - l + 1, beta)
  idx <- as.integer(floor(log(z, base = beta)))
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

    if (cur + stats::rnorm(1, mean = 0, sd = sd2) > t_noisy) break
    if (i > i_max) break
  }

  val <- tryCatch(beta^i + l - 1, error = function(e) .Machine$double.xmax)
  if (!is.finite(val)) .Machine$double.xmax else val
}

#' Estimate a private quantile by log-binning
#'
#' Internal wrapper around `unbounded_quantile_upper()` for both lower- and
#' upper-tail quantiles on a finite target interval `[l, u]`. Lower-tail
#' quantiles are handled by sign-flipping the data and calling the upper-tail
#' routine.
#'
#' @param data Numeric vector.
#' @param l Finite lower truncation bound.
#' @param u Finite upper truncation bound.
#' @param beta Log-binning base for the geometric grid. Must be greater than
#'   `1`.
#' @param q Quantile level in `(0, 1)`.
#' @param eps_1 Positive number defining the `epsilon` privacy parameter for the
#'   noisy threshold.
#' @param delta_1 Number in `(0, 1)` defining the `delta` privacy parameter for
#'   the noisy threshold.
#' @param eps_2 Positive number defining the `epsilon` privacy parameter for the
#'   noisy cumulative scan.
#' @param delta_2 Number in `(0, 1)` defining the `delta` privacy parameter for
#'   the noisy cumulative scan.
#' @param max_extra_bins Nonnegative number of additional bins to search beyond
#'   the largest occupied bin.
#'
#' @return Numeric scalar giving a private quantile estimate.
#' @noRd
unbounded_quantile <- function(data, l, u, beta, q,
                               eps_1, delta_1,
                               eps_2, delta_2,
                               max_extra_bins = 1000) {
  data <- as.numeric(data)

  if (!is.finite(l) || !is.finite(u) || l > u) stop("Need finite l <= u.")
  if (!is.finite(beta) || beta <= 1) stop("beta must be > 1.")
  if (!is.finite(q) || q <= 0 || q >= 1) stop("q must be in (0, 1).")

  if (q < 0.5) {
    est <- -unbounded_quantile_upper(
      data = -data,
      l = -u,
      beta = beta,
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
    beta = beta,
    q = q,
    eps_1 = eps_1,
    delta_1 = delta_1,
    eps_2 = eps_2,
    delta_2 = delta_2,
    max_extra_bins = max_extra_bins
  )

  min(est, u)
}

#' Estimate private scree values with private modified winsorized means
#'
#' Internal implementation of the private modified winsorized mean (PMWM) scree
#' estimator. The method preprocesses the data, computes non-private or private
#' principal component directions, privately estimates lower and upper
#' winsorization bounds for squared centered scores, and releases a noisy
#' winsorized mean for each component.
#'
#' If `split_mode = TRUE`, one subset is used for private quantile estimation
#' and the other subset is used for the winsorized mean step. If
#' `split_mode = FALSE`, the full sample is reused in both steps. If
#' `mono = TRUE`, the final scree vector is post-processed to be nonnegative and
#' nonincreasing.
#'
#' @param X Numeric data matrix with observations in rows.
#' @param k Number of leading principal components.
#' @param eps Positive number defining the total `epsilon` privacy
#'   parameter for the scree routine.
#' @param delta Number in `(0, 1)` defining the total `delta` privacy
#'   parameter for the scree routine.
#' @param g_dppca A logical value indicating whether to compute private
#'   principal component directions.
#' @param cpp.option A logical value passed to `dp_pc_dir()` when private
#'   directions are computed.
#' @param split_mode A logical value indicating whether to split the sample into
#'   quantile and mean subsets.
#' @param center A logical value indicating whether to center the columns of `X`
#'   before computing principal component directions.
#' @param standardize A logical value indicating whether to scale the columns of
#'   `X` by their sample standard deviations after optional centering.
#' @param beta Log-binning base used by the private quantile estimator. Must be
#'   greater than `1`. The default is `1.001`.
#' @param a Finite lower support bound supplied to the private quantile routine.
#' @param b Finite upper support bound supplied to the private quantile routine.
#' @param trim_const Positive constant controlling the practical clipping
#'   proportion `max(trim_const / n_q, eta)`.
#' @param eta Lower bound in the practical clipping proportion. Must lie in
#'   `(0, 0.5)`.
#' @param mono A logical value indicating whether to enforce a nonnegative and
#'   nonincreasing scree sequence by post-processing.
#' @param max_extra_bins Nonnegative number of additional bins to search beyond
#'   the largest occupied bin.
#'
#' @return A list with components `scree`, `scree_np`, and `pve`.
#' @noRd
dp_scree_pmwm <- function(X, k, eps, delta,
                          a, b, trim_const, eta,
                          beta = 1.001,
                          g_dppca = FALSE, cpp.option = FALSE,
                          split_mode = TRUE,
                          center = TRUE, standardize = TRUE,
                          mono = TRUE, max_extra_bins = 1000) {
  validate_scree_inputs(
    X = X,
    k = k,
    eps = eps,
    delta = delta
  )

  X <- as.matrix(X)
  n <- nrow(X)
  k <- as.integer(k)

  if (missing(a) || missing(b) || missing(trim_const) || missing(eta)) {
    stop("a, b, trim_const, and eta must be supplied.")
  }
  if (!is.finite(beta) || beta <= 1) stop("beta must be > 1.")
  if (!is.finite(a) || !is.finite(b) || a > b) stop("Need finite a <= b.")
  if (!is.numeric(trim_const) || length(trim_const) != 1 ||
      !is.finite(trim_const) || trim_const <= 0) {
    stop("trim_const must be a single positive number.")
  }
  if (!is.numeric(eta) || length(eta) != 1 ||
      !is.finite(eta) || eta <= 0 || eta >= 0.5) {
    stop("eta must be in (0, 0.5).")
  }

  if (isTRUE(g_dppca)) {
    eps_dir <- eps / 2
    delta_dir <- delta / 2
    eps_scree <- eps / 2
    delta_scree <- delta / 2
  } else {
    eps_dir <- NULL
    delta_dir <- NULL
    eps_scree <- eps
    delta_scree <- delta
  }

  eps_ell <- eps_scree / k
  delta_ell <- delta_scree / k

  eps_Q <- eps_ell / 4
  delta_Q <- delta_ell / 4
  eps_M <- eps_ell / 2
  delta_M <- delta_ell / 2

  eps_q1 <- eps_Q / 2
  delta_q1 <- delta_Q / 2
  eps_q2 <- eps_Q / 2
  delta_q2 <- delta_Q / 2

  X_proc <- prep_matrix_for_pca(
    X = X,
    center = center,
    standardize = standardize
  )

  V_used <- dp_pc_dir(
    X = X,
    k = k,
    g_dppca = g_dppca,
    eps_dir = eps_dir,
    delta_dir = delta_dir,
    center = center,
    standardize = standardize,
    cpp.option = cpp.option
  )

  Y <- X_proc %*% V_used

  if (isTRUE(split_mode)) {
    m <- floor(n / 2)
    if (m < 1 || (n - m) < 1) {
      stop("split_mode = TRUE requires at least 2 observations.")
    }
    idx_q <- seq_len(m)
    idx_m <- seq.int(m + 1, n)
  } else {
    idx_q <- seq_len(n)
    idx_m <- seq_len(n)
  }

  n_q <- length(idx_q)
  n_m <- length(idx_m)

  trim_param <- min(max(trim_const / n_q, eta), 0.49)

  scree_np <- numeric(k)
  scree <- numeric(k)

  for (ell in seq_len(k)) {
    y <- Y[, ell]
    ybar <- mean(y)
    w <- (y - ybar)^2

    L <- unbounded_quantile(
      data = w[idx_q],
      l = a,
      u = b,
      beta = beta,
      q = trim_param,
      eps_1 = eps_q1,
      delta_1 = delta_q1,
      eps_2 = eps_q2,
      delta_2 = delta_q2,
      max_extra_bins = max_extra_bins
    )

    U <- unbounded_quantile(
      data = w[idx_q],
      l = a,
      u = b,
      beta = beta,
      q = 1 - trim_param,
      eps_1 = eps_q1,
      delta_1 = delta_q1,
      eps_2 = eps_q2,
      delta_2 = delta_q2,
      max_extra_bins = max_extra_bins
    )

    if (!is.finite(L) || !is.finite(U) || U < L) {
      L <- a
      U <- b
    }

    w_win <- pmin(pmax(w[idx_m], L), U)
    mu_hat <- mean(w_win)

    scree_np[ell] <- (n / (n - 1)) * mu_hat

    Delta_ell <- (n / (n - 1)) * (U - L) / n_m
    sd_noise <- Delta_ell * sqrt(2 * log(1.25 / delta_M)) / eps_M

    scree[ell] <- max(
      scree_np[ell] + stats::rnorm(1, mean = 0, sd = sd_noise),
      0
    )
  }

  if (isTRUE(mono)) {
    scree <- scree_post_processing(scree)
  }

  list(
    scree = scree,
    scree_np = scree_np,
    pve = scree_to_pve(scree)
  )
}
