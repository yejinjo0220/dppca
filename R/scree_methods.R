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
  if (
    !is.numeric(k) || length(k) != 1 || !is.finite(k) ||
    k != as.integer(k) || k < 1 || k > d
  ) {
    stop("`k` must be an integer in {1, ..., ncol(X)}.", call. = FALSE)
  }
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
#' @return A list with components `scree` and `pve`.
#' @noRd
dp_scree_clipped <- function(X, k, eps, delta,
                             center = TRUE, standardize = FALSE,
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
    eps_pc <- eps / 2
    delta_pc <- delta / 2
    eps_scree <- eps / 2
    delta_scree <- delta / 2
  } else {
    eps_pc <- NULL
    delta_pc <- NULL
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
    eps = eps_pc,
    delta = delta_pc,
    center = center,
    standardize = standardize,
    cpp.option = cpp.option
  )

  Y <- X_proc %*% V_used

  eps_ell <- eps_scree / k
  delta_ell <- delta_scree / k

  scree_base <- numeric(k)
  scree <- numeric(k)

  for (ell in seq_len(k)) {
    y <- Y[, ell]
    ybar <- mean(y)
    w <- (y - ybar)^2

    w_clip <- pmin(w, C_clip)
    mu_hat <- mean(w_clip)

    scree_base[ell] <- (n / (n - 1)) * mu_hat

    Delta_ell <- C_clip / (n - 1)
    sd_noise <- Delta_ell * sqrt(2 * log(1.25 / delta_ell)) / eps_ell

    scree[ell] <- max(
      scree_base[ell] + stats::rnorm(1, mean = 0, sd = sd_noise),
      0
    )
  }

  if (isTRUE(mono)) {
    scree <- scree_post_processing(scree)
  }

  list(
    scree = scree,
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
#' @return A list with components `scree` and `pve`.
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
    eps_pc <- eps / 2
    delta_pc <- delta / 2
    eps_scree <- eps / 2
    delta_scree <- delta / 2
  } else {
    eps_pc <- NULL
    delta_pc <- NULL
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
    eps = eps_pc,
    delta = delta_pc,
    center = center,
    standardize = standardize,
    cpp.option = cpp.option
  )

  Y <- X_proc %*% V_used

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

    scree[ell] <- (n / (n - 1)) * max(muT, 0)
  }

  if (isTRUE(mono)) {
    scree <- scree_post_processing(scree)
  }

  list(
    scree = scree,
    pve = scree_to_pve(scree)
  )
}


#' Estimate a pure-DP upper quantile with a public lower bound
#'
#' Implements the one-sided unbounded quantile mechanism of Durfee (2023). For
#' data bounded below by the public value `lower`, the routine searches the
#' geometric grid `beta^i + lower - 1`, for `i = 1, ..., max_steps`, until a
#' noisy strict-below count reaches a noisy target count `q * length(x)`.
#'
#' The noisy target and every noisy query independently use one-sided
#' exponential noise with privacy parameter `epsilon / 2`. Because the count
#' queries are monotone and have sensitivity one, one call is pure
#' `epsilon`-DP under fixed-size replacement adjacency.
#'
#'
#' @param x Nonempty numeric vector containing only finite values.
#' @param q Quantile level in `(0, 1)`.
#' @param epsilon Positive total privacy parameter for one quantile release.
#' @param lower Finite public lower bound satisfying `x >= lower`. The default
#'   is `0`.
#' @param beta Geometric-grid base greater than `1`. The default is `1.01`.
#' @param max_steps Positive integer limiting the number of grid queries. The
#'   default is `5000L`.
#'
#' @return A numeric scalar giving the private quantile estimate. If no query
#'   crosses the noisy target within `max_steps`, the function warns and returns
#'   the final grid value.
#' @noRd
unbounded_quantile_upper <- function(x, q, epsilon, lower = 0,
                                     beta = 1.01, max_steps = 5000L) {
  x <- as.numeric(x)
  if (length(x) < 1L || anyNA(x) || any(!is.finite(x))) {
    stop("Durfee quantile input must contain finite values.", call. = FALSE)
  }
  if (
    !is.numeric(lower) || length(lower) != 1L || !is.finite(lower)
  ) {
    stop("`lower` must be a finite number.", call. = FALSE)
  }
  if (any(x < lower)) {
    stop("Durfee's one-sided estimator requires a valid public lower bound.",
         call. = FALSE)
  }
  if (
    !is.numeric(q) || length(q) != 1L ||
    !is.finite(q) || q <= 0 || q >= 1
  ) {
    stop("`q` must lie in (0, 1).", call. = FALSE)
  }
  if (
    !is.numeric(epsilon) || length(epsilon) != 1L ||
    !is.finite(epsilon) || epsilon <= 0
  ) {
    stop("Durfee quantile epsilon must be positive.", call. = FALSE)
  }
  if (
    !is.numeric(beta) || length(beta) != 1L ||
    !is.finite(beta) || beta <= 1
  ) {
    stop("`beta` must be greater than 1.", call. = FALSE)
  }
  if (
    !is.numeric(max_steps) || length(max_steps) != 1L ||
    !is.finite(max_steps) || max_steps < 1 ||
    max_steps > .Machine$integer.max || max_steps != floor(max_steps)
  ) {
    stop("`max_steps` must be a positive integer.", call. = FALSE)
  }
  max_steps <- as.integer(max_steps)

  # Monotone AboveThreshold with sensitivity one. Exponential noise and
  # epsilon_1 = epsilon_2 = epsilon / 2 give pure epsilon-DP under
  # fixed-size replacement adjacency.
  eps_threshold <- epsilon / 2
  eps_query <- epsilon / 2
  noisy_threshold <- q * length(x) + stats::rexp(1L, rate = eps_threshold)
  sorted_x <- sort(x)
  estimate <- lower
  halted <- FALSE

  for (i in seq_len(max_steps)) {
    estimate <- beta^i + lower - 1
    if (!is.finite(estimate)) {
      estimate <- .Machine$double.xmax
    }
    count_below <- findInterval(estimate, sorted_x, left.open = TRUE)
    noisy_count <- count_below + stats::rexp(1L, rate = eps_query)
    if (noisy_count >= noisy_threshold) {
      halted <- TRUE
      break
    }
  }
  if (!halted) {
    warning(
      "Durfee quantile search reached `max_steps`; consider a smaller `beta` only with a larger step limit, or increase the step limit.",
      call. = FALSE
    )
  }

  return(estimate)
}


#' Estimate a fully unbounded pure-DP quantile
#'
#' Implements the fully unbounded quantile mechanism in Algorithm 4 of Durfee
#' (2023). The routine applies `unbounded_quantile_upper()` once in the positive
#' direction at quantile `q` and once to the sign-flipped data at quantile
#' `1 - q`. A public zero-anchored transformation maps Algorithm 4 search step
#' `k` to step `k + 1` of the one-sided helper, so no public lower or upper data
#' bound is required.
#'
#' Each one-sided call receives `epsilon / 2` and internally divides that budget
#' equally between its noisy target and query counts. Sequential composition of
#' the two calls therefore gives pure `epsilon`-DP under fixed-size replacement
#' adjacency.
#'
#' @param x Nonempty numeric vector containing only finite values.
#' @param q Quantile level in `(0, 1)`.
#' @param epsilon Positive total privacy parameter for the two-sided release.
#' @param beta Geometric-grid base greater than `1`. The default is `1.01`.
#' @param max_steps Positive integer limiting each one-sided search. The default
#'   is `5000L`.
#'
#' @return A numeric scalar giving the private quantile estimate. If either
#'   one-sided search reaches `max_steps`, the function warns and does not select
#'   the capped result from that side.
#' @noRd
unbounded_quantile <- function(x, q, epsilon,
                               beta = 1.01, max_steps = 5000L) {
  x <- as.numeric(x)
  if (length(x) < 1L || anyNA(x) || any(!is.finite(x))) {
    stop("Durfee quantile input must contain finite values.", call. = FALSE)
  }
  if (
    !is.numeric(q) || length(q) != 1L ||
    !is.finite(q) || q <= 0 || q >= 1
  ) {
    stop("`q` must lie in (0, 1).", call. = FALSE)
  }
  if (
    !is.numeric(epsilon) || length(epsilon) != 1L ||
    !is.finite(epsilon) || epsilon <= 0
  ) {
    stop("Durfee quantile epsilon must be positive.", call. = FALSE)
  }
  if (
    !is.numeric(beta) || length(beta) != 1L ||
    !is.finite(beta) || beta <= 1
  ) {
    stop("`beta` must be greater than 1.", call. = FALSE)
  }
  if (
    !is.numeric(max_steps) || length(max_steps) != 1L ||
    !is.finite(max_steps) || max_steps < 1 ||
    max_steps > .Machine$integer.max || max_steps != floor(max_steps)
  ) {
    stop("`max_steps` must be a positive integer.", call. = FALSE)
  }
  max_steps <- as.integer(max_steps)

  # Map the zero-anchored Algorithm 4 query at step k to the one-sided
  # helper's query at step k + 1 while satisfying its public lower bound.
  positive_search_values <- pmin(
    pmax(beta * x + beta - 1, 0),
    .Machine$double.xmax
  )
  negative_search_values <- pmin(
    pmax(-beta * x + beta - 1, 0),
    .Machine$double.xmax
  )

  positive_hit_step_limit <- FALSE
  positive_raw_estimate <- withCallingHandlers(
    unbounded_quantile_upper(
      x = positive_search_values,
      q = q,
      epsilon = epsilon / 2,
      lower = 0,
      beta = beta,
      max_steps = max_steps
    ),
    warning = function(w) {
      if (grepl(
        "Durfee quantile search reached `max_steps`",
        conditionMessage(w),
        fixed = TRUE
      )) {
        positive_hit_step_limit <<- TRUE
        invokeRestart("muffleWarning")
      }
    }
  )

  negative_hit_step_limit <- FALSE
  negative_raw_estimate <- withCallingHandlers(
    unbounded_quantile_upper(
      x = negative_search_values,
      q = 1 - q,
      epsilon = epsilon / 2,
      lower = 0,
      beta = beta,
      max_steps = max_steps
    ),
    warning = function(w) {
      if (grepl(
        "Durfee quantile search reached `max_steps`",
        conditionMessage(w),
        fixed = TRUE
      )) {
        negative_hit_step_limit <<- TRUE
        invokeRestart("muffleWarning")
      }
    }
  )

  positive_step <- max(
    0,
    round(log1p(positive_raw_estimate) / log(beta)) - 1
  )
  negative_step <- max(
    0,
    round(log1p(negative_raw_estimate) / log(beta)) - 1
  )

  positive_estimate <- beta^positive_step - 1
  negative_estimate <- beta^negative_step - 1
  if (!is.finite(positive_estimate)) {
    positive_estimate <- .Machine$double.xmax
  }
  if (!is.finite(negative_estimate)) {
    negative_estimate <- .Machine$double.xmax
  }

  estimate <- if (!positive_hit_step_limit && positive_step > 0) {
    positive_estimate
  } else if (!negative_hit_step_limit && negative_step > 0) {
    -negative_estimate
  } else {
    0
  }

  if (positive_hit_step_limit || negative_hit_step_limit) {
    warning(
      "Durfee signed quantile search reached `max_steps`; consider a smaller `beta` only with a larger step limit, or increase the step limit.",
      call. = FALSE
    )
  }

  estimate
}


#' Estimate private scree values with private modified winsorized means
#'
#' Internal implementation of the private modified winsorized mean (PMWM) scree
#' estimator. The method preprocesses the data, computes non-private or private
#' principal component directions, privately estimates lower and upper
#' winsorization bounds for squared centered scores using a pure-DP
#' exponential-noise unbounded quantile routine, and releases a Gaussian-noised
#' winsorized mean for each component.
#'
#' If `split_mode = TRUE`, one subset is used for private quantile estimation
#' and the other subset is used for the winsorized mean step. If
#' `split_mode = FALSE`, the full sample is reused in both steps. The private
#' quantile releases consume only `epsilon`, while `delta` is reserved for the
#' Gaussian winsorized-mean release. If `mono = TRUE`, the final scree vector is
#' post-processed to be nonnegative and nonincreasing.
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
#' @param a Finite public lower post-processing bound for the private
#'   winsorization cutoffs.
#' @param b Finite public upper post-processing bound for the private
#'   winsorization cutoffs.
#' @param trim_const Positive constant controlling the practical clipping
#'   proportion `max(trim_const / n_q, eta)`.
#' @param eta Lower bound in the practical clipping proportion. Must lie in
#'   `[0, 0.5)`.
#' @param mono A logical value indicating whether to enforce a nonnegative and
#'   nonincreasing scree sequence by post-processing.
#' @return A list with components `scree` and `pve`.
#' @noRd
dp_scree_pmwm <- function(X, k, eps, delta,
                          a, b, trim_const, eta,
                          beta = 1.001,
                          g_dppca = FALSE, cpp.option = FALSE,
                          split_mode = TRUE,
                          center = TRUE, standardize = FALSE,
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
      !is.finite(eta) || eta < 0 || eta >= 0.5) {
    stop("eta must be in [0, 0.5).")
  }

  if (isTRUE(g_dppca)) {
    eps_pc <- eps / 2
    delta_pc <- delta / 2
    eps_scree <- eps / 2
    delta_scree <- delta / 2
  } else {
    eps_pc <- NULL
    delta_pc <- NULL
    eps_scree <- eps
    delta_scree <- delta
  }

  eps_ell <- eps_scree / k
  delta_ell <- delta_scree / k

  # Each fully unbounded private quantile receives eps_ell / 4. The helper
  # manages its internal split across the two one-sided searches. Because the
  # quantile releases are pure DP, they do not consume delta. The Gaussian
  # winsorized-mean release receives the remaining eps_ell / 2 and delta_ell.
  eps_Q <- eps_ell / 4
  eps_M <- eps_ell / 2
  delta_M <- delta_ell

  X_proc <- prep_matrix_for_pca(
    X = X,
    center = center,
    standardize = standardize
  )

  V_used <- dp_pc_dir(
    X = X,
    k = k,
    g_dppca = g_dppca,
    eps = eps_pc,
    delta = delta_pc,
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

  scree_base <- numeric(k)
  scree <- numeric(k)

  for (ell in seq_len(k)) {
    y <- Y[, ell]
    ybar <- mean(y)
    w <- (y - ybar)^2

    L <- unbounded_quantile(
      x = w[idx_q],
      q = trim_param,
      epsilon = eps_Q,
      beta = beta
    )

    U <- unbounded_quantile(
      x = w[idx_q],
      q = 1 - trim_param,
      epsilon = eps_Q,
      beta = beta
    )

    L <- min(max(L, a), b)
    U <- min(max(U, a), b)

    if (!is.finite(L) || !is.finite(U) || U < L) {
      L <- a
      U <- b
    }

    w_win <- pmin(pmax(w[idx_m], L), U)
    mu_hat <- mean(w_win)

    scree_base[ell] <- (n / (n - 1)) * mu_hat

    Delta_ell <- (n / (n - 1)) * (U - L) / n_m
    sd_noise <- Delta_ell * sqrt(2 * log(1.25 / delta_M)) / eps_M

    scree[ell] <- max(
      scree_base[ell] + stats::rnorm(1, mean = 0, sd = sd_noise),
      0
    )
  }

  if (isTRUE(mono)) {
    scree <- scree_post_processing(scree)
  }

  list(
    scree = scree,
    pve = scree_to_pve(scree)
  )
}
