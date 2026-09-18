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
  if (!is.numeric(X) || is.complex(X) ||
      any(!is.finite(X))) {
    stop(
      "X must contain finite real numeric values.",
      call. = FALSE
    )
  }
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

#' Squared differences of disjoint score pairs
#'
#' Randomly permutes the rows of `Y` and returns
#' \eqn{W_{jl}=(Y_{b_j,l}-Y_{a_j,l})^2/2} for
#' \eqn{m=\lfloor n/2\rfloor} disjoint pairs. One row is unused when `n`
#' is odd. The permutation is independent of the score values.
#'
#' @param Y Projected score matrix with observations in rows.
#' @return An `m` by `ncol(Y)` matrix of squared pair differences.
#' @noRd
paired_score_squares <- function(Y) {
  Y <- as.matrix(Y)

  m <- nrow(Y) %/% 2L
  idx <- sample.int(nrow(Y))

  a <- idx[2L * seq_len(m) - 1L]
  b <- idx[2L * seq_len(m)]

  ( (Y[b, , drop = FALSE] - Y[a, , drop = FALSE]) / sqrt(2) )^2
}


#' Estimate private scree values with clipped means
#'
#' Uses `paired_score_squares()` to form squared differences from disjoint
#' score pairs, clips them at `C_clip`, and perturbs the vector of means.
#' For `m = floor(nrow(Y) / 2)` and `k = ncol(Y)`, Gaussian noise has standard
#' deviation `sqrt(k) * C_clip / (m * gamma)`, where `gamma` is obtained from
#' the method's epsilon and delta by `mu_from_eps_delta()`. There is no
#' additional sample-size rescaling. See [dp_scree()] for the assumptions
#' underlying the sensitivity calculation.
#'
#' @param Y Projected score matrix with at least two observations in rows and
#'   selected private principal components in columns.
#' @param eps Positive epsilon budget for the joint scree-vector release,
#'   after the shared loading allocation.
#' @param delta Delta budget in `(0, 1)` for the joint scree-vector release.
#' @param C_clip Positive public clipping threshold for squared pair differences.
#' @param mono Whether to enforce a nonnegative, nonincreasing final sequence.
#'
#' @return A list with `scree` and `pve` vectors.
#' @noRd
dp_scree_clipped <- function(
    Y, eps, delta, C_clip, mono = TRUE
) {
  Y <- as.matrix(Y)
  n <- nrow(Y)
  k <- ncol(Y)

  validate_scree_inputs(
    X = Y,
    k = k,
    eps = eps,
    delta = delta
  )

  n_pairs <- n %/% 2L

  if (!is.numeric(C_clip) || length(C_clip) != 1 ||
      !is.finite(C_clip) || C_clip <= 0) {
    stop("C_clip must be a single positive number.")
  }

  eps_scree <- eps
  delta_scree <- delta

  W <- paired_score_squares(Y)

  scree_base <- colMeans(pmin(W, C_clip))

  gamma_C <- mu_from_eps_delta(eps_scree, delta_scree)
  sd_noise <- sqrt(k) * C_clip / (n_pairs * gamma_C)

  scree <- pmax(
    scree_base + stats::rnorm(k, mean = 0, sd = sd_noise),
    0
  )

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

  noisy <- as.numeric(counts) + VGAM::rlaplace(length(counts), scale = 2 / eps_m2)
  noisy <- pmax(noisy, 0)

  2^as.integer(names(counts)[which.max(noisy)])
}

#' Estimate a private scalar scale proxy
#'
#' Randomly permutes `w`, forms squared differences from disjoint pairs of
#' its entries, summarizes them by block medians, and calls `dp_hist_m2()`.
#' For Huber scree estimation, `w` already contains squared score-pair
#' differences; this helper performs an additional pairing for scale estimation.
#'
#' @param w Numeric vector of length at least four; in scree estimation, a
#'   column returned by `paired_score_squares()`.
#' @param eps_m2 Positive epsilon budget for the pure-DP scale-proxy step.
#' @param k_min_m2 Integer lower dyadic-bin index.
#' @param k_max_m2 Integer upper dyadic-bin index.
#' @param M Optional number of blocks. If `NULL`, uses
#'   `floor(sqrt(length(w) / 2))`; values above the number of available
#'   input pairs are reduced to that number.
#'
#' @return A positive scalar giving the private scale proxy.
#' @noRd
dp_m2 <- function(w, eps_m2, k_min_m2, k_max_m2, M = NULL) {
  w <- as.numeric(w)
  n <- length(w)
  w <- w[sample.int(n)]  # permute the sample

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
#' @param eps_tau Positive threshold-tuning parameter. Huber scree estimation
#'   passes the joint gradient-descent epsilon budget divided by `k`; this
#'   calculation does not itself release data or spend an additional budget.
#' @param n_tau Effective sample size; the number of score pairs for scree.
#'
#' @return Nonnegative numeric scalar Huber threshold.
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


#' Convert a Gaussian-DP parameter to delta
#'
#' Evaluates the trade-off conversion from a positive Gaussian-DP parameter
#' `mu` and a positive `epsilon` value to the corresponding nonnegative
#' `delta` value.
#'
#' @param mu Positive Gaussian-DP parameter.
#' @param eps Positive `epsilon` privacy parameter.
#'
#' @return A nonnegative numeric scalar giving the corresponding `delta` value.
#' @noRd
gdp_delta <- function(mu, eps) {
  z1 <- -eps / mu + mu / 2
  z2 <- -eps / mu - mu / 2

  p1 <- stats::pnorm(z1)
  p2 <- exp(eps + stats::pnorm(z2, log.p = TRUE))

  max(p1 - p2, 0)
}

#' Recover a Gaussian-DP parameter from epsilon and delta
#'
#' Numerically inverts `gdp_delta()` on a log scale.
#'
#' @param eps Positive `epsilon` privacy parameter.
#' @param delta Number in `(0, 1)` defining the target `delta` privacy
#'   parameter.
#' @param tol Positive root-finding tolerance passed to [stats::uniroot()].
#'
#' @return A positive numeric scalar giving the Gaussian-DP parameter `mu`.
#' @noRd
mu_from_eps_delta <- function(eps, delta, tol = 1e-12) {
  if (!is.finite(eps) || eps <= 0) {
    stop("eps must be > 0.")
  }
  if (!is.finite(delta) || delta <= 0 || delta >= 1) {
    stop("delta must be in (0, 1).")
  }

  log_mu <- stats::uniroot(
    function(log_mu) {
      gdp_delta(exp(log_mu), eps) - delta
    },
    interval = c(log(1e-12), log(1e3)),
    tol = tol
  )$root

  exp(log_mu)
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
  w <- w[sample.int(n)]  # permute the sample

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

  # privacy budget allocation based on GDP <-> (eps, delta)-DP transformation
  gdp_mu <- mu_from_eps_delta(eps_gd, delta_gd)
  # eps_step <- eps_gd / iter_tot
  # del_step <- delta_gd / iter_tot
  gdp_mu_step <- gdp_mu / sqrt(T)

  for (t in 0:(T - 1L)) {
    r <- w - mu
    psi <- winsorization(r, -tau, tau)
    g <- mean(psi)

    Delta_step <- (2 * eta0 * tau) / n
    # sd_noise <- Delta_step * sqrt(2 * log(1.25 / del_step)) / eps_step
    sd_noise <- Delta_step  / gdp_mu_step

    mu <- mu + eta0 * g + stats::rnorm(1, mean = 0, sd = sd_noise)
  }

  mu
}


#' Estimate private scree values with Huber-type private means
#'
#' Forms `m = floor(nrow(Y) / 2)` squared score-pair differences per component.
#' A private scale from `dp_m2()` is converted to a threshold by
#' `tau_from_m2()`. All components are then updated together by the noisy
#' Huber-gradient loop in this function; `dp_huber_noisy_gd()` is not called.
#' At least four pairs, or eight observations, are required.
#'
#' The pure-DP scale step uses `m2_frac * eps / k` per component. The gradient
#' step uses `(1 - m2_frac) * eps` and all delta. Component `l` receives
#' Gaussian noise with standard deviation
#' `2 * eta0 * tau[l] * sqrt(k * T) / (m * gamma_gd)` at each iteration.
#' Each iterate is projected to nonnegative values. There is no additional
#' sample-size rescaling; optional monotone adjustment is applied last.
#' See [dp_scree()] for preprocessing and accounting conditions.
#'
#' @param Y Projected score matrix with observations in rows and selected
#'   private principal components in columns.
#' @param eps Positive total epsilon budget for this scree method after
#'   the shared loading allocation.
#' @param delta Delta budget in `(0, 1)`, used entirely by gradient descent.
#' @param k_min_m2,k_max_m2 Public dyadic-bin bounds used by `dp_hist_m2()`.
#' @param m2_frac Fraction in `(0, 1)` of epsilon used for scale estimation.
#' @param mu0 Finite public initial value for noisy gradient descent.
#' @param eta0 Positive fixed step size.
#' @param T Optional iteration count; defaults to `ceiling(log(m))`.
#' @param M Optional scale block count; defaults to `floor(sqrt(m / 2))`.
#' @param mono Whether to enforce a nonnegative, nonincreasing final sequence.
#'
#' @return A list with `scree` and `pve` vectors.
#' @noRd
dp_scree_huber <- function(
    Y, eps, delta,
    k_min_m2, k_max_m2, m2_frac,
    mu0 = 0, eta0 = 1, T = NULL, M = NULL,
    mono = TRUE
) {
  Y <- as.matrix(Y)
  n <- nrow(Y)
  k <- ncol(Y)

  validate_scree_inputs(
    X = Y,
    k = k,
    eps = eps,
    delta = delta
  )

  n_pairs <- n %/% 2L

  if (n_pairs < 4L) {
    stop("Huber scale estimation requires at least 4 pairs (n >= 8).")
  }
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

  if (is.null(T)) T <- ceiling(log(n_pairs))
  T <- max(1L, as.integer(T))

  if (is.null(M)) M <- floor(sqrt(n_pairs / 2))
  M <- max(1L, as.integer(M))

  eps_scree <- eps
  delta_scree <- delta

  eps_m2_ell <- m2_frac * eps_scree / k

  eps_gd <- (1 - m2_frac) * eps_scree
  delta_gd <- delta_scree

  eps_tau_tune <- eps_gd / k

  W <- paired_score_squares(Y)

  if (length(T) != 1L || !is.finite(T) || T < 1L) {
    stop("T must be one positive integer.")
  }
  if (length(mu0) != 1L || !is.finite(mu0)) {
    stop("mu0 must be one finite public value.")
  }
  if (length(eta0) != 1L || !is.finite(eta0) || eta0 <= 0) {
    stop("eta0 must be one finite positive value.")
  }

  tau <- numeric(k)
  for (ell in seq_len(k)) {
    m2_hat <- dp_m2(
      w = W[, ell],
      eps_m2 = eps_m2_ell,
      k_min_m2 = k_min_m2,
      k_max_m2 = k_max_m2,
      M = M
    )

    tau[ell] <- tau_from_m2(
      m2_hat = m2_hat,
      eps_tau = eps_tau_tune,
      n_tau = n_pairs
    )
  }

  if (any(!is.finite(tau)) || any(tau <= 0)) {
    stop("Huber thresholds must be finite and positive.")
  }

  gamma_gd <- mu_from_eps_delta(eps_gd, delta_gd)

  sd_noise <- 2 * eta0 * tau * sqrt(k * T) / (n_pairs * gamma_gd)

  if (any(!is.finite(sd_noise))) {
    stop("Non-finite Gaussian noise scale.")
  }

  theta <- rep(mu0, k)

  for (t in seq_len(T)) {
    g <- vapply(seq_len(k), function(ell) {
      r <- W[, ell] - theta[ell]
      mean(pmin(pmax(r, -tau[ell]), tau[ell]))
    }, numeric(1))

    theta <- pmax(
      theta + eta0 * g +
        stats::rnorm(k, mean = 0, sd = sd_noise),
      0
    )
  }

  scree <- theta

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
#' Forms squared differences from `m = floor(nrow(Y) / 2)` disjoint score
#' pairs. For each component, lower and upper cutoffs are estimated with
#' `unbounded_quantile_upper()`, using public lower bound zero, then
#' truncated to `[a, b]`. The pair values are winsorized to those cutoffs.
#'
#' Each of the `2 * k` pure-DP quantiles receives `eps / (4 * k)`. The joint
#' winsorized-mean release uses `eps / 2` and all delta. Gaussian noise for
#' component `l` has standard deviation
#' `sqrt(k) * (U[l] - L[l]) / (n_m * gamma_M)`, where `n_m` is the number of
#' pair rows used for the mean. There is no additional sample-size rescaling.
#' See [dp_scree()] for preprocessing and accounting conditions.
#'
#' With `split_mode = TRUE`, pair rows are randomly split between quantile
#' and mean estimation, requiring at least two pairs or four observations.
#' The default `FALSE` reuses all pair rows in both steps. Negative final
#' estimates are truncated at zero, with optional monotone adjustment.
#'
#' @param Y Projected score matrix with at least two observations in rows and
#'   selected private principal components in columns.
#' @param eps Positive epsilon budget for this method after the shared
#'   loading allocation; divided between quantiles and the joint mean release.
#' @param delta Delta budget in `(0, 1)`, used entirely by the joint mean release.
#' @param a,b Finite public post-processing bounds for private cutoffs.
#' @param trim_const Positive public constant in the trimming proportion
#'   `min(max(trim_const / n_q, eta), 0.49)`, where `n_q` counts quantile pairs.
#' @param eta Public lower bound on the trimming proportion, in `[0, 0.5)`.
#' @param beta Geometric search-grid base greater than one. Default `1.001`.
#' @param split_mode Whether to split pair rows into quantile and mean subsets.
#'   Default `FALSE`.
#' @param mono Whether to enforce a nonnegative, nonincreasing final sequence.
#'
#' @return A list with `scree` and `pve` vectors.
#' @noRd
dp_scree_pmwm <- function(
    Y, eps, delta,
    a, b, trim_const, eta,
    beta = 1.001,
    split_mode = FALSE,
    mono = TRUE
) {
  Y <- as.matrix(Y)
  n <- nrow(Y)
  k <- ncol(Y)

  validate_scree_inputs(
    X = Y,
    k = k,
    eps = eps,
    delta = delta
  )

  n_pairs <- n %/% 2L

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

  eps_scree <- eps
  delta_scree <- delta

  # for 2k quantile estimates
  eps_Q <- eps_scree / (4 * k)

  # for mean estimate
  eps_M <- eps_scree / 2
  delta_M <- delta_scree

  W <- paired_score_squares(Y)

  if (isTRUE(split_mode)) {
    if (n_pairs < 2L) {
      stop("split_mode = TRUE requires at least 2 pairs (n >= 4).")
    }

    n_q <- n_pairs %/% 2L
    idx <- sample.int(n_pairs)

    idx_q <- idx[seq_len(n_q)]
    idx_m <- idx[-seq_len(n_q)]
  } else {
    idx_q <- seq_len(n_pairs)
    idx_m <- seq_len(n_pairs)
  }

  n_q <- length(idx_q)
  n_m <- length(idx_m)

  trim_param <- min(max(trim_const / n_q, eta), 0.49)

  scree_base <- numeric(k)
  cutoff_width <- numeric(k)

  for (ell in seq_len(k)) {
    w <- W[, ell]

    L <- unbounded_quantile_upper(
      x = w[idx_q],
      q = trim_param,
      epsilon = eps_Q,
      beta = beta
    )

    U <- unbounded_quantile_upper(
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

    scree_base[ell] <- mean(w_win)
    cutoff_width[ell] <- U - L
  }

  gamma_M <- mu_from_eps_delta(eps_M, delta_M)
  sd_noise <- sqrt(k) * cutoff_width / (n_m * gamma_M)

  scree <- pmax(
    scree_base + stats::rnorm(k, mean = 0, sd = sd_noise),
    0
  )

  if (isTRUE(mono)) {
    scree <- scree_post_processing(scree)
  }

  list(
    scree = scree,
    pve = scree_to_pve(scree)
  )
}
