#' Differentially private eigenvalue estimation via PMWM
#'
#' This function estimates leading PCA eigenvalues by first obtaining a PCA
#' basis, then computing projected squared scores, privately estimating lower
#' and upper trimming bounds via an unbounded quantile mechanism, winsorizing
#' the projected squared scores, and finally privatizing the winsorized mean
#' with Gaussian noise.
#'
#' @param X Numeric data matrix with observations in rows.
#' @param k Number of leading components.
#' @param eps_total Total privacy budget epsilon.
#' @param delta_total Total privacy budget delta.
#' @param dp_v_flag Logical; whether to privatize the PCA basis.
#' @param cpp.option Logical passed to \code{mech_tau_sph()}.
#' @param split_mode Logical; whether to split the sample into quantile and mean
#'   subsets.
#' @param center Logical; whether to center columns before PCA.
#' @param standardize Logical; whether to standardize columns before PCA.
#' @param beta Log-binning base; must be greater than 1.
#' @param a Lower support bound for the quantile estimator.
#' @param b Upper support bound for the quantile estimator.
#' @param trim_const Positive constant controlling the trimming proportion.
#' @param eta Lower bound for the trimming proportion.
#' @param mono Logical; whether to apply monotone post-processing.
#' @param max_extra_bins Integer buffer used by the unbounded quantile scan.
#'
#' @return A list containing private eigenvalue estimates and auxiliary objects.
#' @export
dp_lambda_pmwm <- function(
    X, k,
    eps_total, delta_total,
    dp_v_flag = FALSE,
    cpp.option = FALSE,
    split_mode = TRUE,
    center = TRUE,
    standardize = TRUE,
    beta = 1.01,
    a, b,
    trim_const = 10,
    eta = 0.01,
    mono = TRUE,
    max_extra_bins = 1000
) {
  validate_dp_pca_inputs(
    X = X,
    k = k,
    eps_total = eps_total,
    delta_total = delta_total
  )

  X <- as.matrix(X)
  n <- nrow(X)
  d <- ncol(X)
  k <- as.integer(k)

  if (!is.logical(dp_v_flag) || length(dp_v_flag) != 1) {
    stop("dp_v_flag must be TRUE/FALSE.")
  }
  if (!is.logical(split_mode) || length(split_mode) != 1) {
    stop("split_mode must be TRUE/FALSE.")
  }
  if (!is.logical(mono) || length(mono) != 1) {
    stop("mono must be TRUE/FALSE.")
  }

  if (!is.finite(beta) || beta <= 1) stop("beta must be > 1.")
  if (!is.finite(a) || !is.finite(b) || a > b) stop("Need finite a <= b.")
  if (!is.numeric(trim_const) || length(trim_const) != 1 || trim_const <= 0) {
    stop("trim_const must be a single positive number.")
  }
  if (!is.numeric(eta) || length(eta) != 1 || eta <= 0 || eta >= 0.5) {
    stop("eta must be in (0, 0.5).")
  }

  max_extra_bins <- as.integer(max_extra_bins)
  if (!is.finite(max_extra_bins) || max_extra_bins < 0) {
    stop("max_extra_bins must be a nonnegative integer.")
  }

  # privacy split: basis vs eigenvalues
  if (isTRUE(dp_v_flag)) {
    eps_vec   <- eps_total / 2
    del_vec   <- delta_total / 2
    eps_lam_T <- eps_total / 2
    del_lam_T <- delta_total / 2
  } else {
    eps_vec   <- 0
    del_vec   <- 0
    eps_lam_T <- eps_total
    del_lam_T <- delta_total
  }

  # per-component budget
  eps_ell <- eps_lam_T / k
  del_ell <- del_lam_T / k

  # split component budget: quantiles + mean
  eps_Q <- eps_ell / 4
  del_Q <- del_ell / 4
  eps_M <- eps_ell / 2
  del_M <- del_ell / 2

  eps_q1 <- eps_Q / 2
  del_q1 <- del_Q / 2
  eps_q2 <- eps_Q / 2
  del_q2 <- del_Q / 2

  # preprocess + basis
  X_proc <- prep_matrix_for_pca(
    X = X,
    center = center,
    standardize = standardize
  )

  basis <- compute_pca_bases(
    X_proc = X_proc,
    k = k,
    dp_v_flag = dp_v_flag,
    eps_vec = eps_vec,
    delta_vec = del_vec,
    cpp.option = cpp.option
  )

  V_np <- basis$V_np
  V_used <- basis$V_used

  # projected scores
  Y <- X_proc %*% V_used

  # split for quantiles / mean
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
  n_m <- length(idx_m)

  # trimming ratio
  cap <- max(0.025, eta)
  trim_param <- min(max(trim_const / n, eta), cap)

  lambda_hat <- numeric(k)
  lambda_dp_raw <- numeric(k)
  L_vec <- numeric(k)
  U_vec <- numeric(k)

  for (ell in seq_len(k)) {
    y <- Y[, ell]
    ybar <- mean(y)
    w <- (y - ybar)^2

    L <- unbounded_quantile(
      data = w[idx_q],
      l = a, u = b, b = beta,
      q = trim_param,
      eps_1 = eps_q1, delta_1 = del_q1,
      eps_2 = eps_q2, delta_2 = del_q2,
      max_extra_bins = max_extra_bins
    )

    U <- unbounded_quantile(
      data = w[idx_q],
      l = a, u = b, b = beta,
      q = 1 - trim_param,
      eps_1 = eps_q1, delta_1 = del_q1,
      eps_2 = eps_q2, delta_2 = del_q2,
      max_extra_bins = max_extra_bins
    )

    if (!is.finite(L) || !is.finite(U) || U < L) {
      L <- a
      U <- b
    }

    w_clip <- pmin(pmax(w[idx_m], L), U)
    mu_hat <- mean(w_clip)

    Delta <- (U - L) / n_m
    sd_mean <- Delta * sqrt(2 * log(1.25 / del_M)) / eps_M
    mu_dp <- mu_hat + stats::rnorm(1, mean = 0, sd = sd_mean)

    L_vec[ell] <- L
    U_vec[ell] <- U
    lambda_hat[ell] <- mu_hat
    lambda_dp_raw[ell] <- max(mu_dp, 0)
  }

  if (isTRUE(mono)) {
    lambda_dp <- lambda_mono_decreasing_nonneg(lambda_dp_raw)
  } else {
    lambda_dp <- lambda_dp_raw
  }

  evr <- safe_evr(lambda_dp)

  list(
    lambda_dp = lambda_dp,
    lambda_dp_raw = lambda_dp_raw,
    lambda_hat = lambda_hat,
    evr = evr,
    V_np = V_np,
    V_used = V_used,
    bounds = data.frame(
      comp = seq_len(k),
      L = L_vec,
      U = U_vec
    ),
    settings = list(
      n = n,
      d = d,
      k = k,
      eps_total = eps_total,
      delta_total = delta_total,
      dp_v_flag = dp_v_flag,
      eps_vec = eps_vec,
      delta_vec = del_vec,
      eps_lam_total = eps_lam_T,
      delta_lam_total = del_lam_T,
      eps_component = eps_ell,
      delta_component = del_ell,
      split_mode = split_mode,
      center = center,
      standardize = standardize,
      beta = beta,
      a = a,
      b = b,
      trim_const = trim_const,
      eta = eta,
      trim_param = trim_param,
      mono = mono,
      max_extra_bins = max_extra_bins
    )
  )
}
