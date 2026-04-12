#' Differentially private eigenvalue estimation via Huber-type robust mean
#'
#' This function estimates leading PCA eigenvalues by first obtaining a PCA
#' basis, then applying a scalar robust private mean estimator to the projected
#' squared scores along each component.
#'
#' @param X Numeric data matrix with observations in rows.
#' @param k Number of leading components.
#' @param eps_total Total privacy budget epsilon.
#' @param delta_total Total privacy budget delta.
#' @param dp_v_flag Logical; whether to privatize the PCA basis.
#' @param cpp.option Logical passed to \code{mech_tau_sph()}.
#' @param center Logical; whether to center columns before PCA.
#' @param standardize Logical; whether to standardize columns before PCA.
#' @param mu0 Initial point for the scalar noisy GD updates.
#' @param eta0 Base step size for noisy GD.
#' @param T Number of GD iterations. If \code{NULL}, defaults to \code{ceiling(log(n))}.
#' @param M Number of blocks for the scale proxy routine. If \code{NULL},
#'   defaults to \code{floor(sqrt(n/2))}.
#' @param k_min_star,k_max_star Integer range of dyadic histogram bins.
#' @param include_zero_bin Logical; whether to include a separate zero bin.
#' @param m2_frac Fraction of the eigenvalue budget allocated to the scale proxy.
#' @param step_schedule Either \code{"fixed"} or \code{"1/sqrt(t)"}.
#' @param mono Logical; whether to apply monotone post-processing.
#'
#' @return A list containing private eigenvalue estimates and auxiliary objects.
#' @export
dp_lambda_huber <- function(
    X, k,
    eps_total, delta_total,
    dp_v_flag = FALSE,
    cpp.option = FALSE,
    center = TRUE,
    standardize = TRUE,
    mu0 = 0,
    eta0 = 1,
    T = NULL,
    M = NULL,
    k_min_star = -40L,
    k_max_star = 40L,
    include_zero_bin = TRUE,
    m2_frac = 1 / 4,
    step_schedule = c("fixed", "1/sqrt(t)"),
    mono = TRUE
) {
  validate_dp_pca_inputs(
    X = X,
    k = k,
    eps_total = eps_total,
    delta_total = delta_total
  )

  step_schedule <- match.arg(step_schedule)

  X <- as.matrix(X)
  n <- nrow(X)
  d <- ncol(X)
  k <- as.integer(k)

  if (!is.logical(dp_v_flag) || length(dp_v_flag) != 1) {
    stop("dp_v_flag must be TRUE/FALSE.")
  }
  if (!is.logical(mono) || length(mono) != 1) {
    stop("mono must be TRUE/FALSE.")
  }
  if (!is.finite(eta0) || eta0 <= 0) stop("eta0 must be > 0.")
  if (!is.finite(m2_frac) || m2_frac <= 0 || m2_frac >= 1) {
    stop("m2_frac must be in (0,1).")
  }

  k_min_star <- as.integer(k_min_star)
  k_max_star <- as.integer(k_max_star)
  if (k_min_star > k_max_star) stop("Need k_min_star <= k_max_star.")

  if (is.null(T)) T <- ceiling(log(n))
  T <- as.integer(T)
  if (T < 1) T <- 1L

  if (is.null(M)) M <- floor(sqrt(n / 2))
  M <- as.integer(M)
  if (M < 1) M <- 1L

  # privacy split: basis vs eigenvalues
  if (isTRUE(dp_v_flag)) {
    eps_V <- eps_total / 2
    del_V <- delta_total / 2
    eps_L <- eps_total / 2
    del_L <- delta_total / 2
  } else {
    eps_V <- 0
    del_V <- 0
    eps_L <- eps_total
    del_L <- delta_total
  }

  # within lambda budget: m2 vs GD
  eps_m2 <- eps_L * m2_frac
  del_m2 <- del_L * m2_frac
  eps_gd <- eps_L * (1 - m2_frac)
  del_gd <- del_L * (1 - m2_frac)

  eps_m2_ell <- eps_m2 / k
  del_m2_ell <- del_m2 / k
  eps_gd_ell <- eps_gd / k
  del_gd_ell <- del_gd / k

  X_proc <- prep_matrix_for_pca(
    X = X,
    center = center,
    standardize = standardize
  )

  basis <- compute_pca_bases(
    X_proc = X_proc,
    k = k,
    dp_v_flag = dp_v_flag,
    eps_vec = eps_V,
    delta_vec = del_V,
    cpp.option = cpp.option
  )

  V_np <- basis$V_np
  V_used <- basis$V_used

  Y <- X_proc %*% V_used

  lambda_dp_raw <- numeric(k)
  m2_hat_vec <- numeric(k)
  tau_vec <- numeric(k)
  lambda_hat <- numeric(k)

  for (ell in seq_len(k)) {
    y <- Y[, ell]
    ybar <- mean(y)
    w <- (y - ybar)^2

    alg2 <- dp_m2_alg2_scalar(
      w = w,
      eps_m2 = eps_m2_ell,
      M = M,
      k_min_star = k_min_star,
      k_max_star = k_max_star,
      include_zero_bin = include_zero_bin
    )

    m2_hat <- alg2$m2_hat
    m2_hat_vec[ell] <- m2_hat

    tau_ell <- tau_from_m2_scalar(
      m2_hat = m2_hat,
      eps_gd_ell = eps_gd_ell,
      n_eff = n
    )
    if (!is.finite(tau_ell) || tau_ell <= 0) tau_ell <- 1
    tau_vec[ell] <- tau_ell

    muT <- dp_huber_noisy_gd_scalar(
      w = w,
      eps_gd = eps_gd_ell,
      delta_gd = del_gd_ell,
      tau = tau_ell,
      T = T,
      mu0 = mu0,
      eta0 = eta0,
      step_schedule = step_schedule
    )

    lambda_hat[ell] <- mean(w)
    lambda_dp_raw[ell] <- max(muT, 0)
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
    m2_hat = m2_hat_vec,
    tau = tau_vec,
    settings = list(
      n = n,
      d = d,
      k = k,
      eps_total = eps_total,
      delta_total = delta_total,
      dp_v_flag = dp_v_flag,
      eps_V = eps_V,
      del_V = del_V,
      eps_lambda = eps_L,
      del_lambda = del_L,
      m2_frac = m2_frac,
      eps_m2_ell = eps_m2_ell,
      del_m2_ell = del_m2_ell,
      eps_gd_ell = eps_gd_ell,
      del_gd_ell = del_gd_ell,
      mu0 = mu0,
      eta0 = eta0,
      T = T,
      M = M,
      k_min_star = k_min_star,
      k_max_star = k_max_star,
      include_zero_bin = include_zero_bin,
      step_schedule = step_schedule,
      mono = mono
    )
  )
}
