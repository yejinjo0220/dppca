#' Differentially private eigenvalue estimation via clipped projected variances
#'
#' This function estimates the leading PCA eigenvalues using projected score
#' variances after rowwise L2 clipping, followed by Gaussian perturbation.
#'
#' @param X Numeric data matrix with observations in rows.
#' @param k Number of leading components.
#' @param eps_total Total privacy budget epsilon.
#' @param delta_total Total privacy budget delta.
#' @param center Logical; whether to center columns before PCA.
#' @param standardize Logical; whether to standardize columns before PCA.
#' @param B_mult Positive clipping multiplier; clipping radius is
#'   \code{B_mult * sqrt(d)}.
#' @param dp_v_flag Logical; whether to privatize the PCA basis.
#' @param cpp.option Logical passed to \code{mech_tau_sph()}.
#' @param mono Logical; whether to apply monotone post-processing.
#'
#' @return A list containing private eigenvalue estimates and auxiliary objects.
#' @export
dp_lambda_clipped <- function(
    X,
    k,
    eps_total,
    delta_total,
    center = TRUE,
    standardize = TRUE,
    B_mult = 2,
    dp_v_flag = FALSE,
    cpp.option = FALSE,
    mono = TRUE
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

  if (!is.numeric(B_mult) || length(B_mult) != 1 || !is.finite(B_mult) || B_mult <= 0) {
    stop("B_mult must be a single positive number.")
  }
  if (!is.logical(dp_v_flag) || length(dp_v_flag) != 1) {
    stop("dp_v_flag must be TRUE/FALSE.")
  }
  if (!is.logical(mono) || length(mono) != 1) {
    stop("mono must be TRUE/FALSE.")
  }

  # privacy split: basis vs eigenvalues
  if (isTRUE(dp_v_flag)) {
    eps_vec   <- eps_total / 2
    delta_vec <- delta_total / 2
    eps_lam   <- eps_total / 2
    delta_lam <- delta_total / 2
  } else {
    eps_vec   <- 0
    delta_vec <- 0
    eps_lam   <- eps_total
    delta_lam <- delta_total
  }

  # preprocess
  X_proc <- prep_matrix_for_pca(
    X = X,
    center = center,
    standardize = standardize
  )

  # PCA bases
  basis <- compute_pca_bases(
    X_proc = X_proc,
    k = k,
    dp_v_flag = dp_v_flag,
    eps_vec = eps_vec,
    delta_vec = delta_vec,
    cpp.option = cpp.option
  )

  V_np <- basis$V_np
  V_used <- basis$V_used

  # rowwise L2 clipping
  B <- B_mult * sqrt(d)
  row_norm <- sqrt(rowSums(X_proc^2))
  scale_fac <- pmin(1, B / pmax(row_norm, 1e-12))
  X_clip <- X_proc * scale_fac

  # projected variances
  Y <- X_clip %*% V_used
  lambda_hat <- apply(Y, 2, stats::var)

  # Gaussian mechanism for eigenvalues
  Delta1 <- (4 * B^2) / n
  eps_each <- eps_lam / k
  delta_each <- delta_lam / k
  sd_noise <- Delta1 * sqrt(2 * log(1.25 / delta_each)) / eps_each

  lambda_dp_raw <- lambda_hat + stats::rnorm(k, mean = 0, sd = sd_noise)
  lambda_dp_raw <- pmax(lambda_dp_raw, 0)

  # monotone post-processing
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
    settings = list(
      n = n,
      d = d,
      k = k,
      eps_total = eps_total,
      delta_total = delta_total,
      center = center,
      standardize = standardize,
      B_mult = B_mult,
      B = B,
      dp_v_flag = dp_v_flag,
      eps_vec = eps_vec,
      delta_vec = delta_vec,
      eps_lam = eps_lam,
      delta_lam = delta_lam,
      mono = mono
    )
  )
}
