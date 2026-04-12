#' Compute differentially private scree estimates
#'
#' User-facing computation function for DP PCA eigenvalue estimation.
#' This function returns both the selected DP scree estimate and the
#' corresponding non-private PCA reference values.
#'
#' @param X Numeric data matrix with observations in rows.
#' @param k Number of leading components.
#' @param method One of \code{"clipped"}, \code{"pmwm"}, or \code{"huber"}.
#' @param eps_total Total privacy budget epsilon.
#' @param delta_total Total privacy budget delta.
#' @param center Logical; whether to center columns before PCA.
#' @param standardize Logical; whether to standardize columns before PCA.
#' @param dp_v_flag Logical; whether to privatize the PCA basis.
#' @param cpp.option Logical passed to \code{mech_tau_sph()}.
#' @param mono Logical; whether to apply monotone postprocessing.
#' @param B_mult Clipping multiplier for the clipped estimator.
#' @param beta Log-binning base for PMWM.
#' @param a,b Support bounds for PMWM.
#' @param trim_const,eta PMWM trimming parameters.
#' @param split_mode Logical; whether PMWM splits the sample.
#' @param mu0,eta0,T,M Huber estimator parameters.
#' @param k_min_star,k_max_star,include_zero_bin Huber scale-estimation parameters.
#' @param m2_frac Huber privacy split parameter.
#' @param step_schedule Huber GD step schedule.
#'
#' @return A list containing:
#' \describe{
#'   \item{method}{Selected DP method.}
#'   \item{lambda_np}{Non-private PCA eigenvalues.}
#'   \item{evr_np}{Non-private explained variance ratios.}
#'   \item{lambda_dp}{DP eigenvalue estimates.}
#'   \item{evr}{DP explained variance ratios.}
#'   \item{result}{Full method-specific output from the selected estimator.}
#'   \item{settings}{Common input settings.}
#' }
#' @export
dp_scree <- function(
    X,
    k,
    method = c("clipped", "pmwm", "huber"),
    eps_total,
    delta_total,
    center = TRUE,
    standardize = TRUE,
    dp_v_flag = FALSE,
    cpp.option = FALSE,
    mono = TRUE,
    B_mult = 2,
    beta = 1.01,
    a = NULL,
    b = NULL,
    trim_const = 10,
    eta = 0.01,
    split_mode = TRUE,
    mu0 = 0,
    eta0 = 1,
    T = NULL,
    M = NULL,
    k_min_star = -40L,
    k_max_star = 40L,
    include_zero_bin = TRUE,
    m2_frac = 1 / 4,
    step_schedule = c("fixed", "1/sqrt(t)")
) {
  method <- match.arg(method)
  step_schedule <- match.arg(step_schedule)

  X <- as.matrix(X)
  validate_dp_pca_inputs(X, k, eps_total, delta_total)

  if (!requireNamespace("rARPACK", quietly = TRUE)) {
    stop("Package 'rARPACK' is required.")
  }

  if (method == "pmwm" && (is.null(a) || is.null(b))) {
    stop("For method = 'pmwm', you must provide 'a' and 'b'.")
  }

  X_proc <- prep_matrix_for_pca(
    X,
    center = center,
    standardize = standardize
  )

  S_np <- stats::cov(X_proc)
  eig_np <- rARPACK::eigs_sym(S_np, k = k)
  lambda_np <- as.numeric(eig_np$values)
  evr_np <- safe_evr(lambda_np)

  result <- switch(
    method,
    clipped = dp_lambda_clipped(
      X = X,
      k = k,
      eps_total = eps_total,
      delta_total = delta_total,
      center = center,
      standardize = standardize,
      B_mult = B_mult,
      dp_v_flag = dp_v_flag,
      cpp.option = cpp.option,
      mono = mono
    ),
    pmwm = dp_lambda_pmwm(
      X = X,
      k = k,
      eps_total = eps_total,
      delta_total = delta_total,
      dp_v_flag = dp_v_flag,
      cpp.option = cpp.option,
      split_mode = split_mode,
      center = center,
      standardize = standardize,
      beta = beta,
      a = a,
      b = b,
      trim_const = trim_const,
      eta = eta,
      mono = mono
    ),
    huber = dp_lambda_huber(
      X = X,
      k = k,
      eps_total = eps_total,
      delta_total = delta_total,
      dp_v_flag = dp_v_flag,
      cpp.option = cpp.option,
      center = center,
      standardize = standardize,
      mu0 = mu0,
      eta0 = eta0,
      T = T,
      M = M,
      k_min_star = k_min_star,
      k_max_star = k_max_star,
      include_zero_bin = include_zero_bin,
      m2_frac = m2_frac,
      step_schedule = step_schedule,
      mono = mono
    )
  )

  list(
    method = method,
    lambda_np = lambda_np,
    evr_np = evr_np,
    lambda_dp = result$lambda_dp,
    evr = result$evr,
    result = result,
    settings = list(
      n = nrow(X),
      d = ncol(X),
      k = k,
      eps_total = eps_total,
      delta_total = delta_total,
      center = center,
      standardize = standardize,
      dp_v_flag = dp_v_flag,
      mono = mono
    )
  )
}
