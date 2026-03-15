#' Differentially Private PCA
#'
#' @description
#' Performs (optionally differentially private) principal component analysis (PCA)
#' on a numeric data matrix. When \code{dp = TRUE}, a tau-spherical DP mechanism
#' is applied to the (scaled) data before computing the eigendecomposition.
#'
#' @param X A numeric matrix or an object coercible to a matrix with
#'   \eqn{n} rows (observations) and \eqn{p} columns (variables).
#' @param center Logical; should the variables be centered before PCA?
#'   Passed to \code{scale()}.
#' @param scale. Logical; should the variables be scaled to unit variance?
#'   Passed to \code{scale()}.
#' @param axes Integer vector indicating which principal components to return
#'   (e.g., \code{c(1, 2)} for the first two PCs). The maximum index in
#'   \code{axes} must not exceed the number of variables \code{p}.
#' @param dp Logical; if \code{TRUE}, perform differentially private PCA using
#'   a tau-spherical mechanism applied to the empirical covariance/kernel.
#' @param eps Privacy budget \eqn{\varepsilon} to be used when \code{dp = TRUE}.
#'   Must be strictly positive.
#' @param delta Privacy parameter \eqn{\delta} to be used when \code{dp = TRUE}.
#'   Must satisfy \eqn{0 < \delta < 1}.
#' @param cpp.option Logical; if \code{TRUE}, use the C++ implementation inside
#'   \code{mech_tau_sph()} when available, otherwise use a pure R version.
#'
#' @details
#' The data matrix \code{X} is first centered and/or scaled using \code{scale()}
#' according to \code{center} and \code{scale.}. For the non-DP case
#' (\code{dp = FALSE}), the empirical Gram matrix \eqn{G = X^\top X} is formed,
#' and the leading eigenvalues/eigenvectors are obtained via
#' \code{rARPACK::eigs_sym()}.
#'
#' For the DP case (\code{dp = TRUE}), a Gaussian noise level
#' \code{sig_sph} is computed from \code{eps} and \code{delta}, and the
#' tau-spherical mechanism \code{mech_tau_sph()} is applied to the processed
#' data \code{X_proc} to obtain a privatized kernel \code{tilde_Ksph}. The
#' leading eigenvalues/eigenvectors of this privatized kernel are then used to
#' construct differentially private scores and loadings.
#'
#' @return A list with the following components:
#' \itemize{
#'   \item \code{X_pca}: An \eqn{n \times |\mathrm{axes}|} matrix of PCA scores
#'         for the requested components (PC axes).
#'   \item \code{loadings}: A \eqn{p \times |\mathrm{axes}|} matrix of loadings
#'         (eigenvectors) corresponding to the selected components.
#'   \item \code{eigvals}: A numeric vector of the selected eigenvalues.
#'   \item \code{dp}: Logical flag indicating whether DP-PCA was used.
#'   \item \code{eps}: Privacy budget \eqn{\varepsilon} (only for \code{dp = TRUE}).
#'   \item \code{delta}: Privacy parameter \eqn{\delta} (only for \code{dp = TRUE}).
#'   \item \code{sigma}: Noise scale used in the tau-spherical mechanism
#'         (only for \code{dp = TRUE}).
#' }
#'
#' @seealso \code{\link{dp_score_plot}}, \code{\link{dp_score_plot_group}}
#'
#' @importFrom rARPACK eigs_sym
#' @export
dp_pca <- function(X, center = TRUE, scale. = FALSE, axes = c(1, 2),
                   dp = FALSE, eps = NULL, delta = NULL,
                   cpp.option = FALSE) {

  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)

  if (length(axes) < 2) {
    stop("axes must contain at least two principal component indices (e.g. c(1, 2)).")
  }
  k_max <- max(axes)
  if (k_max > p) {
    stop("Requested principal component index in axes exceeds the number of variables in X.")
  }

  if (dp) {
    if (is.null(eps) || is.null(delta)) {
      stop("When 'dp = TRUE', both eps and delta must be specified.")
    }
  }

  X_proc <- scale(X, center = center, scale = scale.)

  # Non-DP PCA ---------------------------------------------------------
  if (!dp) {
    G <- t(X_proc) %*% X_proc

    eig_res <- rARPACK::eigs_sym(G, k = k_max)
    U <- eig_res$vectors
    eigvals <- eig_res$values

    loadings <- U[, axes, drop = FALSE]
    scores   <- X_proc %*% loadings

    colnames(scores) <- paste0("PC", axes)
    rownames(scores) <- rownames(X)

    return(list(
      X_pca   = scores,
      loadings = loadings,
      eigvals = eigvals[axes],
      dp      = FALSE
    ))
  }

  # DP PCA -------------------------------------------------------------
  sig_sph <- 4 * sqrt(2 * log(1.25 / (delta / 2))) / (n * (eps / 2))

  tilde_Ksph <- mech_tau_sph(X_proc, sig_sph, cpp.option = cpp.option)

  eig_dp <- rARPACK::eigs_sym(tilde_Ksph, k = k_max)
  V_dp <- eig_dp$vectors
  eigvals_dp <- eig_dp$values

  loadings_dp <- V_dp[, axes, drop = FALSE]
  scores_dp   <- X_proc %*% loadings_dp

  colnames(scores_dp) <- paste0("PC", axes)
  rownames(scores_dp) <- rownames(X)

  return(list(
    X_pca   = scores_dp,
    loadings = loadings_dp,
    eigvals = eigvals_dp[axes],
    dp      = TRUE,
    eps     = eps,
    delta   = delta,
    sigma   = sig_sph
  ))
}
