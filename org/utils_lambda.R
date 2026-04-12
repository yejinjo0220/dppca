#' Validate common PCA inputs
#'
#' @param X Numeric matrix-like object.
#' @param k Number of components.
#' @param eps_total Privacy epsilon.
#' @param delta_total Privacy delta.
#'
#' @return Invisibly returns TRUE.
#' @keywords internal
validate_dp_pca_inputs <- function(X, k, eps_total, delta_total) {
  X <- as.matrix(X)
  n <- nrow(X)
  d <- ncol(X)

  if (n < 2) stop("Need n >= 2.")
  if (d < 1) stop("Need ncol(X) >= 1.")

  if (!is.numeric(k) || length(k) != 1 || !is.finite(k) || k < 1 || k > d) {
    stop("k must be an integer in {1, ..., ncol(X)}.")
  }

  if (!is.numeric(eps_total) || length(eps_total) != 1 || !is.finite(eps_total) || eps_total <= 0) {
    stop("eps_total must be a single positive number.")
  }

  if (!is.numeric(delta_total) || length(delta_total) != 1 ||
      !is.finite(delta_total) || delta_total <= 0 || delta_total >= 1) {
    stop("delta_total must be a single number in (0, 1).")
  }

  invisible(TRUE)
}

#' Preprocess matrix for PCA
#'
#' @param X Numeric matrix-like object.
#' @param center Logical; whether to center columns.
#' @param standardize Logical; whether to scale columns by their SD.
#'
#' @return A numeric matrix.
#' @keywords internal
prep_matrix_for_pca <- function(X, center = TRUE, standardize = TRUE) {
  X <- as.matrix(X)
  storage.mode(X) <- "double"

  X_proc <- X

  if (isTRUE(center)) {
    X_proc <- scale(X_proc, center = TRUE, scale = FALSE)
  }

  if (isTRUE(standardize)) {
    sds <- apply(X_proc, 2, stats::sd)
    if (any(!is.finite(sds)) || any(sds == 0)) {
      stop("Some columns have sd = 0 or non-finite sd.")
    }
    X_proc <- sweep(X_proc, 2, sds, "/")
  }

  X_proc
}

#' Re-orthonormalize basis vectors
#'
#' @param V Numeric matrix.
#' @param k Number of columns to retain.
#'
#' @return A d x k orthonormal matrix.
#' @keywords internal
orthonormalize_basis <- function(V, k) {
  V <- as.matrix(V)
  k <- as.integer(k)

  Q <- qr.Q(qr(V))
  Q[, 1:k, drop = FALSE]
}

#' Compute non-private PCA basis
#'
#' @param X_proc Preprocessed numeric matrix.
#' @param k Number of leading components.
#'
#' @return A d x k orthonormal matrix.
#' @keywords internal
compute_np_pca_basis <- function(X_proc, k) {
  if (!requireNamespace("rARPACK", quietly = TRUE)) {
    stop("Package 'rARPACK' is required.")
  }

  X_proc <- as.matrix(X_proc)
  d <- ncol(X_proc)

  S_np <- stats::cov(X_proc)
  V_np <- rARPACK::eigs_sym(S_np, k = k)$vectors
  V_np <- as.matrix(V_np)

  if (nrow(V_np) != d || ncol(V_np) != k) {
    stop("V_np returned by eigs_sym has wrong dimension.")
  }

  orthonormalize_basis(V_np, k = k)
}

#' Compute optionally private PCA basis
#'
#' @param X_proc Preprocessed numeric matrix.
#' @param k Number of leading components.
#' @param dp_v_flag Logical; whether to privatize the basis.
#' @param eps_vec Privacy epsilon for basis release.
#' @param delta_vec Privacy delta for basis release.
#' @param cpp.option Logical passed to \code{mech_tau_sph()}.
#'
#' @return A list with \code{V_np} and \code{V_used}.
#' @keywords internal
compute_pca_bases <- function(X_proc, k,
                              dp_v_flag = FALSE,
                              eps_vec = 0,
                              delta_vec = 0,
                              cpp.option = FALSE) {
  X_proc <- as.matrix(X_proc)
  n <- nrow(X_proc)
  d <- ncol(X_proc)

  V_np <- compute_np_pca_basis(X_proc, k = k)
  V_used <- V_np

  if (isTRUE(dp_v_flag)) {
    if (!requireNamespace("rARPACK", quietly = TRUE)) {
      stop("Package 'rARPACK' is required.")
    }

    sigma_sph <- 2 * sqrt(2 * log(1.25 / delta_vec)) / (n * eps_vec)
    tilde_K <- mech_tau_sph(X_proc, sig = sigma_sph, cpp.option = cpp.option)
    V_dp <- rARPACK::eigs_sym(tilde_K, k = k)$vectors
    V_dp <- as.matrix(V_dp)

    if (nrow(V_dp) != d || ncol(V_dp) != k) {
      stop("V_dp returned by eigs_sym has wrong dimension.")
    }

    V_used <- orthonormalize_basis(V_dp, k = k)
  }

  list(
    V_np = V_np,
    V_used = V_used
  )
}

#' Enforce nonnegative nonincreasing sequence
#'
#' @param x Numeric vector.
#'
#' @return A numeric vector of the same length.
#' @keywords internal
lambda_mono_decreasing_nonneg <- function(x) {
  y <- pmax(as.numeric(x), 0)
  fit <- stats::isoreg(seq_along(y), -y)
  pmax(-fit$yf, 0)
}

#' Compute explained variance ratios safely
#'
#' @param lambda Numeric vector of eigenvalue estimates.
#'
#' @return Numeric vector of the same length.
#' @keywords internal
safe_evr <- function(lambda) {
  lambda <- as.numeric(lambda)
  s <- sum(lambda)

  if (!is.finite(s) || s <= 0) {
    return(rep(0, length(lambda)))
  }

  lambda / s
}
