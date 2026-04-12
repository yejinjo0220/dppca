# ============================================================
# helper_pca.R
# Internal helper functions for PCA-related routines
# ============================================================

#' Preprocess matrix for PCA
#'
#' @description
#' Internal helper that preprocesses a numeric matrix before PCA by applying 
#' optional centering and standardization.
#'
#' @param X Numeric matrix-like object.
#' @param center Logical; whether to center columns.
#' @param standardize Logical; whether to scale columns by their standard deviations.
#'
#' @return A numeric matrix after preprocessing.
#'
#' @keywords internal
prep_matrix_for_pca <- function(X, center = TRUE, standardize = FALSE) {
  X <- as.matrix(X)
  storage.mode(X) <- "double"

  if (!is.numeric(X)) {
    stop("X must be numeric or coercible to a numeric matrix.")
  }
  if (nrow(X) < 2) {
    stop("Need nrow(X) >= 2.")
  }
  if (ncol(X) < 1) {
    stop("Need ncol(X) >= 1.")
  }

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

#' Re-orthonormalize principal component directions
#'
#' @description
#' Internal helper that re-orthonormalizes candidate principal component
#' directions using a QR decomposition and keeps the first \code{k} columns.
#'
#' @param V Numeric matrix.
#' @param k Number of columns to retain.
#'
#' @return A direction matrix with orthonormal columns.
#'
#' @keywords internal
orthonormalize_dir <- function(V, k) {
  V <- as.matrix(V)
  k <- as.integer(k)

  if (!is.numeric(V)) {
    stop("V must be numeric.")
  }
  if (ncol(V) < k || k < 1) {
    stop("k must satisfy 1 <= k <= ncol(V).")
  }

  Q <- qr.Q(qr(V))
  Q[, 1:k, drop = FALSE]
}

#' Compute principal component directions
#'
#' @description
#' Internal helper that returns the principal component directions actually used
#' in downstream routines. By default, it returns the usual non-private
#' principal component directions. When \code{g_dppca = TRUE}, it returns
#' differentially private principal component directions based on the spherical
#' Kendall mechanism.
#'
#' @param X_proc Preprocessed numeric matrix.
#' @param k Number of leading components.
#' @param g_dppca Logical; whether to privatize the principal component
#'   directions.
#' @param eps_dir Privacy epsilon for releasing principal component directions.
#' @param delta_dir Privacy delta for releasing principal component directions.
#' @param cpp.option Logical passed to \code{mech_tau_sph()}.
#'
#' @return A direction matrix whose columns are the principal component directions used.
#'
#' @keywords internal
compute_pc_dir <- function(X_proc, k,
                           g_dppca = FALSE,
                           eps_dir = NULL,
                           delta_dir = NULL,
                           cpp.option = FALSE) {
  X_proc <- as.matrix(X_proc)
  k <- as.integer(k)

  if (!is.numeric(X_proc)) {
    stop("X_proc must be numeric.")
  }
  if (nrow(X_proc) < 2) {
    stop("Need nrow(X_proc) >= 2.")
  }
  if (ncol(X_proc) < 1) {
    stop("Need ncol(X_proc) >= 1.")
  }
  if (k < 1 || k > ncol(X_proc)) {
    stop("k must satisfy 1 <= k <= ncol(X_proc).")
  }
  if (!requireNamespace("rARPACK", quietly = TRUE)) {
    stop("Package 'rARPACK' is required.")
  }

  n <- nrow(X_proc)
  d <- ncol(X_proc)

  # ------------------------------------------------------------
  # Non-private principal component directions
  # ------------------------------------------------------------
  S_np <- stats::cov(X_proc)
  V_np <- rARPACK::eigs_sym(S_np, k = k)$vectors
  V_np <- as.matrix(V_np)

  if (nrow(V_np) != d || ncol(V_np) != k) {
    stop("V_np returned by eigs_sym has wrong dimension.")
  }

  V_used <- orthonormalize_dir(V_np, k = k)

  # ------------------------------------------------------------
  # DP principal component directions
  # ------------------------------------------------------------
  if (isTRUE(g_dppca)) {
    if (is.null(eps_dir) || is.null(delta_dir)) {
      stop("eps_dir and delta_dir must be supplied when g_dppca = TRUE.")
    }
    if (!is.numeric(eps_dir) || length(eps_dir) != 1 ||
        !is.finite(eps_dir) || eps_dir <= 0) {
      stop("eps_dir must be a single positive number.")
    }
    if (!is.numeric(delta_dir) || length(delta_dir) != 1 ||
        !is.finite(delta_dir) || delta_dir <= 0 || delta_dir >= 1) {
      stop("delta_dir must be a single number in (0, 1).")
    }

    sigma_sph <- 2 * sqrt(2 * log(1.25 / delta_dir)) / (n * eps_dir)
    tilde_K <- mech_tau_sph(X_proc, sig = sigma_sph, cpp.option = cpp.option)

    V_dp <- rARPACK::eigs_sym(tilde_K, k = k)$vectors
    V_dp <- as.matrix(V_dp)

    if (nrow(V_dp) != d || ncol(V_dp) != k) {
      stop("V_dp returned by eigs_sym has wrong dimension.")
    }

    V_used <- orthonormalize_dir(V_dp, k = k)
  }

  V_used
}

#' Euclidean norm
#'
#' @description
#' Internal helper that computes the Euclidean norm of a numeric vector.
#'
#' @param x A numeric vector.
#'
#' @return A nonnegative numeric scalar.
#'
#' @keywords internal
norm2 <- function(x) {
  x <- as.numeric(x)
  sqrt(sum(x^2))
}

#' Normalize a vector to unit norm
#'
#' @description
#' Internal helper that normalizes a vector to have Euclidean norm 1.
#'
#' @param x A numeric vector.
#'
#' @return A numeric vector with Euclidean norm 1.
#'
#' @keywords internal
normalize <- function(x) {
  x <- as.numeric(x)
  nx <- norm2(x)

  if (!is.finite(nx) || nx <= 0) {
    stop("Cannot normalize a vector with nonpositive norm.")
  }

  x / nx
}

#' Convert a packed symmetric vector into a matrix
#'
#' @description
#' Internal helper that converts a packed vector representation of a symmetric
#' matrix into the full matrix form.
#'
#' @param v Numeric vector of length \eqn{p(p+1)/2}.
#' @param p Matrix dimension.
#'
#' @return A symmetric \eqn{p \times p} matrix.
#'
#' @keywords internal
vec2mat <- function(v, p) {
  v <- as.numeric(v)
  p <- as.integer(p)

  if (length(p) != 1 || !is.finite(p) || p < 1) {
    stop("p must be a positive integer.")
  }

  q <- length(v)
  if (q != p * (p + 1) / 2) {
    stop("length(v) is not compatible with p.")
  }

  out <- matrix(0, p, p)
  diag(out) <- v[1:p]

  if (q > p) {
    out[lower.tri(out)] <- v[(p + 1):q] / sqrt(2)
    out[upper.tri(out)] <- t(out)[upper.tri(out)]
  }

  out
}

#' Spherical Kendall matrix
#'
#' @description
#' Internal helper that computes the spherical Kendall matrix from pairwise
#' normalized differences of observations.
#'
#' @param X Numeric matrix.
#' @param cpp.option Logical; currently only \code{FALSE} is supported.
#'
#' @return A symmetric matrix.
#'
#' @keywords internal
tau_sph <- function(X, cpp.option = FALSE) {
  X <- as.matrix(X)
  n <- nrow(X)
  d <- ncol(X)

  if (!is.numeric(X)) stop("X must be numeric.")
  if (n < 2) stop("Need n >= 2.")
  if (d < 1) stop("Need ncol(X) >= 1.")

  if (isTRUE(cpp.option)) {
    stop(
      "cpp.option = TRUE is not yet supported in this package. ",
      "Please use cpp.option = FALSE."
    )
  }

  hK <- matrix(0, d, d)

  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      a <- X[j, ] - X[i, ]
      a <- normalize(a)
      hK <- hK + a %*% t(a)
    }
  }

  hK * (2 / (n * (n - 1)))
}

#' Gaussian mechanism applied to the spherical Kendall matrix
#'
#' @description
#' Internal helper that adds Gaussian noise to the spherical Kendall matrix to
#' construct a noisy symmetric matrix for DP principal component direction
#' estimation.
#'
#' @param X Numeric matrix.
#' @param sig Noise standard deviation.
#' @param cpp.option Logical; currently only \code{FALSE} is supported.
#'
#' @return A noisy symmetric matrix.
#'
#' @keywords internal
mech_tau_sph <- function(X, sig, cpp.option = FALSE) {
  X <- as.matrix(X)
  n <- nrow(X)
  d <- ncol(X)

  if (!is.numeric(X)) stop("X must be numeric.")
  if (n < 2) stop("Need n >= 2.")
  if (d < 1) stop("Need ncol(X) >= 1.")
  if (!is.numeric(sig) || length(sig) != 1 || !is.finite(sig) || sig < 0) {
    stop("sig must be a single nonnegative number.")
  }

  q <- d * (d + 1) / 2
  hK <- tau_sph(X, cpp.option = cpp.option)
  zeta <- stats::rnorm(q, mean = 0, sd = sig)
  zeta_mat <- vec2mat(zeta, d)

  hK + zeta_mat
}
