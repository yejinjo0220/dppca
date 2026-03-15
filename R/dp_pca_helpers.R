# helper functions for DP-PCA ----------------------------------------------------

#' Euclidean norm
#'
#' Compute the Euclidean norm of a numeric vector.
#'
#' @param x A numeric vector.
#'
#' @return A nonnegative numeric scalar.
#' @keywords internal
norm2 <- function(x) {
  x <- as.numeric(x)
  sqrt(sum(x^2))
}

#' Normalize a vector to unit norm
#'
#' @param x A numeric vector.
#'
#' @return A vector with Euclidean norm 1.
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
#' The input vector is assumed to contain the diagonal entries first,
#' followed by the lower-triangular entries scaled by sqrt(2), matching
#' the convention used in the original research code.
#'
#' @param v Numeric vector of length p * (p + 1) / 2.
#' @param p Matrix dimension.
#'
#' @return A symmetric p x p matrix.
#' @keywords internal
vec2mat <- function(v, p) {
  v <- as.numeric(v)
  p <- as.integer(p)

  if (length(p) != 1 || !is.finite(p) || p < 1) {
    stop("p must be a positive integer.")
  }

  q <- length(v)
  if (q != p * (p + 1) / 2) {
    stop("dimension is not compatible")
  }

  out <- matrix(0, p, p)
  diag(out) <- v[1:p]

  if (q > p) {
    out[lower.tri(out)] <- v[(p + 1):q] / sqrt(2)
    out[upper.tri(out)] <- t(out)[upper.tri(out)]
  }

  out
}

#' Sample multivariate Kendall's tau matrix under spherical transformation
#'
#' @param X Numeric matrix.
#' @param cpp.option Logical; currently only FALSE is supported.
#'
#' @return A symmetric d x d matrix.
#' @keywords internal
tau_sph <- function(X, cpp.option = FALSE) {
  X <- as.matrix(X)
  n <- nrow(X)
  d <- ncol(X)

  if (n < 2) stop("Need n >= 2.")
  if (d < 1) stop("Need ncol(X) >= 1.")

  if (isTRUE(cpp.option)) {
    stop(
      "cpp.option = TRUE is not yet supported in the dppca package. ",
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

#' Gaussian mechanism applied to spherical Kendall matrix
#'
#' @param X Numeric matrix.
#' @param sig Noise standard deviation.
#' @param cpp.option Logical; currently only FALSE is supported.
#'
#' @return A noisy symmetric d x d matrix.
#' @keywords internal
mech_tau_sph <- function(X, sig, cpp.option = FALSE) {
  X <- as.matrix(X)
  n <- nrow(X)
  d <- ncol(X)

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
