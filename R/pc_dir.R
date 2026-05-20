# ============================================================
# dp_pc_dir.R
# Differentially private principal component direction estimation
# ============================================================


#' Estimate principal component directions
#'
#' This function returns the principal component directions used by the plotting
#' and estimation functions in this package. By default, it computes the usual
#' non-private directions from the sample covariance matrix. If `g_dppca = TRUE`,
#' it computes private directions using the spherical-transformation version of
#' g-DPPCA proposed by \insertCite{kim2025robustdppca;textual}{dppca}: it adds
#' Gaussian noise to the spherical Kendall matrix and then extracts its leading
#' eigenvectors.
#'
#' @param X A numeric matrix or data frame. Rows correspond to observations and
#'  columns correspond to variables.
#' @param k Number of principal component directions to return. Must be an
#'   integer between `1` and the number of columns in `X`.
#' @param center A logical value indicating whether to center the columns of `X`
#'   before computing principal component directions. The default is `TRUE`.
#' @param standardize A logical value indicating whether to scale the columns of
#'   `X` by their sample standard deviations after optional centering. The
#'   default is `FALSE`.
#' @param g_dppca Whether to compute private principal component directions using
#'   the spherical Kendall mechanism based on the g-DPPCA method. If `FALSE`,
#'   the usual non-private directions are computed from the sample covariance
#'   matrix. The default is `FALSE`.
#' @param eps_dir Positive number defining the `epsilon` privacy parameter for
#'   private principal component directions. Required when `g_dppca = TRUE`.
#' @param delta_dir Number in `(0, 1)` defining the `delta` privacy parameter for
#'   private principal component directions. Required when `g_dppca = TRUE`.
#' @param cpp.option A logical value reserved for a future C++ implementation of
#'   the spherical Kendall matrix. Currently only `FALSE` is supported.
#'
#' @return A numeric matrix with `ncol(X)` rows and `k` columns. The columns are
#'   orthonormal principal component directions.
#'
#' @details
#' The non-private option computes leading eigenvectors of the sample covariance
#' matrix of the preprocessed data. The private option is based on the spherical
#' Kendall mechanism of \insertCite{kim2025robustdppca;textual}{dppca}: it first
#' forms the spherical Kendall matrix from pairwise normalized differences, adds
#' symmetric Gaussian noise, and then computes leading eigenvectors. The final
#' eigenvector matrix is re-orthonormalized by QR decomposition.
#' For a detailed procedure and mathematical formulations,
#' refer \url{https://yejinjo0220.github.io/dppca/articles/pc_direction}.
#'
#' @references
#' \insertRef{kim2025robustdppca}{dppca}
#'
#' @examples
#' data(gau, package = "dppca")
#'
#' # Use a small subset to keep the example fast.
#' X <- head(gau, 50)
#'
#' # Non-private principal component directions
#' V <- dp_pc_dir(X, k = 2)
#' head(V)
#'
#' # Private principal component directions
#' set.seed(123)
#' V_private <- dp_pc_dir(
#'   X,
#'   k = 2,
#'   g_dppca = TRUE,
#'   eps_dir = 1,
#'   delta_dir = 1e-2
#' )
#' head(V_private)
#'
#' # Generate a small low-rank dataset.
#' n <- 50
#' z1 <- rnorm(n)
#' z2 <- rnorm(n)
#' X <- cbind(
#'   x1 = z1 + 0.2 * rnorm(n),
#'   x2 = 0.8 * z1 + 0.2 * rnorm(n),
#'   x3 = z2 + 0.2 * rnorm(n),
#'   x4 = 0.5 * z1 - 0.4 * z2 + 0.2 * rnorm(n)
#' )
#'
#' # Non-private principal component directions
#' V <- dp_pc_dir(X, k = 2)
#' head(V)
#'
#' # Private principal component directions
#' V_private <- dp_pc_dir(
#'   X,
#'   k = 2,
#'   g_dppca = TRUE,
#'   eps_dir = 1,
#'   delta_dir = 1e-2
#' )
#' head(V_private)
#'
#' @export
#'
dp_pc_dir <- function(X,
                      k,
                      center = TRUE,
                      standardize = FALSE,
                      g_dppca = FALSE,
                      eps_dir = NULL,
                      delta_dir = NULL,
                      cpp.option = FALSE) {
  X_proc <- prep_matrix_for_pca(X, center = center,standardize = standardize)

  k <- validate_pc_dir_k(k, ncol(X_proc))
  validate_logical_value(g_dppca, "g_dppca")
  validate_logical_value(cpp.option, "cpp.option")

  n <- nrow(X_proc)

  if (g_dppca) {
    eps_dir <- validate_positive_number(eps_dir, "eps_dir")
    delta_dir <- validate_probability(delta_dir, "delta_dir")

    sigma_sph <- 2 * sqrt(2 * log(1.25 / delta_dir)) / (n * eps_dir)
    pca_matrix <- mech_tau_sph(
      X_proc,
      sig = sigma_sph,
      cpp.option = cpp.option
    )
  } else {
    pca_matrix <- stats::cov(X_proc)
  }

  V <- leading_eigenvectors(pca_matrix, k = k)
  orthonormalize_dir(V, k = k)
}

# Internal helpers ------------------------------------------------------------

#' Preprocess a matrix for principal component analysis
#'
#' @param X A numeric matrix or data frame.
#' @param center A logical value indicating whether to center columns.
#' @param standardize A logical value indicating whether to scale columns.
#'
#' @return A numeric matrix after preprocessing.
#' @noRd
prep_matrix_for_pca <- function(X, center = TRUE, standardize = FALSE) {
  validate_logical_value(center, "center")
  validate_logical_value(standardize, "standardize")

  if (is.data.frame(X)) {
    numeric_cols <- vapply(X, is.numeric, logical(1))
    if (!all(numeric_cols)) {
      stop("All columns of `X` must be numeric.", call. = FALSE)
    }
  }

  X <- as.matrix(X)

  if (!is.numeric(X)) {
    stop("`X` must be a numeric matrix or data frame.", call. = FALSE)
  }

  storage.mode(X) <- "double"

  if (nrow(X) < 2) {
    stop("`X` must have at least two rows.", call. = FALSE)
  }
  if (ncol(X) < 1) {
    stop("`X` must have at least one column.", call. = FALSE)
  }
  if (anyNA(X) || any(!is.finite(X))) {
    stop("`X` must contain only finite, non-missing values.", call. = FALSE)
  }

  X_proc <- X

  if (center) {
    X_proc <- scale(X_proc, center = TRUE, scale = FALSE)
  }

  if (standardize) {
    sds <- apply(X_proc, 2, stats::sd)
    if (any(!is.finite(sds)) || any(sds <= 0)) {
      stop(
        "All columns of `X` must have positive finite standard deviations ",
        "when `standardize = TRUE`.",
        call. = FALSE
      )
    }
    X_proc <- sweep(X_proc, 2, sds, "/")
  }

  X_proc
}

#' Validate the number of requested principal component directions
#'
#' @param k Number of directions.
#' @param d Number of variables.
#'
#' @return A positive integer scalar.
#' @noRd
validate_pc_dir_k <- function(k, d) {
  if (!is.numeric(k) || length(k) != 1 || !is.finite(k)) {
    stop("`k` must be a single positive integer.", call. = FALSE)
  }

  k_int <- as.integer(k)
  if (k != k_int) {
    stop("`k` must be a single positive integer.", call. = FALSE)
  }
  if (k_int < 1 || k_int > d) {
    stop("`k` must satisfy `1 <= k <= ncol(X)` after preprocessing.", call. = FALSE)
  }

  k_int
}

#' Validate a logical value
#'
#' @param x Object to validate.
#' @param arg Argument name.
#'
#' @return A logical value.
#' @noRd
validate_logical_value <- function(x, arg) {
  if (!is.logical(x) || length(x) != 1 || is.na(x)) {
    stop("`", arg, "` must be a single `TRUE` or `FALSE` value.", call. = FALSE)
  }

  x
}

#' Validate a positive numeric scalar
#'
#' @param x Object to validate.
#' @param arg Argument name.
#'
#' @return A positive numeric scalar.
#' @noRd
validate_positive_number <- function(x, arg) {
  if (!is.numeric(x) || length(x) != 1 || !is.finite(x) || x <= 0) {
    stop("`", arg, "` must be a single positive number.", call. = FALSE)
  }

  x
}

#' Validate a probability strictly between zero and one
#'
#' @param x Object to validate.
#' @param arg Argument name.
#'
#' @return A numeric scalar in `(0, 1)`.
#' @noRd
validate_probability <- function(x, arg) {
  if (!is.numeric(x) || length(x) != 1 || !is.finite(x) || x <= 0 || x >= 1) {
    stop("`", arg, "` must be a single number in `(0, 1)`.", call. = FALSE)
  }

  x
}

#' Compute leading eigenvectors of a symmetric matrix
#'
#' @param A A numeric symmetric matrix.
#' @param k Number of leading eigenvectors.
#'
#' @return A numeric matrix with `k` columns.
#' @noRd
leading_eigenvectors <- function(A, k) {
  A <- as.matrix(A)
  d <- nrow(A)

  if (ncol(A) != d) {
    stop("`A` must be a square matrix.", call. = FALSE)
  }
  if (!is.numeric(A) || anyNA(A) || any(!is.finite(A))) {
    stop("`A` must be a finite numeric matrix.", call. = FALSE)
  }

  # `rARPACK::eigs_sym()` is useful for partial decompositions, but base
  # `eigen()` is safer for small matrices and when all eigenvectors are needed.
  if (k >= d - 1) {
    eig <- eigen(A, symmetric = TRUE)
    return(as.matrix(eig$vectors[, seq_len(k), drop = FALSE]))
  }

  if (requireNamespace("rARPACK", quietly = TRUE)) {
    out <- rARPACK::eigs_sym(A, k = k, which = "LA")$vectors
    return(as.matrix(out))
  }

  eig <- eigen(A, symmetric = TRUE)
  as.matrix(eig$vectors[, seq_len(k), drop = FALSE])
}

#' Re-orthonormalize principal component directions
#'
#' @param V A numeric matrix of candidate directions.
#' @param k Number of columns to retain.
#'
#' @return A numeric matrix with orthonormal columns.
#' @noRd
orthonormalize_dir <- function(V, k) {
  V <- as.matrix(V)

  if (!is.numeric(V) || anyNA(V) || any(!is.finite(V))) {
    stop("`V` must be a finite numeric matrix.", call. = FALSE)
  }
  if (ncol(V) < k || k < 1) {
    stop("`k` must satisfy `1 <= k <= ncol(V)`.", call. = FALSE)
  }

  Q <- qr.Q(qr(V))
  Q[, seq_len(k), drop = FALSE]
}

#' Compute the Euclidean norm
#'
#' @param x A numeric vector.
#'
#' @return A nonnegative numeric scalar.
#' @noRd
norm2 <- function(x) {
  x <- as.numeric(x)
  sqrt(sum(x^2))
}

#' Normalize a vector to unit Euclidean norm
#'
#' @param x A numeric vector.
#'
#' @return A numeric vector with unit Euclidean norm.
#' @noRd
normalize <- function(x) {
  x <- as.numeric(x)
  nx <- norm2(x)

  if (!is.finite(nx) || nx <= sqrt(.Machine$double.eps)) {
    stop(
      "Cannot normalize a vector with zero or non-finite norm. ",
      "The spherical Kendall mechanism is undefined for duplicate observations.",
      call. = FALSE
    )
  }

  x / nx
}

#' Convert a packed symmetric vector into a symmetric matrix
#'
#' @param v A numeric vector of length `p * (p + 1) / 2`.
#' @param p Matrix dimension.
#'
#' @return A symmetric `p` by `p` matrix.
#' @noRd
vec2mat <- function(v, p) {
  v <- as.numeric(v)
  p <- as.integer(p)

  if (length(p) != 1 || !is.finite(p) || p < 1) {
    stop("`p` must be a positive integer.", call. = FALSE)
  }

  q <- length(v)
  if (q != p * (p + 1) / 2) {
    stop("`length(v)` is not compatible with `p`.", call. = FALSE)
  }

  out <- matrix(0, p, p)
  diag(out) <- v[seq_len(p)]

  if (q > p) {
    out[lower.tri(out)] <- v[(p + 1):q] / sqrt(2)
    out[upper.tri(out)] <- t(out)[upper.tri(out)]
  }

  out
}

#' Compute the spherical Kendall matrix
#'
#' @param X A numeric matrix.
#' @param cpp.option A logical value. Currently only `FALSE` is supported.
#'
#' @return A symmetric matrix.
#' @noRd
tau_sph <- function(X, cpp.option = FALSE) {
  X <- as.matrix(X)
  validate_logical_value(cpp.option, "cpp.option")

  n <- nrow(X)
  d <- ncol(X)

  if (!is.numeric(X) || anyNA(X) || any(!is.finite(X))) {
    stop("`X` must be a finite numeric matrix.", call. = FALSE)
  }
  if (n < 2) {
    stop("`X` must have at least two rows.", call. = FALSE)
  }
  if (d < 1) {
    stop("`X` must have at least one column.", call. = FALSE)
  }

  if (cpp.option) {
    stop(
      "`cpp.option = TRUE` is not yet supported. ",
      "Please use `cpp.option = FALSE`.",
      call. = FALSE
    )
  }

  hK <- matrix(0, d, d)

  for (i in seq_len(n - 1)) {
    for (j in (i + 1):n) {
      a <- X[j, ] - X[i, ]
      a <- normalize(a)
      hK <- hK + tcrossprod(a)
    }
  }

  hK * (2 / (n * (n - 1)))
}

#' Add Gaussian noise to the spherical Kendall matrix
#'
#' @param X A numeric matrix.
#' @param sig Noise standard deviation.
#' @param cpp.option A logical value. Currently only `FALSE` is supported.
#'
#' @return A noisy symmetric matrix.
#' @noRd
mech_tau_sph <- function(X, sig, cpp.option = FALSE) {
  X <- as.matrix(X)
  sig <- validate_positive_number(sig, "sig")

  d <- ncol(X)
  q <- d * (d + 1) / 2

  hK <- tau_sph(X, cpp.option = cpp.option)
  zeta <- stats::rnorm(q, mean = 0, sd = sig)
  zeta_mat <- vec2mat(zeta, d)

  hK + zeta_mat
}
