# ============================================================
# helper_scree.R
# Internal helper functions for scree-related routines
# ============================================================

#' Validate common inputs for scree-related routines
#'
#' @description
#' Internal helper used at the beginning of scree-related routines to validate
#' the most common user inputs. In particular, this function checks that:
#' \itemize{
#'   \item \code{X} can be coerced to a numeric matrix,
#'   \item the sample size and ambient dimension are valid,
#'   \item \code{k} is a valid number of leading components,
#'   \item \code{eps_total} is a positive scalar, and
#'   \item \code{delta_total} is a scalar in \code{(0, 1)}.
#' }
#' This helper does not return a transformed object for downstream use; its role
#' is simply to fail early with informative error messages before private scree
#' estimation begins.
#'
#' @param X A matrix-like object containing the data, with observations in rows
#'   and variables in columns.
#' @param k Integer number of leading principal components for which scree values
#'   will be estimated.
#' @param eps_total Total privacy budget \eqn{\varepsilon} allocated to the
#'   scree-estimation routine.
#' @param delta_total Total privacy budget \eqn{\delta} allocated to the
#'   scree-estimation routine.
#'
#' @return Invisibly returns \code{TRUE} if all checks pass.
#' @keywords internal
validate_scree_inputs <- function(X, k, eps_total, delta_total) {
  X <- as.matrix(X)
  n <- nrow(X); d <- ncol(X)
  if (n < 2) stop("Need n >= 2.")
  if (d < 1) stop("Need ncol(X) >= 1.")
  if (!is.numeric(k) || length(k) != 1 || !is.finite(k) || k < 1 || k > d) stop("k must be an integer in {1, ..., ncol(X)}.")
  if (!is.numeric(eps_total) || length(eps_total) != 1 || !is.finite(eps_total) || eps_total <= 0) stop("eps_total must be a single positive number.")
  if (!is.numeric(delta_total) || length(delta_total) != 1 || !is.finite(delta_total) || delta_total <= 0 || delta_total >= 1) stop("delta_total must be a single number in (0, 1).")
  invisible(TRUE)
}

#' Enforce a nonnegative and nonincreasing scree sequence
#'
#' @description
#' Internal post-processing helper for scree estimates. Since scree values are
#' variances, they should be nonnegative, and in many workflows it is desirable
#' to enforce a monotone decreasing profile across components. This function
#' first truncates the input at zero and then applies isotonic regression to the
#' negated sequence so that the final output is nonincreasing.
#'
#' @param x Numeric vector of raw scree estimates.
#'
#' @return Numeric vector of the same length as \code{x}, containing a
#'   nonnegative and nonincreasing version of the input.
#' @keywords internal
scree_post_processing <- function(x) {
  y <- pmax(as.numeric(x), 0)
  fit <- stats::isoreg(seq_along(y), -y)
  pmax(-fit$yf, 0)
}

#' Convert scree values to explained variance ratios safely
#'
#' @description
#' Internal helper that converts a vector of scree estimates to explained
#' variance ratios (EVR). If the total scree is positive and finite, this returns
#' the normalized vector \code{scree / sum(scree)}. If the total scree is
#' nonpositive or nonfinite, this helper returns a zero vector instead of
#' propagating problematic values.
#'
#' @param scree Numeric vector of scree estimates.
#'
#' @return Numeric vector of the same length as \code{scree}, interpreted as
#'   explained variance ratios.
#' @keywords internal
scree_to_evr <- function(scree) {
  scree <- as.numeric(scree)
  s <- sum(scree)
  if (!is.finite(s) || s <= 0) return(rep(0, length(scree)))
  scree / s
}

#' Clamp numeric values to a fixed interval
#'
#' @description
#' Internal helper that truncates numeric values to the closed interval
#' \code{[lo, hi]}. This is used in several places where intermediate values must
#' be kept inside valid numerical or algorithmic bounds.
#'
#' @param x Numeric vector.
#' @param lo Finite lower bound.
#' @param hi Finite upper bound satisfying \code{lo <= hi}.
#'
#' @return Numeric vector of the same length as \code{x}, with all entries
#'   truncated to \code{[lo, hi]}.
#' @keywords internal
winsorization <- function(x, lo, hi) {
  if (!is.finite(lo) || !is.finite(hi) || lo > hi) stop("Need finite lo <= hi.")
  pmin(pmax(x, lo), hi)
}


#' Differentially private scree estimation via clipped mean
#'
#' @description
#' Internal helper that estimates the leading scree values using a clipped-mean
#' differentially private estimator applied to the centered squared scores of
#' each principal component.
#'
#' The procedure is as follows:
#' \enumerate{
#'   \item preprocess the data matrix for PCA,
#'   \item compute a non-private or private direction matrix,
#'   \item form the projected score matrix,
#'   \item for each component \eqn{\ell}, compute centered squared scores
#'         \eqn{w_{i\ell} = (y_{i\ell} - \bar y_\ell)^2},
#'   \item clip each \eqn{w_{i\ell}} at \code{C_clip},
#'   \item release the clipped mean with Gaussian noise,
#'   \item apply the \eqn{n/(n-1)} correction to obtain the final scree estimate.
#' }
#'
#' If \code{mono = TRUE}, the final scree vector is post-processed so that it is
#' nonnegative and nonincreasing.
#'
#' @param X Numeric data matrix with observations in rows.
#' @param k Integer number of leading principal components.
#' @param eps_total Total privacy epsilon allocated to the scree routine.
#' @param delta_total Total privacy delta allocated to the scree routine.
#' @param center Logical; whether to center columns before PCA.
#' @param standardize Logical; whether to standardize columns before PCA.
#' @param C_clip Positive clipping threshold applied to centered squared scores.
#' @param g_dppca Logical; whether to privatize the PCA direction matrix.
#' @param cpp.option Logical passed to \code{mech_tau_sph()} when private
#'   directions are computed.
#' @param mono Logical; whether to post-process the final scree vector so that
#'   it is nonnegative and nonincreasing.
#'
#' @return A list with components:
#' \describe{
#'   \item{\code{scree}}{Private scree estimates.}
#'   \item{\code{scree_np}}{Non-private clipped scree estimates.}
#'   \item{\code{evr}}{Explained variance ratios computed from the private scree
#'     estimates.}
#' }
#' @keywords internal
dp_scree_clipped <- function(X, k, eps_total, delta_total,
                             center = TRUE, standardize = TRUE,
                             C_clip,
                             g_dppca = FALSE, cpp.option = FALSE,
                             mono = TRUE) {
  validate_scree_inputs(
    X = X,
    k = k,
    eps_total = eps_total,
    delta_total = delta_total
  )
  
  X <- as.matrix(X)
  n <- nrow(X)
  k <- as.integer(k)
  
  if (!is.numeric(C_clip) || length(C_clip) != 1 ||
      !is.finite(C_clip) || C_clip <= 0) {
    stop("C_clip must be a single positive number.")
  }
  
  if (isTRUE(g_dppca)) {
    eps_dir <- eps_total / 2
    delta_dir <- delta_total / 2
    eps_scree <- eps_total / 2
    delta_scree <- delta_total / 2
  } else {
    eps_dir <- NULL
    delta_dir <- NULL
    eps_scree <- eps_total
    delta_scree <- delta_total
  }
  
  X_proc <- prep_matrix_for_pca(
    X = X,
    center = center,
    standardize = standardize
  )
  
  V_used <- compute_pc_dir(
    X_proc = X_proc,
    k = k,
    g_dppca = g_dppca,
    eps_dir = eps_dir,
    delta_dir = delta_dir,
    cpp.option = cpp.option
  )
  
  Y <- X_proc %*% V_used
  
  eps_ell <- eps_scree / k
  delta_ell <- delta_scree / k
  
  scree_np <- numeric(k)
  scree <- numeric(k)
  
  for (ell in seq_len(k)) {
    y <- Y[, ell]
    ybar <- mean(y)
    w <- (y - ybar)^2
    
    w_clip <- pmin(w, C_clip)
    mu_hat <- mean(w_clip)
    
    scree_np[ell] <- (n / (n - 1)) * mu_hat
    
    Delta_ell <- C_clip / (n - 1)
    sd_noise <- Delta_ell * sqrt(2 * log(1.25 / delta_ell)) / eps_ell
    
    scree[ell] <- max(
      scree_np[ell] + stats::rnorm(1, mean = 0, sd = sd_noise),
      0
    )
  }
  
  if (isTRUE(mono)) {
    scree <- scree_post_processing(scree)
  }
  
  list(
    scree = scree,
    scree_np = scree_np,
    evr = scree_to_evr(scree)
  )
}


#' Private scale proxy from a noisy dyadic histogram
#'
#' @description
#' Internal helper used in the Huber scree routine. Starting from a nonnegative
#' scalar sample \code{u}, this function places values into dyadic bins indexed by
#' powers of two, perturbs the resulting histogram counts with Laplace noise, and
#' returns the dyadic scale corresponding to the largest noisy bin.
#'
#' @param u Numeric nonnegative vector whose scale is to be summarized.
#' @param eps_m2 Positive privacy epsilon used for the noisy histogram step.
#' @param k_min_m2 Integer lower bound for the dyadic bin index.
#' @param k_max_m2 Integer upper bound for the dyadic bin index.
#'
#' @return A nonnegative numeric scalar corresponding to the selected dyadic
#'   scale level.
#' @keywords internal
dp_hist_m2 <- function(u, eps_m2, k_min_m2 = -20, k_max_m2 = 40) {
  u <- as.numeric(u)
  
  if (length(u) < 1) stop("u must have length >= 1.")
  if (!is.finite(eps_m2) || eps_m2 <= 0) stop("eps_m2 must be > 0.")
  
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
  
  noisy <- as.numeric(counts) + VGAM::rlaplace(length(counts), scale = 1 / eps_m2)
  noisy <- pmax(noisy, 0)
  
  2^as.integer(names(counts)[which.max(noisy)])
}

#' Private scalar scale proxy from paired differences
#'
#' @description
#' Internal helper used by the Huber scree estimator to obtain a private scale
#' proxy for a one-dimensional sample. The input vector \code{w} is first paired,
#' transformed into squared paired differences, aggregated by block medians, and
#' then fed into a noisy dyadic histogram.
#'
#' @param w Numeric vector, typically the centered squared projected scores for a
#'   single principal component.
#' @param eps_m2 Positive privacy epsilon used for the scale-proxy step.
#' @param M Optional integer number of blocks used in the block-median step. If
#'   \code{NULL}, a default based on \code{sqrt(n / 2)} is used.
#' @param k_min_m2 Integer lower bound for the dyadic bin index.
#' @param k_max_m2 Integer upper bound for the dyadic bin index.
#'
#' @return A nonnegative numeric scalar private scale proxy.
#' @keywords internal
dp_m2 <- function(w, eps_m2, M = NULL, k_min_m2 = -20, k_max_m2 = 40) {
  w <- as.numeric(w)
  n <- length(w)
  
  if (n < 4) stop("Need length(w) >= 4.")
  if (!is.finite(eps_m2) || eps_m2 <= 0) stop("eps_m2 must be > 0.")
  
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

#' Convert a private scale proxy into a Huber threshold
#'
#' @description
#' Internal helper used by the Huber scree routine. Given a private scale proxy
#' \code{m2_hat}, this function produces a scalar Huber robustification threshold.
#'
#' @param m2_hat Nonnegative private scale proxy.
#' @param eps_tau Positive privacy epsilon allocated to the Huber
#'   noisy-gradient-descent step for one component.
#' @param n_tau Effective sample size used in the threshold calculation.
#'
#' @return Positive numeric scalar Huber threshold.
#' @keywords internal
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

#' DP Huber noisy gradient descent for a scalar mean
#'
#' @description
#' Internal helper that implements a scalar Huber-type differentially private
#' mean estimator via noisy gradient descent. The optional argument
#' \code{scale_factor} is useful when the released quantity should be a scaled
#' version of the mean, such as the \eqn{n/(n-1)}-adjusted scree convention.
#'
#' @param w Numeric vector representing the one-dimensional sample.
#' @param eps_gd Positive privacy epsilon allocated to the noisy-gradient-descent
#'   stage.
#' @param delta_gd Positive privacy delta allocated to the noisy-gradient-descent
#'   stage.
#' @param tau Positive Huber threshold.
#' @param T Positive integer number of gradient-descent iterations.
#' @param mu0 Initial value for the iterative procedure.
#' @param eta0 Base step size.
#' @param step_schedule Character string specifying the step-size schedule.
#'   Allowed values are \code{"fixed"} and \code{"1/sqrt(t)"}.
#' @param scale_factor Positive multiplicative factor applied to the released
#'   estimate and to the noise calibration.
#'
#' @return Numeric scalar representing the final private estimate.
#' @keywords internal
dp_huber_noisy_gd <- function(w, eps_gd, delta_gd, tau, T, mu0 = 0, eta0 = 1) {
  w <- as.numeric(w)
  n <- length(w)
  
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
  
  eps_step <- eps_gd / T
  del_step <- delta_gd / T
  
  for (t in 0:(T - 1L)) {
    r <- w - mu
    psi <- winsorization(r, -tau, tau)
    g <- mean(psi)
    
    Delta_step <- (2 * eta0 * tau) / n
    sd_noise <- Delta_step * sqrt(2 * log(1.25 / del_step)) / eps_step
    
    mu <- mu + eta0 * g + stats::rnorm(1, mean = 0, sd = sd_noise)
  }
  
  mu
}


#' Differentially private scree estimation via Huber-type robust mean
#'
#' @description
#' Internal helper that estimates the leading scree values using a Huber-type
#' differentially private mean estimator applied to the centered squared scores
#' of each principal component.
#'
#' The procedure is as follows:
#' \enumerate{
#'   \item preprocess the data matrix for PCA,
#'   \item compute a non-private or private direction matrix,
#'   \item form the projected score matrix,
#'   \item for each component \eqn{\ell}, compute centered squared scores
#'         \eqn{w_{i\ell} = (y_{i\ell} - \bar y_\ell)^2},
#'   \item privately estimate a scale proxy from \eqn{w_{1\ell}, \dots, w_{n\ell}}
#'         using \code{dp_m2()},
#'   \item convert that scale proxy into a Huber threshold using
#'         \code{tau_from_m2()},
#'   \item run noisy gradient descent via \code{dp_huber_noisy_gd()} to obtain
#'         a private mean estimate of \eqn{w_{i\ell}},
#'   \item apply the \eqn{n/(n-1)} correction to obtain the final scree estimate.
#' }
#'
#' If \code{mono = TRUE}, the final scree vector is post-processed so that it is
#' nonnegative and nonincreasing.
#'
#' @param X Numeric data matrix with observations in rows.
#' @param k Integer number of leading principal components.
#' @param eps_total Total privacy epsilon allocated to the scree routine.
#' @param delta_total Total privacy delta allocated to the scree routine.
#' @param g_dppca Logical; whether to privatize the PCA direction matrix.
#' @param cpp.option Logical passed to \code{mech_tau_sph()} when private
#'   directions are computed.
#' @param center Logical; whether to center columns before PCA.
#' @param standardize Logical; whether to standardize columns before PCA.
#' @param mu0 Initial value for noisy gradient descent.
#' @param eta0 Fixed step size used in noisy gradient descent.
#' @param T Optional integer number of gradient descent iterations. If
#'   \code{NULL}, a default based on \code{ceiling(log(n))} is used.
#' @param M Optional integer number of blocks used in \code{dp_m2()}. If
#'   \code{NULL}, a default based on \code{floor(sqrt(n / 2))} is used.
#' @param k_min_m2 Integer lower bound for dyadic histogram bins used in
#'   \code{dp_hist_m2()}.
#' @param k_max_m2 Integer upper bound for dyadic histogram bins used in
#'   \code{dp_hist_m2()}.
#' @param m2_frac Fraction of the scree privacy budget allocated to the private
#'   scale-proxy step \code{dp_m2()}.
#' @param mono Logical; whether to post-process the final scree vector so that
#'   it is nonnegative and nonincreasing.
#'
#' @return A list with components:
#' \describe{
#'   \item{\code{scree}}{Private scree estimates.}
#'   \item{\code{scree_np}}{Non-private scree estimates computed as
#'     \eqn{(n/(n-1)) \times mean(w_{i\ell})}.}
#'   \item{\code{evr}}{Explained variance ratios computed from the private scree
#'     estimates.}
#' }
#' @keywords internal
dp_scree_huber <- function(X, k, eps_total, delta_total,
                           g_dppca = FALSE, cpp.option = FALSE,
                           center = TRUE, standardize = FALSE,
                           mu0 = 0, eta0 = 1, T = NULL, M = NULL,
                           k_min_m2 = -20, k_max_m2 = 40,
                           m2_frac = 1/4,
                           mono = TRUE) {
  validate_scree_inputs(X = X, k = k, eps_total = eps_total, delta_total = delta_total)
  
  X <- as.matrix(X)
  n <- nrow(X)
  k <- as.integer(k)
  
  if (!is.finite(eta0) || eta0 <= 0) {
    stop("eta0 must be > 0.")
  }
  if (!is.finite(m2_frac) || m2_frac <= 0 || m2_frac >= 1) {
    stop("m2_frac must be in (0, 1).")
  }
  
  if (is.null(T)) T <- ceiling(log(n))
  T <- max(1L, as.integer(T))
  
  if (is.null(M)) M <- floor(sqrt(n / 2))
  M <- max(1L, as.integer(M))
  
  if (isTRUE(g_dppca)) {
    eps_dir <- eps_total / 2
    delta_dir <- delta_total / 2
    eps_scree <- eps_total / 2
    delta_scree <- delta_total / 2
  } else {
    eps_dir <- NULL
    delta_dir <- NULL
    eps_scree <- eps_total
    delta_scree <- delta_total
  }
  
  eps_m2 <- eps_scree * m2_frac
  delta_m2 <- delta_scree * m2_frac
  
  eps_gd <- eps_scree * (1 - m2_frac)
  delta_gd <- delta_scree * (1 - m2_frac)
  
  eps_m2_ell <- eps_m2 / k
  delta_m2_ell <- delta_m2 / k
  
  eps_gd_ell <- eps_gd / k
  delta_gd_ell <- delta_gd / k
  
  X_proc <- prep_matrix_for_pca(X = X, center = center, standardize = standardize)
  
  V_used <- compute_pc_dir(
    X_proc = X_proc,
    k = k,
    g_dppca = g_dppca,
    eps_dir = eps_dir,
    delta_dir = delta_dir,
    cpp.option = cpp.option
  )
  
  Y <- X_proc %*% V_used
  
  scree_np <- numeric(k)
  scree <- numeric(k)
  
  for (ell in seq_len(k)) {
    y <- Y[, ell]
    ybar <- mean(y)
    w <- (y - ybar)^2
    
    m2_hat <- dp_m2(
      w = w,
      eps_m2 = eps_m2_ell,
      M = M,
      k_min_m2 = k_min_m2,
      k_max_m2 = k_max_m2
    )
    
    tau_ell <- tau_from_m2(
      m2_hat = m2_hat,
      eps_tau = eps_gd_ell,
      n_tau = n
    )
    
    if (!is.finite(tau_ell) || tau_ell <= 0) {
      tau_ell <- 1
    }
    
    muT <- dp_huber_noisy_gd(
      w = w,
      eps_gd = eps_gd_ell,
      delta_gd = delta_gd_ell,
      tau = tau_ell,
      T = T,
      mu0 = mu0,
      eta0 = eta0
    )
    
    scree_np[ell] <- (n / (n - 1)) * mean(w)
    scree[ell] <- (n / (n - 1)) * max(muT, 0)
  }
  
  if (isTRUE(mono)) {
    scree <- scree_post_processing(scree)
  }
  
  list(
    scree = scree,
    scree_np = scree_np,
    evr = scree_to_evr(scree)
  )
}


#' Unbounded private upper quantile via log-binning
#'
#' @description
#' Internal helper that computes a private upper-tail quantile estimate using a
#' logarithmic grid together with noisy thresholding and noisy cumulative counts.
#'
#' This function is used in the PMW scree routine to estimate an upper clipping
#' bound privately. The search is carried out over a geometric grid determined by
#' the lower anchor \code{l} and the log-binning base \code{beta}.
#'
#' @param data Numeric vector.
#' @param l Finite lower anchor for the search grid.
#' @param beta Log-binning base used in the geometric grid. Must satisfy
#'   \code{beta > 1}.
#' @param q Quantile level in \code{(0, 1)}.
#' @param eps_1 Privacy epsilon for the noisy threshold step.
#' @param delta_1 Privacy delta for the noisy threshold step.
#' @param eps_2 Privacy epsilon for the noisy cumulative scan.
#' @param delta_2 Privacy delta for the noisy cumulative scan.
#' @param max_extra_bins Nonnegative integer giving the extra number of bins
#'   searched past the largest occupied bin.
#'
#' @return Numeric scalar private upper-quantile estimate.
#' @keywords internal
unbounded_quantile_upper <- function(data, l, beta, q,
                                     eps_1, delta_1,
                                     eps_2, delta_2,
                                     max_extra_bins = 1000) {
  data <- as.numeric(data)
  n <- length(data)
  
  if (n < 1) stop("data must have length >= 1.")
  if (!is.finite(l)) stop("l must be finite.")
  if (!is.finite(beta) || beta <= 1) stop("beta must be > 1.")
  if (!is.finite(q) || q <= 0 || q >= 1) stop("q must be in (0, 1).")
  if (!is.finite(eps_1) || eps_1 <= 0) stop("eps_1 must be > 0.")
  if (!is.finite(delta_1) || delta_1 <= 0 || delta_1 >= 1) {
    stop("delta_1 must be in (0, 1).")
  }
  if (!is.finite(eps_2) || eps_2 <= 0) stop("eps_2 must be > 0.")
  if (!is.finite(delta_2) || delta_2 <= 0 || delta_2 >= 1) {
    stop("delta_2 must be in (0, 1).")
  }
  
  max_extra_bins <- as.integer(max_extra_bins)
  if (!is.finite(max_extra_bins) || max_extra_bins < 0) {
    stop("max_extra_bins must be a nonnegative integer.")
  }
  
  sd1 <- sqrt(2 * log(1.25 / delta_1)) / eps_1
  sd2 <- sqrt(2 * log(1.25 / delta_2)) / eps_2
  
  z <- pmax(data - l + 1, beta)
  idx <- as.integer(floor(log(z, base = beta)))
  counts <- table(idx)
  
  get_count <- function(i) {
    key <- as.character(i)
    if (key %in% names(counts)) as.integer(counts[[key]]) else 0L
  }
  
  t_noisy <- q * n + stats::rnorm(1, mean = 0, sd = sd1)
  
  cur <- 0L
  i <- 0L
  max_bin <- if (length(idx) > 0) max(idx) else 0L
  i_max <- max_bin + max_extra_bins
  
  repeat {
    cur <- cur + get_count(i)
    i <- i + 1L
    
    if (cur + stats::rnorm(1, mean = 0, sd = sd2) > t_noisy) break
    if (i > i_max) break
  }
  
  val <- tryCatch(beta^i + l - 1, error = function(e) .Machine$double.xmax)
  if (!is.finite(val)) .Machine$double.xmax else val
}

#' Unbounded private quantile via log-binning
#'
#' @description
#' Internal wrapper around \code{unbounded_quantile_upper()} that handles both
#' upper and lower quantiles on a finite target interval \code{[l, u]}.
#'
#' If \code{q < 0.5}, the lower-tail quantile is computed by sign-flipping the
#' data and calling the upper-tail routine.
#'
#' @param data Numeric vector.
#' @param l Finite lower truncation bound.
#' @param u Finite upper truncation bound.
#' @param beta Log-binning base used in the geometric grid. Must satisfy
#'   \code{beta > 1}.
#' @param q Quantile level in \code{(0, 1)}.
#' @param eps_1 Privacy epsilon for the noisy threshold step.
#' @param delta_1 Privacy delta for the noisy threshold step.
#' @param eps_2 Privacy epsilon for the noisy cumulative scan.
#' @param delta_2 Privacy delta for the noisy cumulative scan.
#' @param max_extra_bins Nonnegative integer giving the extra number of bins
#'   searched past the largest occupied bin.
#'
#' @return Numeric scalar private quantile estimate.
#' @keywords internal
unbounded_quantile <- function(data, l, u, beta, q,
                               eps_1, delta_1,
                               eps_2, delta_2,
                               max_extra_bins = 1000) {
  data <- as.numeric(data)
  
  if (!is.finite(l) || !is.finite(u) || l > u) stop("Need finite l <= u.")
  if (!is.finite(beta) || beta <= 1) stop("beta must be > 1.")
  if (!is.finite(q) || q <= 0 || q >= 1) stop("q must be in (0, 1).")
  
  if (q < 0.5) {
    est <- -unbounded_quantile_upper(
      data = -data,
      l = -u,
      beta = beta,
      q = 1 - q,
      eps_1 = eps_1,
      delta_1 = delta_1,
      eps_2 = eps_2,
      delta_2 = delta_2,
      max_extra_bins = max_extra_bins
    )
    return(max(est, l))
  }
  
  est <- unbounded_quantile_upper(
    data = data,
    l = l,
    beta = beta,
    q = q,
    eps_1 = eps_1,
    delta_1 = delta_1,
    eps_2 = eps_2,
    delta_2 = delta_2,
    max_extra_bins = max_extra_bins
  )
  
  min(est, u)
}

#' Differentially private scree estimation via PMW
#'
#' @description
#' Internal helper that estimates the leading scree values using a PMW-style
#' procedure based on:
#' \enumerate{
#'   \item PCA preprocessing and direction estimation,
#'   \item centered squared scores for each component,
#'   \item private lower and upper quantile estimation,
#'   \item winsorization of the squared scores,
#'   \item Gaussian release of the winsorized mean,
#'   \item the \eqn{n/(n-1)} scree correction.
#' }
#'
#' If \code{split_mode = TRUE}, one subset of the sample is used for private
#' quantile estimation and the other subset is used for the winsorized mean step.
#' If \code{split_mode = FALSE}, the full sample is reused in both steps.
#'
#' The clipping proportion is the practical choice
#' \code{max(trim_const / n_q, eta)}, truncated above at \code{0.49}.
#'
#' @param X Numeric data matrix with observations in rows.
#' @param k Integer number of leading principal components.
#' @param eps_total Total privacy epsilon allocated to the scree routine.
#' @param delta_total Total privacy delta allocated to the scree routine.
#' @param g_dppca Logical; whether to privatize the PCA direction matrix.
#' @param cpp.option Logical passed to \code{mech_tau_sph()} when private
#'   directions are computed.
#' @param split_mode Logical; if \code{TRUE}, split the sample into quantile and
#'   mean subsets.
#' @param center Logical; whether to center columns before PCA.
#' @param standardize Logical; whether to standardize columns before PCA.
#' @param beta Log-binning base used in the private quantile estimator.
#' @param a Finite lower support bound supplied to the private quantile routine.
#' @param b Finite upper support bound supplied to the private quantile routine.
#' @param trim_const Positive constant controlling the practical clipping
#'   proportion \code{max(trim_const / n_q, eta)}.
#' @param eta Lower bound in the practical clipping proportion.
#' @param mono Logical; whether to post-process the final scree vector so that
#'   it is nonnegative and nonincreasing.
#' @param max_extra_bins Nonnegative integer controlling how far the unbounded
#'   quantile scan extends past the largest occupied bin.
#'
#' @return A list with components:
#' \describe{
#'   \item{\code{scree}}{Private scree estimates.}
#'   \item{\code{scree_np}}{Non-private winsorized scree estimates.}
#'   \item{\code{evr}}{Explained variance ratios computed from the private
#'   scree estimates.}
#' }
#' @keywords internal
dp_scree_pmw <- function(X, k, eps_total, delta_total,
                         g_dppca = FALSE, cpp.option = FALSE,
                         split_mode = TRUE,
                         center = TRUE, standardize = TRUE,
                         beta = 1.01, a, b,
                         trim_const = 10, eta = 0.01,
                         mono = TRUE, max_extra_bins = 1000) {
  validate_scree_inputs(
    X = X,
    k = k,
    eps_total = eps_total,
    delta_total = delta_total
  )
  
  X <- as.matrix(X)
  n <- nrow(X)
  k <- as.integer(k)
  
  if (!is.finite(beta) || beta <= 1) stop("beta must be > 1.")
  if (!is.finite(a) || !is.finite(b) || a > b) stop("Need finite a <= b.")
  if (!is.numeric(trim_const) || length(trim_const) != 1 || trim_const <= 0) {
    stop("trim_const must be a single positive number.")
  }
  if (!is.numeric(eta) || length(eta) != 1 || eta <= 0 || eta >= 0.5) {
    stop("eta must be in (0, 0.5).")
  }
  
  if (isTRUE(g_dppca)) {
    eps_dir <- eps_total / 2
    delta_dir <- delta_total / 2
    eps_scree <- eps_total / 2
    delta_scree <- delta_total / 2
  } else {
    eps_dir <- NULL
    delta_dir <- NULL
    eps_scree <- eps_total
    delta_scree <- delta_total
  }
  
  eps_ell <- eps_scree / k
  delta_ell <- delta_scree / k
  
  eps_Q <- eps_ell / 4
  delta_Q <- delta_ell / 4
  eps_M <- eps_ell / 2
  delta_M <- delta_ell / 2
  
  eps_q1 <- eps_Q / 2
  delta_q1 <- delta_Q / 2
  eps_q2 <- eps_Q / 2
  delta_q2 <- delta_Q / 2
  
  X_proc <- prep_matrix_for_pca(
    X = X,
    center = center,
    standardize = standardize
  )
  
  V_used <- compute_pc_dir(
    X_proc = X_proc,
    k = k,
    g_dppca = g_dppca,
    eps_dir = eps_dir,
    delta_dir = delta_dir,
    cpp.option = cpp.option
  )
  
  Y <- X_proc %*% V_used
  
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
  
  n_q <- length(idx_q)
  n_m <- length(idx_m)
  
  trim_param <- min(max(trim_const / n_q, eta), 0.49)
  
  scree_np <- numeric(k)
  scree <- numeric(k)
  
  for (ell in seq_len(k)) {
    y <- Y[, ell]
    ybar <- mean(y)
    w <- (y - ybar)^2
    
    L <- unbounded_quantile(
      data = w[idx_q],
      l = a,
      u = b,
      beta = beta,
      q = trim_param,
      eps_1 = eps_q1,
      delta_1 = delta_q1,
      eps_2 = eps_q2,
      delta_2 = delta_q2,
      max_extra_bins = max_extra_bins
    )
    
    U <- unbounded_quantile(
      data = w[idx_q],
      l = a,
      u = b,
      beta = beta,
      q = 1 - trim_param,
      eps_1 = eps_q1,
      delta_1 = delta_q1,
      eps_2 = eps_q2,
      delta_2 = delta_q2,
      max_extra_bins = max_extra_bins
    )
    
    if (!is.finite(L) || !is.finite(U) || U < L) {
      L <- a
      U <- b
    }
    
    w_win <- pmin(pmax(w[idx_m], L), U)
    mu_hat <- mean(w_win)
    
    scree_np[ell] <- (n / (n - 1)) * mu_hat
    
    Delta_ell <- (n / (n - 1)) * (U - L) / n_m
    sd_noise <- Delta_ell * sqrt(2 * log(1.25 / delta_M)) / eps_M
    
    scree[ell] <- max(
      scree_np[ell] + stats::rnorm(1, mean = 0, sd = sd_noise),
      0
    )
  }
  
  if (isTRUE(mono)) {
    scree <- scree_post_processing(scree)
  }
  
  list(
    scree = scree,
    scree_np = scree_np,
    evr = scree_to_evr(scree)
  )
}