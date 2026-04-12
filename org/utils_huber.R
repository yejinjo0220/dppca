#' Generate Laplace noise
#'
#' Internal helper for Laplace(0, scale) sampling.
#'
#' @param n Number of samples.
#' @param scale Positive scale parameter.
#'
#' @return Numeric vector of length \code{n}.
#' @keywords internal
rlaplace <- function(n, scale) {
  n <- as.integer(n)
  if (!is.finite(n) || n < 0) stop("n must be a nonnegative integer.")
  if (!is.finite(scale) || scale <= 0) stop("scale must be > 0.")

  u <- stats::runif(n, -0.5, 0.5)
  -scale * sign(u) * log(1 - 2 * abs(u))
}

#' Clamp values to an interval
#'
#' @param x Numeric vector.
#' @param lo Lower bound.
#' @param hi Upper bound.
#'
#' @return Numeric vector with all values truncated to \code{[lo, hi]}.
#' @keywords internal
clamp <- function(x, lo, hi) {
  if (!is.finite(lo) || !is.finite(hi) || lo > hi) {
    stop("Need finite lo <= hi.")
  }
  pmin(pmax(x, lo), hi)
}

#' DP histogram argmax on log2 bins
#'
#' Internal helper used to estimate a scale proxy \code{m2_hat} from block
#' medians via a Laplace-perturbed histogram on dyadic bins.
#'
#' @param u Numeric nonnegative vector of block medians.
#' @param eps_hist Positive privacy epsilon for the histogram release.
#' @param k_min_star,k_max_star Integer range of dyadic bins.
#' @param include_zero_bin Logical; whether to include a separate zero bin.
#'
#' @return A list containing \code{k_star}, \code{m2_hat}, \code{noisy}, and \code{counts}.
#' @keywords internal
dp_hist_argmax_laplace_log2bins <- function(
    u,
    eps_hist,
    k_min_star = -20L,
    k_max_star = 40L,
    include_zero_bin = TRUE
) {
  u <- as.numeric(u)
  M <- length(u)

  if (M < 1) stop("u must have length >= 1.")
  if (!is.finite(eps_hist) || eps_hist <= 0) stop("eps_hist must be > 0.")

  k_min_star <- as.integer(k_min_star)
  k_max_star <- as.integer(k_max_star)
  if (k_min_star > k_max_star) stop("Need k_min_star <= k_max_star.")

  K <- k_max_star - k_min_star + 1L

  if (isTRUE(include_zero_bin)) {
    counts <- integer(K + 1L)
    names(counts) <- c("zero", as.character(k_min_star:k_max_star))
  } else {
    counts <- integer(K)
    names(counts) <- as.character(k_min_star:k_max_star)
  }

  for (val in u) {
    if (!is.finite(val) || val < 0) next

    if (isTRUE(include_zero_bin) && val == 0) {
      counts[["zero"]] <- counts[["zero"]] + 1L
    } else {
      kk <- floor(log(val, base = 2))
      kk <- as.integer(clamp(kk, k_min_star, k_max_star))
      counts[[as.character(kk)]] <- counts[[as.character(kk)]] + 1L
    }
  }

  noisy <- as.numeric(counts) + rlaplace(length(counts), scale = 1 / eps_hist)
  noisy <- pmax(noisy, 0)

  j_star <- which.max(noisy)
  name_star <- names(counts)[j_star]

  if (isTRUE(include_zero_bin) && identical(name_star, "zero")) {
    return(list(
      k_star = NA_integer_,
      m2_hat = 0,
      noisy = noisy,
      counts = counts
    ))
  }

  k_star <- as.integer(name_star)
  m2_hat <- 2^k_star

  list(
    k_star = k_star,
    m2_hat = m2_hat,
    noisy = noisy,
    counts = counts
  )
}

#' Private scale proxy via paired differences and block medians
#'
#' Internal scalar helper corresponding to the Algorithm 2 style routine in the
#' original research script.
#'
#' @param w Numeric vector.
#' @param eps_m2 Positive privacy epsilon.
#' @param M Number of blocks. If \code{NULL}, defaults to \code{floor(sqrt(n/2))}.
#' @param k_min_star,k_max_star Integer range of dyadic histogram bins.
#' @param include_zero_bin Logical; whether to include a separate zero bin.
#'
#' @return A list with \code{m2_hat}, \code{k_star}, \code{u}, \code{M}, and \code{hl}.
#' @keywords internal
dp_m2_alg2_scalar <- function(
    w,
    eps_m2,
    M = NULL,
    k_min_star = -20L,
    k_max_star = 40L,
    include_zero_bin = TRUE
) {
  w <- as.numeric(w)
  n <- length(w)

  if (n < 4) stop("Need length(w) >= 4 for Algorithm 2 (pairing + blocks).")
  if (!is.finite(eps_m2) || eps_m2 <= 0) stop("eps_m2 must be > 0.")

  m_pairs <- floor(n / 2)
  w1 <- w[seq(1, 2 * m_pairs, by = 2)]
  w2 <- w[seq(2, 2 * m_pairs, by = 2)]
  s <- (w1 - w2)^2 / 2

  if (is.null(M)) {
    M <- floor(sqrt(n / 2))
  }
  M <- as.integer(M)

  if (M < 1) stop("M must be >= 1.")
  if (M > length(s)) M <- length(s)

  block_size <- floor(length(s) / M)
  if (block_size < 1) {
    M <- length(s)
    block_size <- 1L
  }

  u <- numeric(M)
  for (b in seq_len(M)) {
    start <- (b - 1L) * block_size + 1L
    end <- if (b == M) length(s) else b * block_size
    u[b] <- stats::median(s[start:end])
  }

  hl <- dp_hist_argmax_laplace_log2bins(
    u = u,
    eps_hist = eps_m2,
    k_min_star = k_min_star,
    k_max_star = k_max_star,
    include_zero_bin = include_zero_bin
  )

  list(
    m2_hat = hl$m2_hat,
    k_star = hl$k_star,
    u = u,
    M = M,
    hl = hl
  )
}

#' Robustification parameter from private scale proxy
#'
#' @param m2_hat Nonnegative scale proxy.
#' @param eps_gd_ell Positive privacy epsilon for the GD part.
#' @param n_eff Effective sample size.
#'
#' @return A positive numeric scalar.
#' @keywords internal
tau_from_m2_scalar <- function(m2_hat, eps_gd_ell, n_eff) {
  m2_hat <- max(as.numeric(m2_hat), 0)

  if (!is.finite(eps_gd_ell) || eps_gd_ell <= 0) {
    stop("eps_gd_ell must be > 0.")
  }
  if (!is.numeric(n_eff) || length(n_eff) != 1 || n_eff < 2) {
    stop("n_eff must be >= 2.")
  }

  ln <- log(n_eff)
  denom <- sqrt((1 + ln) * ln)
  if (!is.finite(denom) || denom <= 0) denom <- 1

  sqrt(m2_hat) * sqrt((eps_gd_ell * n_eff) / denom)
}

#' DP Huber noisy gradient descent for a scalar mean
#'
#' @param w Numeric vector.
#' @param eps_gd,delta_gd Privacy parameters.
#' @param tau Positive Huber threshold.
#' @param T Number of iterations.
#' @param mu0 Initial point.
#' @param eta0 Base step size.
#' @param step_schedule Either \code{"fixed"} or \code{"1/sqrt(t)"}.
#'
#' @return A numeric scalar.
#' @keywords internal
dp_huber_noisy_gd_scalar <- function(
    w,
    eps_gd,
    delta_gd,
    tau,
    T,
    mu0 = 0,
    eta0 = 1,
    step_schedule = c("fixed", "1/sqrt(t)")
) {
  step_schedule <- match.arg(step_schedule)

  w <- as.numeric(w)
  n <- length(w)

  if (n < 2) stop("Need length(w) >= 2.")
  if (!is.finite(eps_gd) || eps_gd <= 0) stop("eps_gd must be > 0.")
  if (!is.finite(delta_gd) || delta_gd <= 0 || delta_gd >= 1) {
    stop("delta_gd must be in (0,1).")
  }
  if (!is.finite(tau) || tau <= 0) stop("tau must be > 0.")
  if (!is.finite(eta0) || eta0 <= 0) stop("eta0 must be > 0.")
  if (!is.numeric(T) || T < 1) stop("T must be >= 1.")

  T <- as.integer(T)
  mu <- as.numeric(mu0)

  eps_step <- eps_gd / T
  del_step <- delta_gd / T

  for (t in 0:(T - 1L)) {
    eta_t <- if (step_schedule == "fixed") eta0 else eta0 / sqrt(t + 1)

    r <- w - mu
    psi <- clamp(r, -tau, tau)
    g <- mean(psi)

    Delta_step <- (2 * eta_t * tau) / n
    sd_noise <- Delta_step * sqrt(2 * log(1.25 / del_step)) / eps_step

    mu <- mu + eta_t * g + stats::rnorm(1, mean = 0, sd = sd_noise)
  }

  mu
}
