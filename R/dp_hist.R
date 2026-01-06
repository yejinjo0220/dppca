
#' Differentially Private 2D Histogram on PCA Scores
#'
#' @description
#' Constructs a (optionally differentially private) 2D histogram on the first
#' two principal component scores of the input data. The function supports
#' differentially private PCA, differentially private frame construction, and
#' two types of DP histogram mechanisms (Gaussian additive noise and sparse
#' Laplace thresholding), as well as the non-DP empirical histogram.
#'
#' @param X A numeric matrix or an object coercible to a matrix with
#'   \eqn{n} rows (observations) and \eqn{p} columns (variables). PCA is
#'   performed on \code{X} and the first two PCs are used for binning.
#' @param center Logical; should the variables be centered before PCA?
#'   Passed to \code{scale()}.
#' @param scale. Logical; should the variables be scaled to unit variance?
#'   Passed to \code{scale()}.
#' @param dp_pca_flag Logical; if \code{TRUE}, apply differentially private PCA
#'   using \code{dp_pca()} with part of the total privacy budget
#'   \code{eps_total}, \code{delta_total}. If \code{FALSE}, use non-DP PCA.
#' @param cpp.option Logical; if \code{TRUE}, allow C++ implementations inside
#'   DP mechanisms (e.g., \code{mech_tau_sph()}) when available.
#' @param axes Integer vector indicating which principal components to return
#'   (e.g., \code{c(1, 2)} for the first two PCs). The maximum index in
#'   \code{axes} must not exceed the number of variables \code{p}.
#'
#' @param eps_total Total privacy budget \eqn{\varepsilon_{\text{total}}} to be
#'   split among DP-PCA (optional), DP frame (\code{dp_frame()}), and DP
#'   histogram mechanisms.
#' @param delta_total Total privacy parameter \eqn{\delta_{\text{total}}} to be
#'   split among the same components.
#' @param eps_ratio Optional numeric vector specifying the relative allocation
#'   of \code{eps_total}. If \code{dp_pca_flag = TRUE}, this must have length 3
#'   and corresponds to \code{(PCA, q, hist)}. If \code{dp_pca_flag = FALSE},
#'   this must have length 2 and corresponds to \code{(q, hist)}. The vector is
#'   normalized to sum to 1. If \code{NULL}, equal weights are used.
#' @param delta_ratio Optional numeric vector specifying the relative allocation
#'   of \code{delta_total}, with the same conventions as \code{eps_ratio}. The
#'   vector is normalized to sum to 1. If \code{NULL}, equal weights are used.
#'
#' @param inflate Non-negative numeric value controlling how much to enlarge the
#'   DP frame relative to the DP min--max interval obtained from
#'   \code{dp_frame()}. A value of \code{inflate = 0.10} enlarges the interval
#'   by 10\%.
#' @param q_frame Optional symmetric quantile level in \eqn{(0,1)} passed to
#'   \code{dp_frame()}. If \code{NULL}, extreme order statistics
#'   \eqn{1/n} and \eqn{(n-1)/n} are used.
#'
#' @param m_x Optional integer specifying the number of bins along the first PC
#'   axis. If \code{NULL}, it is chosen via \code{number_bins()}.
#' @param m_y Optional integer specifying the number of bins along the second PC
#'   axis. If \code{NULL}, it is chosen via \code{number_bins()}.
#' @param bin_method Character string specifying the heuristic used by
#'   \code{number_bins()} when \code{m_x} and/or \code{m_y} are not supplied.
#'   One of \code{"J"} (Jing) or \code{"W"} (Wasserman).
#'
#' @param mechanism Character string indicating which histogram mechanism(s) to
#'   compute:
#'   \itemize{
#'     \item \code{"none"}: non-DP empirical histogram on PCA scores.
#'     \item \code{"add"}: Gaussian additive noise on bin counts with
#'           subsequent normalization.
#'     \item \code{"sparse"}: Laplace noise with thresholding (sparse DP
#'           histogram).
#'     \item \code{"all"}: compute all three versions.
#'   }
#'
#' @param sampling Logical; if \code{TRUE}, additionally generate synthetic
#'   samples from each (available) histogram version using
#'   \code{sample_from_hist()} with \code{k = n}.
#'
#' @details
#' The function proceeds in several steps:
#' \enumerate{
#'   \item \strong{Privacy budget split}: The total budget
#'         \code{(eps_total, delta_total)} is split into components for DP-PCA
#'         (optional), DP frame, and DP histogram according to
#'         \code{eps_ratio} and \code{delta_ratio}.
#'   \item \strong{(DP) PCA}: The data \code{X} are projected to the first two
#'         principal components using \code{dp_pca()} with or without privacy,
#'         depending on \code{dp_pca_flag}.
#'   \item \strong{DP frame}: A differentially private square plotting frame
#'         \code{(xlim, ylim)} is constructed on the 2D scores using
#'         \code{dp_frame()}.
#'   \item \strong{Bin selection}: If \code{m_x} and/or \code{m_y} are missing,
#'         a recommended number of bins per axis is computed by
#'         \code{number_bins()}.
#'   \item \strong{Empirical histogram}: Counts and empirical probabilities
#'         \code{p_hat} are computed on the regular grid induced by
#'         \code{x_breaks} and \code{y_breaks}.
#'   \item \strong{DP histogram mechanisms}:
#'         \itemize{
#'           \item \code{"none"}: returns the empirical probabilities.
#'           \item \code{"add"}: adds Gaussian noise to the bin counts, clips to
#'                 nonnegative values, and renormalizes.
#'           \item \code{"sparse"}: adds Laplace noise to the empirical bin
#'                 probabilities and applies a data-dependent threshold to obtain
#'                 a sparse DP histogram.
#'         }
#'   \item \strong{Optional sampling}: If \code{sampling = TRUE}, each available
#'         histogram version is converted into synthetic samples via
#'         \code{sample_from_hist()}.
#' }
#'
#' @return A list with components:
#' \itemize{
#'   \item \code{pca}: A list containing
#'     \itemize{
#'       \item \code{X_pca}: \eqn{n \times 2} matrix of PCA scores.
#'       \item \code{loadings}: loading matrix from \code{dp_pca()}.
#'       \item \code{eigvals}: selected eigenvalues.
#'       \item \code{dp}: logical flag indicating whether DP-PCA was used.
#'     }
#'   \item \code{frame}: A list with \code{xlim} and \code{ylim} giving the DP
#'     plotting limits.
#'   \item \code{breaks}: A list with \code{x} and \code{y} break vectors.
#'   \item \code{counts}: Vector of empirical bin counts.
#'   \item \code{prob}: Vector of empirical bin probabilities \code{p_hat}.
#'   \item \code{mechanism}: The mechanism argument used.
#'   \item \code{eps}: A list with \code{total}, \code{pca}, \code{q},
#'     and \code{hist} components indicating how \code{eps_total} was split.
#'   \item \code{delta}: A list with \code{total}, \code{pca}, \code{q},
#'     and \code{hist} components indicating how \code{delta_total} was split.
#'   \item \code{none}: If requested (via \code{mechanism}), a list with fields
#'     \code{prob}, \code{coord}, and optionally \code{sample} for the empirical
#'     histogram.
#'   \item \code{add}: If requested, a list with fields \code{prob}, \code{coord},
#'     and optionally \code{sample} for the Gaussian-perturbed DP histogram.
#'   \item \code{sparse}: If requested, a list with fields \code{prob},
#'     \code{coord}, and optionally \code{sample} for the sparse DP histogram.
#' }
#'
#' @seealso
#'   \code{\link{dp_pca}},
#'   \code{\link{dp_frame}},
#'   \code{\link{number_bins}},
#'   \code{\link{sample_from_hist}},
#'   \code{\link{dp_score_plot}},
#'   \code{\link{dp_score_plot_group}}
#'
#' @export
dp_hist <- function(
    X,
    # PCA option
    center = TRUE,
    scale. = FALSE,
    dp_pca_flag = FALSE,
    cpp.option = TRUE,
    axes = c(1, 2),

    # Privacy budget
    eps_total,
    delta_total,
    eps_ratio   = NULL,
    delta_ratio = NULL,

    # Frame option
    inflate = 0.10,
    q_frame = NULL,

    # bin count
    m_x = NULL,
    m_y = NULL,
    bin_method = c("J", "W"),

    mechanism = c("all", "none", "add", "sparse"),

    # Sampling
    sampling  = FALSE
) {
  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)
  if (n <= 1) stop("X must have at least two rows.")
  if (p < 2) stop("X must have at least two columns for PCA.")

  if (!is.finite(eps_total) || eps_total <= 0) {
    stop("'eps_total' must be a positive finite number.")
  }
  if (!is.finite(delta_total) || delta_total <= 0) {
    stop("'delta_total' must be a positive finite number.")
  }

  bin_method <- match.arg(bin_method)
  mechanism  <- match.arg(mechanism)

  # Divide privacy budget ------------------------------------------------------
  if (dp_pca_flag) {
    # Divide into 3: (PCA, q, hist)
    if (is.null(eps_ratio)) {
      eps_ratio <- c(1, 1, 1)
    } else {
      if (!is.numeric(eps_ratio) || length(eps_ratio) != 3)
        stop("When dp_pca_flag = TRUE, 'eps_ratio' must be numeric of length 3 (PCA, q, hist).")
      if (any(eps_ratio < 0)) stop("'eps_ratio' must be nonnegative.")
    }
    if (is.null(delta_ratio)) {
      delta_ratio <- c(1, 1, 1)
    } else {
      if (!is.numeric(delta_ratio) || length(delta_ratio) != 3)
        stop("When dp_pca_flag = TRUE, 'delta_ratio' must be numeric of length 3 (PCA, q, hist).")
      if (any(delta_ratio < 0)) stop("'delta_ratio' must be nonnegative.")
    }

    eps_ratio <- eps_ratio / sum(eps_ratio)
    delta_ratio <- delta_ratio / sum(delta_ratio)

    eps_pca <- eps_total * eps_ratio[1]
    eps_q <- eps_total * eps_ratio[2]
    eps_hist <- eps_total * eps_ratio[3]

    delta_pca <- delta_total * delta_ratio[1]
    delta_q <- delta_total * delta_ratio[2]
    delta_hist <- delta_total * delta_ratio[3]

  } else {
    # Divide into 2: (q, hist)
    if (is.null(eps_ratio)) {
      eps_ratio <- c(1, 1)
    } else {
      if (!is.numeric(eps_ratio) || length(eps_ratio) != 2)
        stop("When dp_pca_flag = FALSE, 'eps_ratio' must be numeric of length 2 (q, hist).")
      if (any(eps_ratio < 0)) stop("'eps_ratio' must be nonnegative.")
    }
    if (is.null(delta_ratio)) {
      delta_ratio <- c(1, 1)
    } else {
      if (!is.numeric(delta_ratio) || length(delta_ratio) != 2)
        stop("When dp_pca_flag = FALSE, 'delta_ratio' must be numeric of length 2 (q, hist).")
      if (any(delta_ratio < 0)) stop("'delta_ratio' must be nonnegative.")
    }
    eps_ratio <- eps_ratio / sum(eps_ratio)
    delta_ratio <- delta_ratio / sum(delta_ratio)

    eps_pca <- NULL
    delta_pca <- NULL

    eps_q <- eps_total * eps_ratio[1]
    eps_hist <- eps_total * eps_ratio[2]

    delta_q <- delta_total * delta_ratio[1]
    delta_hist <- delta_total * delta_ratio[2]
  }


  # DP-PCA ---------------------------------------------------
  pca_res <- dp_pca(
    X,
    center = center,
    scale. = scale.,
    axes   = axes,
    dp     = dp_pca_flag,
    eps    = eps_pca,
    delta  = delta_pca,
    cpp.option = cpp.option
  )
  X_pca <- as.matrix(pca_res$X_pca) # n x 2
  colnames(X_pca) <- c("PC1", "PC2")

  # DP-Frame ------------------------------------------------
  frame_dp <- dp_frame(
    X       = X_pca,
    eps_q   = eps_q,
    delta_q = delta_q,
    inflate = inflate,
    q       = q_frame
  )
  xlim <- frame_dp$xlim
  ylim <- frame_dp$ylim

  if (length(xlim) != 2L || length(ylim) != 2L) {
    stop("dp_frame must return numeric xlim, ylim of length 2.")
  }

  # number of bins ------------------------------------------------------
  if (is.null(m_x) && is.null(m_y)) {
    m_default <- number_bins(X_pca, method = bin_method)
    m_x <- m_default
    m_y <- m_default
  } else {
    if (is.null(m_x)) {
      m_x <- number_bins(X_pca, method = bin_method)
    }
    if (is.null(m_y)) {
      m_y <- number_bins(X_pca, method = bin_method)
    }
  }
  if (m_x < 1 || m_y < 1) stop("Both 'm_x' and 'm_y' must be at least 1.")

  # uniform breaks ---------------------------------------------------------
  x_breaks <- seq(xlim[1], xlim[2], length.out = m_x + 1)
  y_breaks <- seq(ylim[1], ylim[2], length.out = m_y + 1)

  m_x <- length(x_breaks) - 1
  m_y <- length(y_breaks) - 1
  m   <- m_x * m_y

  # calculate count --------------------------------------------------------------
  bx <- cut(X_pca[, 1], breaks = x_breaks, include.lowest = TRUE, labels = FALSE)
  by <- cut(X_pca[, 2], breaks = y_breaks, include.lowest = TRUE, labels = FALSE)

  if (anyNA(bx)) {
    idx <- is.na(bx)
    bx[idx] <- pmax(1, pmin(m_x, findInterval(X_pca[idx, 1], x_breaks, all.inside = TRUE)))
  }
  if (anyNA(by)) {
    idx <- is.na(by)
    by[idx] <- pmax(1, pmin(m_y, findInterval(X_pca[idx, 2], y_breaks, all.inside = TRUE)))
  }

  bidx   <- (by - 1) * m_x + bx
  counts <- as.numeric(table(factor(bidx, levels = 1:m)))
  p_hat  <- counts / n

  # DP-Histogram ---------------------------------------------------------------
  # 1. None
  q_none <- p_hat

  # 2 Add
  q_add <- NULL
  if (mechanism %in% c("all", "add")) {

    sigma <- sqrt(2) * sqrt(2 * log(1.25 / delta_hist)) / eps_hist
    c_tilde <- pmax(counts + stats::rnorm(m, 0, sigma), 0)

    if (sum(c_tilde) <= 0) {
      stop("All privatized bin counts are zero after Gaussian noise and clipping. ",
           "Try using larger 'eps_total' or fewer bins.")
    }
    q_add <- c_tilde / sum(c_tilde)
  }

  # 3. Sparse
  q_sparse <- NULL
  if (mechanism %in% c("all", "sparse")) {
    q_sparse <- numeric(m)

    scale_lap <- 2 / (eps_hist * n)
    thr <- (2 * log(2 / delta_hist)) / (eps_hist * n) + (1 / n)

    for (k in 1:m) {
      if (p_hat[k] == 0) {
        q_sparse[k] <- 0
      } else {
        Zk <- VGAM::rlaplace(1, location = 0, scale = scale_lap)
        qk <- max(p_hat[k] + Zk, 0)
        q_sparse[k] <- if (qk < thr) 0 else qk
      }
    }

    if (sum(q_sparse) > 0) {
      q_sparse <- q_sparse / sum(q_sparse)
    } else {
    }
  }


  # result -------------------------------------------------------
  base_coord <- do.call(
    rbind,
    lapply(1:m_y, function(j) {
      data.frame(
        xmin = x_breaks[-length(x_breaks)],
        xmax = x_breaks[-1],
        ymin = y_breaks[j],
        ymax = y_breaks[j + 1]
      )
    })
  )

  result <- list(
    # PCA
    pca = list(
      X_pca   = X_pca,
      loadings = pca_res$loadings,
      eigvals  = pca_res$eigvals,
      dp       = pca_res$dp
    ),
    # frame
    frame = list(
      xlim = xlim,
      ylim = ylim
    ),
    # bin / count
    breaks = list(x = x_breaks, y = y_breaks),
    counts = counts,
    prob   = p_hat,
    mechanism = mechanism,
    eps = list(
      total = eps_total,
      pca   = eps_pca,
      q     = eps_q,
      hist  = eps_hist
    ),
    delta = list(
      total = delta_total,
      pca   = delta_pca,
      q     = delta_q,
      hist  = delta_hist
    )
  )

  # none
  if (mechanism %in% c("all", "none")) {
    coord_none <- base_coord
    coord_none$prob <- q_none
    coord_none$mechanism <- "none"

    if (sampling) {
      sample_none <- sample_from_hist(coord_none, k = n)
      result$none <- list(prob = q_none, coord = coord_none, sample = sample_none)
    } else {
      result$none <- list(prob = q_none, coord = coord_none)
    }
  }

  # add
  if (mechanism %in% c("all", "add")) {
    coord_add <- base_coord
    coord_add$prob <- q_add
    coord_add$mechanism <- "add"

    if (sampling) {
      sample_add <- sample_from_hist(coord_add, k = n)
      result$add <- list(prob = q_add, coord = coord_add, sample = sample_add)
    } else {
      result$add <- list(prob = q_add, coord = coord_add)
    }
  }

  # sparse
  if (mechanism %in% c("all", "sparse")) {
    coord_sparse <- base_coord
    coord_sparse$prob <- q_sparse
    coord_sparse$mechanism <- "sparse"

    if (sampling && sum(q_sparse) > 0) {
      sample_sparse <- sample_from_hist(coord_sparse, k = n)
      result$sparse <- list(prob = q_sparse, coord = coord_sparse, sample = sample_sparse)
    } else {
      result$sparse <- list(prob = q_sparse, coord = coord_sparse)
    }
  }

  result
}

