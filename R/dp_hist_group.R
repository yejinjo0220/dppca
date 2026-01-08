
#' Group-wise DP 2D Histograms on PCA Scores
#'
#' @description
#' Constructs (optionally differentially private) 2D histograms on the first two
#' principal component scores, separately for each group. The function supports
#' differentially private PCA, a shared differentially private frame across
#' groups, and group-wise DP histogram mechanisms (Gaussian additive and sparse
#' Laplace thresholding).
#'
#' @param X A numeric matrix or data frame with \eqn{n} rows (observations) and
#'   \eqn{p} columns (variables and possibly a group label column). PCA is
#'   performed on the feature columns, and the first two PCs are used for
#'   binning.
#' @param G Group labels. Either
#'   \itemize{
#'     \item a vector of length \code{n} giving the group for each row of
#'           \code{X}, or
#'     \item a single character string giving the name of a column in \code{X}
#'           that contains the group labels.
#'   }
#'
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
#'   split among DP-PCA (optional), DP frame (\code{dp_frame()}), and group-wise
#'   DP histograms.
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
#'   axis. If \code{NULL}, it is chosen via \code{number_bins()} based on the
#'   pooled PCA scores.
#' @param m_y Optional integer specifying the number of bins along the second PC
#'   axis. If \code{NULL}, it is chosen via \code{number_bins()}.
#' @param bin_method Character string specifying the heuristic used by
#'   \code{number_bins()} when \code{m_x} and/or \code{m_y} are not supplied.
#'   One of \code{"J"} (Jing) or \code{"W"} (Wasserman).
#'
#' @param mechanism Character string indicating which histogram mechanism(s) to
#'   compute for each group:
#'   \itemize{
#'     \item \code{"none"}: (not implemented here; only DP mechanisms are
#'           computed).
#'     \item \code{"add"}: Gaussian additive-noise DP histogram per group.
#'     \item \code{"sparse"}: sparse Laplace-thresholded DP histogram per group.
#'     \item \code{"all"}: compute both \code{"add"} and \code{"sparse"}.
#'   }
#'
#' @param sampling Logical; if \code{TRUE}, additionally generate synthetic
#'   samples from each group's DP histogram(s) using \code{sample_from_hist()}
#'   with \code{k = n_g}, where \code{n_g} is the group size.
#'
#' @details
#' The function first separates the feature matrix and group labels. PCA
#' (possibly DP) is applied to the pooled data to obtain common principal
#' component axes. A single DP frame \code{(xlim, ylim)} is then constructed on
#' the pooled scores using \code{dp_frame()}, ensuring that all groups share the
#' same plotting region and bin grid.
#'
#' For each group, counts are computed on the common grid. Depending on
#' \code{mechanism}, Gaussian-perturbed and/or sparse Laplace-thresholded DP
#' histograms are constructed. Optionally, synthetic samples are drawn from
#' each group's DP histogram using \code{sample_from_hist()}.
#'
#' @return
#' A list with components:
#' \itemize{
#'   \item \code{pca}: A list containing
#'     \itemize{
#'       \item \code{X_pca}: \eqn{n \times 2} matrix of pooled PCA scores.
#'       \item \code{loadings}: loading matrix from \code{dp_pca()}.
#'       \item \code{eigvals}: selected eigenvalues.
#'       \item \code{dp}: logical flag indicating whether DP-PCA was used.
#'     }
#'   \item \code{frame}: A list with \code{xlim} and \code{ylim} giving the
#'     common DP plotting limits.
#'   \item \code{breaks}: A list with \code{x} and \code{y} break vectors.
#'   \item \code{eps}: A list with \code{total}, \code{pca}, \code{q},
#'     and \code{hist} components indicating how \code{eps_total} was split.
#'   \item \code{delta}: A list with \code{total}, \code{pca}, \code{q},
#'     and \code{hist} components indicating how \code{delta_total} was split.
#'   \item \code{groups}: A named list with one entry per group. For each group
#'     \code{g}, the entry contains:
#'     \itemize{
#'       \item \code{n}: group size \code{n_g}.
#'       \item \code{counts}: empirical bin counts.
#'       \item \code{add}: (if requested) a list with \code{p
dp_hist_group <- function(
    X,
    G, # group label

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
    bin_method = c("J", "W", "none"),

    mechanism = c("all", "none", "add", "sparse"),

    # Sampling
    sampling  = FALSE
) {
  bin_method <- match.arg(bin_method)
  # Manual bin: user supplies m_x, m_y
  if (bin_method == "none") {
    if (is.null(m_x) || is.null(m_y)) {
      stop("When bin_method = 'none', you must supply both m_x and m_y.", call. = FALSE)
    }
    if (!is.numeric(m_x) || !is.numeric(m_y) || any(!is.finite(c(m_x, m_y)))) {
      stop("'m_x' and 'm_y' must be finite numeric.", call. = FALSE)
    }
    m_x <- as.integer(m_x)
    m_y <- as.integer(m_y)
    if (m_x < 1 || m_y < 1) stop("Both 'm_x' and 'm_y' must be at least 1.", call. = FALSE)
  }

  mechanism  <- match.arg(mechanism)

  X <- as.data.frame(X)

  if (is.character(G) && length(G) == 1L) {
    if (!G %in% colnames(X)) {
      stop("Column '", G, "' not found in X.")
    }
    group_vec <- X[[G]]
    X_feat    <- X[, setdiff(colnames(X), G), drop = FALSE]
  } else {
    group_vec <- G
    X_feat    <- X
  }

  X_mat <- as.matrix(X_feat)
  n <- nrow(X_mat)
  p <- ncol(X_mat)
  if (length(group_vec) != n) stop("length(G) must equal nrow(X).")
  if (n <= 1) stop("X must have at least two rows.")
  if (p < 2) stop("X must have at least two columns for PCA.")

  groups <- unique(group_vec)

  if (!is.finite(eps_total) || eps_total <= 0) {
    stop("'eps_total' must be a positive finite number.")
  }
  if (!is.finite(delta_total) || delta_total <= 0) {
    stop("'delta_total' must be a positive finite number.")
  }

  # Divide privacy budget
  if (dp_pca_flag) {
    # Divide into 3: (PCA, q, hist)
    if (is.null(eps_ratio)) {
      eps_ratio <- c(1, 1, 1)
    } else {
      if (!is.numeric(eps_ratio) || length(eps_ratio) != 3)
        stop("When dp_pca_flag = TRUE, 'eps_ratio' must be length 3 (PCA, q, hist).")
      if (any(eps_ratio < 0)) stop("'eps_ratio' must be nonnegative.")
    }
    if (is.null(delta_ratio)) {
      delta_ratio <- c(1, 1, 1)
    } else {
      if (!is.numeric(delta_ratio) || length(delta_ratio) != 3)
        stop("When dp_pca_flag = TRUE, 'delta_ratio' must be length 3 (PCA, q, hist).")
      if (any(delta_ratio < 0)) stop("'delta_ratio' must be nonnegative.")
    }

    eps_ratio   <- eps_ratio   / sum(eps_ratio)
    delta_ratio <- delta_ratio / sum(delta_ratio)

    eps_pca  <- eps_total   * eps_ratio[1]
    eps_q    <- eps_total   * eps_ratio[2]
    eps_hist <- eps_total   * eps_ratio[3]

    delta_pca  <- delta_total * delta_ratio[1]
    delta_q    <- delta_total * delta_ratio[2]
    delta_hist <- delta_total * delta_ratio[3]

  } else {
    # Divide into 2: (q, hist)
    if (is.null(eps_ratio)) {
      eps_ratio <- c(1, 1)
    } else {
      if (!is.numeric(eps_ratio) || length(eps_ratio) != 2)
        stop("When dp_pca_flag = FALSE, 'eps_ratio' must be length 2 (q, hist).")
      if (any(eps_ratio < 0)) stop("'eps_ratio' must be nonnegative.")
    }
    if (is.null(delta_ratio)) {
      delta_ratio <- c(1, 1)
    } else {
      if (!is.numeric(delta_ratio) || length(delta_ratio) != 2)
        stop("When dp_pca_flag = FALSE, 'delta_ratio' must be length 2 (q, hist).")
      if (any(delta_ratio < 0)) stop("'delta_ratio' must be nonnegative.")
    }

    eps_ratio   <- eps_ratio   / sum(eps_ratio)
    delta_ratio <- delta_ratio / sum(delta_ratio)

    eps_pca  <- NULL
    delta_pca <- NULL

    eps_q    <- eps_total   * eps_ratio[1]
    eps_hist <- eps_total   * eps_ratio[2]

    delta_q    <- delta_total * delta_ratio[1]
    delta_hist <- delta_total * delta_ratio[2]
  }

  # DP-PCA
  pca_res <- dp_pca(
    X         = X_mat,
    center    = center,
    scale.    = scale.,
    axes      = axes,
    dp        = dp_pca_flag,
    eps       = eps_pca,
    delta     = delta_pca,
    cpp.option = cpp.option
  )
  X_pca <- as.matrix(pca_res$X_pca) # n x 2

  # DP frame
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

  # bin count
  if (bin_method == "none") {
    if (is.null(m_x) || is.null(m_y)) {
      stop("When bin_method = 'none', you must supply both m_x and m_y.", call. = FALSE)
    }
    if (!is.numeric(m_x) || !is.numeric(m_y) || any(!is.finite(c(m_x, m_y)))) {
      stop("'m_x' and 'm_y' must be finite numeric.", call. = FALSE)
    }
    m_x <- as.integer(m_x)
    m_y <- as.integer(m_y)
    if (m_x < 1 || m_y < 1) stop("Both 'm_x' and 'm_y' must be at least 1.", call. = FALSE)
  }

  # bin count
  if (bin_method != "none") {
    if (is.null(m_x) && is.null(m_y)) {
      m_default <- number_bins(X_pca, method = bin_method)
      m_x <- m_default
      m_y <- m_default
    } else {
      if (is.null(m_x)) m_x <- number_bins(X_pca, method = bin_method)
      if (is.null(m_y)) m_y <- number_bins(X_pca, method = bin_method)
    }
  }

  if (m_x < 1 || m_y < 1) stop("Both 'm_x' and 'm_y' must be at least 1.")

  x_breaks <- seq(xlim[1], xlim[2], length.out = m_x + 1)
  y_breaks <- seq(ylim[1], ylim[2], length.out = m_y + 1)
  m   <- m_x * m_y

  # Histogram with group
  result_groups <- list()

  for (g in unique(group_vec)) {
    idx_g <- which(group_vec == g)
    Xg    <- X_pca[idx_g, , drop = FALSE]
    n_g   <- nrow(Xg)
    if (n_g == 0) next

    # bin index
    bx <- cut(Xg[,1], breaks = x_breaks, include.lowest = TRUE, labels = FALSE)
    by <- cut(Xg[,2], breaks = y_breaks, include.lowest = TRUE, labels = FALSE)

    if (anyNA(bx)) {
      idx <- is.na(bx)
      bx[idx] <- pmax(1, pmin(m_x, findInterval(Xg[idx, 1], x_breaks, all.inside = TRUE)))
    }
    if (anyNA(by)) {
      idx <- is.na(by)
      by[idx] <- pmax(1, pmin(m_y, findInterval(Xg[idx, 2], y_breaks, all.inside = TRUE)))
    }

    bidx   <- (by - 1) * m_x + bx
    counts <- as.numeric(table(factor(bidx, levels = 1:m)))

    group_res <- list(
      n      = n_g,
      counts = counts
    )

    # Add
    if (mechanism %in% c("all", "add")) {
      sigma <- sqrt(2) * sqrt(2 * log(1.25 / delta_hist)) / eps_hist
      c_tilde <- pmax(counts + stats::rnorm(m, mean = 0, sd = sigma), 0)
      if (sum(c_tilde) <= 0) {
        stop(sprintf("Group '%s': All privatized counts are zero after Gaussian noise.", g))
      }
      q_add <- c_tilde / sum(c_tilde)

      coord_add <- do.call(
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
      coord_add$prob <- q_add

      if (sampling) {
        samp_add <- sample_from_hist(coord_add, k = n_g)
        group_res$add <- list(prob = q_add, coord = coord_add, sample = samp_add)
      } else {
        group_res$add <- list(prob = q_add, coord = coord_add)
      }
    }

    # Sparse
    if (mechanism %in% c("all", "sparse")) {
      q_sparse <- numeric(m)
      scale_lap <- 2 / eps_hist
      thr <- (2 * log(2 / delta_hist)) / eps_hist + 1

      for (k in 1:m) {
        ck <- counts[k]
        if (ck == 0) {
          q_sparse[k] <- 0
        } else {
          Zk <- VGAM::rlaplace(1, location = 0, scale = scale_lap)
          c_tilde <- max(ck + Zk, 0)
          if (c_tilde >= thr) {
            q_sparse[k] <- c_tilde
          } else {
            q_sparse[k] <- 0
          }
        }
      }

      if (sum(q_sparse) > 0) {
        q_sparse <- q_sparse / sum(q_sparse)
      } else {
      }

      coord_sp <- do.call(
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
      coord_sp$prob <- q_sparse

      if (sampling) {
        if (sum(q_sparse) > 0) {
          samp_sp <- sample_from_hist(coord_sp, k = n_g)
          group_res$sparse <- list(prob = q_sparse, coord = coord_sp, sample = samp_sp)
        } else {
          group_res$sparse <- list(prob = q_sparse, coord = coord_sp, sample = NULL)
        }
      } else {
        group_res$sparse <- list(prob = q_sparse, coord = coord_sp)
      }
    }

    result_groups[[as.character(g)]] <- group_res
  }

  # result
  result <- list(
    pca = list(
      X_pca   = X_pca,
      loadings = pca_res$loadings,
      eigvals  = pca_res$eigvals,
      dp       = pca_res$dp
    ),
    frame = list(
      xlim = xlim,
      ylim = ylim
    ),
    breaks = list(x = x_breaks, y = y_breaks),
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
    ),
    groups = result_groups,
    mechanism = mechanism
  )

  result
}
