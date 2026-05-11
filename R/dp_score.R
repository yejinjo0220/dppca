# ============================================================
# dp_score_final.R
# User-facing functions for DP score histograms
# ============================================================

.split_score_privacy_budget <- function(eps_total, delta_total, g_dppca) {
  if (isTRUE(g_dppca)) {
    list(
      eps_dir = eps_total / 3,
      eps_q = eps_total / 3,
      eps_hist = eps_total / 3,
      delta_dir = delta_total / 3,
      delta_q = delta_total / 3,
      delta_hist = delta_total / 3
    )
  } else {
    list(
      eps_dir = NULL,
      eps_q = eps_total / 2,
      eps_hist = eps_total / 2,
      delta_dir = NULL,
      delta_q = delta_total / 2,
      delta_hist = delta_total / 2
    )
  }
}

#' Differentially private score histograms
#'
#' @description
#' Computes two-dimensional PCA scores from the input data and constructs an
#' empirical histogram and selected differentially private histogram versions on
#' the resulting score space.
#'
#' @param X Numeric matrix or an object coercible to a matrix with rows
#'   corresponding to observations and columns corresponding to variables.
#' @param center Logical; whether the variables should be centered before PCA.
#' @param standardize Logical; whether the variables should be scaled to unit
#'   variance before PCA.
#' @param g_dppca Logical; whether to use a differentially private PCA direction
#'   matrix. If \code{FALSE}, the usual sample PCA directions are used.
#' @param cpp.option Logical passed to \code{mech_tau_sph()} when private
#'   directions are computed.
#' @param axes Integer vector of length 2 indicating which principal components
#'   to use. Default is \code{c(1, 2)}.
#' @param eps_total Total privacy budget.
#' @param delta_total Total privacy parameter.
#' @param method Character vector specifying which DP histogram method(s) to
#'   compute. Options are \code{"add"} for the Gaussian additive DP histogram and
#'   \code{"sparse"} for the sparse Laplace-thresholded DP histogram.
#' @param frame Optional user-specified frame.
#' @param m_x Optional number of bins along the x-axis.
#' @param m_y Optional number of bins along the y-axis.
#' @param bin_method Character string specifying the heuristic used by
#'   \code{recommend_bins()} when \code{m_x} and/or \code{m_y} are not supplied.
#'   One of \code{"WZ"}, \code{"Lei"}, or \code{"none"}.
#'
#' @return A list with components:
#' \itemize{
#'   \item \code{score}: \eqn{n \times 2} score matrix.
#'   \item \code{frame}: a list with \code{xlim} and \code{ylim}.
#'   \item \code{none}: data frame for the empirical histogram.
#'   \item \code{add}: data frame for the Gaussian additive DP histogram, or
#'         \code{NULL} if not requested.
#'   \item \code{sparse}: data frame for the sparse DP histogram, or
#'         \code{NULL} if not requested.
#' }
#'
#' @importFrom VGAM rlaplace
#' @export
dp_score <- function(
    X,
    center = TRUE,
    standardize = FALSE,
    g_dppca = FALSE,
    cpp.option = FALSE,
    axes = c(1, 2),
    eps_total,
    delta_total,
    method = c("add", "sparse"),
    frame = NULL,
    m_x = NULL,
    m_y = NULL,
    bin_method = c("WZ", "Lei", "none")
) {
  X <- as.matrix(X)

  if (!is.numeric(X)) stop("X must be numeric or coercible to a numeric matrix.")
  if (nrow(X) < 2) stop("X must have at least two rows.")
  if (ncol(X) < 2) stop("X must have at least two columns.")
  if (length(axes) != 2) stop("'axes' must be an integer vector of length 2.")
  if (any(!is.finite(axes)) || any(axes <= 0) || any(axes != as.integer(axes))) {
    stop("'axes' must contain positive integers.")
  }
  if (max(axes) > ncol(X)) stop("The requested component index in 'axes' exceeds ncol(X).")
  if (!is.finite(eps_total) || eps_total <= 0) stop("'eps_total' must be a positive finite number.")
  if (!is.finite(delta_total) || delta_total <= 0 || delta_total >= 1) {
    stop("'delta_total' must be in (0, 1).")
  }

  axes <- as.integer(axes)
  k_max <- max(axes)
  method <- match.arg(method, choices = c("add", "sparse"), several.ok = TRUE)
  bin_method <- match.arg(bin_method)
  n <- nrow(X)

  budget <- .split_score_privacy_budget(
    eps_total = eps_total,
    delta_total = delta_total,
    g_dppca = g_dppca
  )

  X_proc <- prep_matrix_for_pca(
    X = X,
    center = center,
    standardize = standardize
  )

  V_all <- dp_pc_dir(
    X = X,
    k = k_max,
    center = center,
    standardize = standardize,
    g_dppca = g_dppca,
    eps_dir = budget$eps_dir,
    delta_dir = budget$delta_dir,
    cpp.option = cpp.option
  )

  V <- V_all[, axes, drop = FALSE]
  X_score <- as.matrix(X_proc %*% V)
  colnames(X_score) <- paste0("PC", axes)

  frame_out <- dp_frame(
    X       = X_score,
    eps_q   = budget$eps_q,
    delta_q = budget$delta_q,
    inflate = 0.20,
    q       = NULL,
    frame   = frame
  )

  xlim <- frame_out$xlim
  ylim <- frame_out$ylim

  if (bin_method == "none") {
    if (is.null(m_x) || is.null(m_y)) {
      stop("When bin_method = 'none', both m_x and m_y must be supplied.")
    }
  } else {
    if (is.null(m_x)) m_x <- recommend_bins(X_score, method = bin_method)
    if (is.null(m_y)) m_y <- recommend_bins(X_score, method = bin_method)
  }

  if (!is.numeric(m_x) || length(m_x) != 1 || m_x < 1) stop("'m_x' must be a positive integer.")
  if (!is.numeric(m_y) || length(m_y) != 1 || m_y < 1) stop("'m_y' must be a positive integer.")

  m_x <- as.integer(m_x)
  m_y <- as.integer(m_y)
  m   <- m_x * m_y

  x_breaks <- seq(xlim[1], xlim[2], length.out = m_x + 1)
  y_breaks <- seq(ylim[1], ylim[2], length.out = m_y + 1)

  bx <- cut(X_score[, 1], breaks = x_breaks, include.lowest = TRUE, labels = FALSE)
  by <- cut(X_score[, 2], breaks = y_breaks, include.lowest = TRUE, labels = FALSE)

  if (anyNA(bx)) {
    idx <- which(is.na(bx))
    bx[idx] <- findInterval(X_score[idx, 1], x_breaks, all.inside = TRUE)
  }
  if (anyNA(by)) {
    idx <- which(is.na(by))
    by[idx] <- findInterval(X_score[idx, 2], y_breaks, all.inside = TRUE)
  }

  bidx   <- (by - 1) * m_x + bx
  counts <- as.numeric(table(factor(bidx, levels = 1:m)))
  p_hat  <- counts / n

  base_coord <- do.call(
    rbind,
    lapply(seq_len(m_y), function(j) {
      data.frame(
        xmin = x_breaks[-length(x_breaks)],
        xmax = x_breaks[-1],
        ymin = y_breaks[j],
        ymax = y_breaks[j + 1]
      )
    })
  )

  hist_none <- base_coord
  hist_none$prob <- p_hat

  hist_add <- NULL
  hist_sparse <- NULL

  eps_hist_method <- budget$eps_hist / length(method)
  delta_hist_method <- budget$delta_hist / length(method)

  if ("add" %in% method) {
    sigma <- sqrt(2) * sqrt(2 * log(1.25 / delta_hist_method)) / eps_hist_method
    c_tilde <- pmax(counts + stats::rnorm(m, mean = 0, sd = sigma), 0)

    if (sum(c_tilde) <= 0) {
      stop("All privatized bin counts are zero after Gaussian noise and clipping. Try larger eps_total or fewer bins.")
    }

    hist_add <- base_coord
    hist_add$prob <- c_tilde / sum(c_tilde)
  }

  if ("sparse" %in% method) {
    q_sparse <- numeric(m)
    scale_lap <- 2 / (eps_hist_method * n)
    thr <- (2 * log(2 / delta_hist_method)) / (eps_hist_method * n) + (1 / n)

    for (kk in seq_len(m)) {
      if (p_hat[kk] == 0) {
        q_sparse[kk] <- 0
      } else {
        z_k <- VGAM::rlaplace(1, location = 0, scale = scale_lap)
        q_k <- max(p_hat[kk] + z_k, 0)
        q_sparse[kk] <- if (q_k < thr) 0 else q_k
      }
    }

    if (sum(q_sparse) > 0) q_sparse <- q_sparse / sum(q_sparse)

    hist_sparse <- base_coord
    hist_sparse$prob <- q_sparse
  }

  list(
    score  = X_score,
    frame  = list(xlim = xlim, ylim = ylim),
    none   = hist_none,
    add    = hist_add,
    sparse = hist_sparse,
    method = method
  )
}

#' Plot differentially private score histograms
#'
#' @description
#' Computes PCA score histograms through \code{dp_score()} and returns both the
#' calculation results and the corresponding plots.
#'
#' @param X Numeric matrix or an object coercible to a matrix.
#' @param center Logical; whether the variables should be centered before PCA.
#' @param standardize Logical; whether the variables should be scaled before PCA.
#' @param g_dppca Logical; whether to use a DP PCA direction matrix.
#' @param cpp.option Logical passed to \code{mech_tau_sph()} when private
#'   directions are computed.
#' @param axes Integer vector of length 2 specifying the principal components
#'   used for score construction.
#' @param eps_total Total privacy budget.
#' @param delta_total Total privacy parameter.
#' @param method Character vector specifying which DP histogram method(s) to
#'   plot. Options are \code{"add"} and \code{"sparse"}.
#' @param frame Optional user-specified frame.
#' @param m_x Optional number of bins along the x-axis.
#' @param m_y Optional number of bins along the y-axis.
#' @param bin_method Character string specifying the bin recommendation rule.
#'
#' @return A list with components:
#' \itemize{
#'   \item \code{score}: the output of \code{dp_score()}.
#'   \item \code{plot}: a list containing the combined plot and each
#'         individual panel.
#' }
#'
#' @export
dp_score_plot <- function(
    X,
    center = TRUE,
    standardize = FALSE,
    g_dppca = FALSE,
    cpp.option = FALSE,
    axes = c(1, 2),
    eps_total,
    delta_total,
    method = c("add", "sparse"),
    frame = NULL,
    m_x = NULL,
    m_y = NULL,
    bin_method = c("WZ", "Lei", "none")
) {
  method <- match.arg(method, choices = c("add", "sparse"), several.ok = TRUE)
  bin_method <- match.arg(bin_method)
  color <- "#6A5ACD"

  score_res <- dp_score(
    X = X,
    center = center,
    standardize = standardize,
    g_dppca = g_dppca,
    cpp.option = cpp.option,
    axes = axes,
    eps_total = eps_total,
    delta_total = delta_total,
    method = method,
    frame = frame,
    m_x = m_x,
    m_y = m_y,
    bin_method = bin_method
  )

  X_score <- as.data.frame(score_res$score)
  colnames(X_score) <- c("pc_x", "pc_y")

  xlim <- score_res$frame$xlim
  ylim <- score_res$frame$ylim
  pc_names <- colnames(score_res$score)

  p_scatter <- ggplot2::ggplot(X_score, ggplot2::aes(x = .data$pc_x, y = .data$pc_y)) +
    ggplot2::geom_point(alpha = 0.6, size = 1.8, color = color) +
    ggplot2::coord_fixed(xlim = xlim, ylim = ylim) +
    ggplot2::scale_x_continuous(expand = c(0, 0), breaks = pretty(xlim, n = 5)) +
    ggplot2::scale_y_continuous(expand = c(0, 0), breaks = pretty(ylim, n = 5)) +
    theme_dp_base() +
    ggplot2::labs(x = pc_names[1], y = pc_names[2])

  p_scatter <- add_title_dp(p_scatter, "Original Scatter")

  p_none <- make_hist_plot_dp(score_res$none, xlim, ylim, color, "Original Hist")
  p_add <- NULL
  p_sparse <- NULL

  plot_panels <- list(p_scatter, p_none)

  if ("add" %in% method) {
    p_add <- make_hist_plot_dp(score_res$add, xlim, ylim, color, "Add DP Hist")
    plot_panels <- c(plot_panels, list(p_add))
  }

  if ("sparse" %in% method) {
    p_sparse <- make_hist_plot_dp(score_res$sparse, xlim, ylim, color, "Sparse DP Hist")
    plot_panels <- c(plot_panels, list(p_sparse))
  }

  p_all <- patchwork::wrap_plots(plot_panels, nrow = 1)

  list(
    score = score_res,
    plot = list(
      scatter = p_scatter,
      none = p_none,
      add = p_add,
      sparse = p_sparse,
      all = p_all
    )
  )
}

#' Group-wise DP score histograms
#'
#' @description
#' Computes pooled two-dimensional PCA scores and constructs group-wise empirical
#' and selected DP histograms on a common frame and grid.
#'
#' @param X Numeric matrix or data frame.
#' @param G Group labels, either a vector of length \code{nrow(X)} or a single
#'   column name in \code{X}.
#' @param center Logical; whether the variables should be centered before PCA.
#' @param standardize Logical; whether the variables should be scaled before PCA.
#' @param g_dppca Logical; whether to use a DP PCA direction matrix.
#' @param cpp.option Logical passed to \code{mech_tau_sph()} when private
#'   directions are computed.
#' @param axes Integer vector of length 2 specifying the principal components
#'   used for score construction.
#' @param eps_total Total privacy budget.
#' @param delta_total Total privacy parameter.
#' @param method Character vector specifying which DP histogram method(s) to
#'   compute. Options are \code{"add"} and \code{"sparse"}.
#' @param frame Optional user-specified frame.
#' @param m_x Optional number of bins along the x-axis.
#' @param m_y Optional number of bins along the y-axis.
#' @param bin_method Character string specifying the bin recommendation rule.
#'
#' @return A list with components \code{score}, \code{frame}, and \code{groups}.
#'
#' @importFrom VGAM rlaplace
#' @export
dp_score_group <- function(
    X,
    G,
    center = TRUE,
    standardize = FALSE,
    g_dppca = FALSE,
    cpp.option = FALSE,
    axes = c(1, 2),
    eps_total,
    delta_total,
    method = c("add", "sparse"),
    frame = NULL,
    m_x = NULL,
    m_y = NULL,
    bin_method = c("WZ", "Lei", "none")
) {
  method <- match.arg(method, choices = c("add", "sparse"), several.ok = TRUE)
  bin_method <- match.arg(bin_method)
  X <- as.data.frame(X)

  if (is.character(G) && length(G) == 1L) {
    if (!G %in% colnames(X)) stop("Column '", G, "' not found in X.")
    group_vec <- X[[G]]
    X_feat    <- X[, setdiff(colnames(X), G), drop = FALSE]
  } else {
    group_vec <- G
    X_feat    <- X
  }

  X_mat <- as.matrix(X_feat)

  if (!is.numeric(X_mat)) stop("Feature columns of X must be numeric.")
  if (nrow(X_mat) < 2) stop("X must have at least two rows.")
  if (ncol(X_mat) < 2) stop("X must have at least two feature columns.")
  if (length(group_vec) != nrow(X_mat)) stop("length(G) must equal nrow(X).")
  if (length(axes) != 2) stop("'axes' must be an integer vector of length 2.")
  if (any(!is.finite(axes)) || any(axes <= 0) || any(axes != as.integer(axes))) {
    stop("'axes' must contain positive integers.")
  }
  if (max(axes) > ncol(X_mat)) stop("The requested component index in 'axes' exceeds ncol(X).")
  if (!is.finite(eps_total) || eps_total <= 0) stop("'eps_total' must be a positive finite number.")
  if (!is.finite(delta_total) || delta_total <= 0 || delta_total >= 1) {
    stop("'delta_total' must be in (0, 1).")
  }

  axes <- as.integer(axes)
  k_max <- max(axes)
  g_levels <- as.character(unique(group_vec))

  budget <- .split_score_privacy_budget(
    eps_total = eps_total,
    delta_total = delta_total,
    g_dppca = g_dppca
  )

  X_proc <- prep_matrix_for_pca(
    X = X_mat,
    center = center,
    standardize = standardize
  )

  V_all <- dp_pc_dir(
    X = X_mat,
    k = k_max,
    center = center,
    standardize = standardize,
    g_dppca = g_dppca,
    eps_dir = budget$eps_dir,
    delta_dir = budget$delta_dir,
    cpp.option = cpp.option
  )

  V <- V_all[, axes, drop = FALSE]
  X_score <- as.matrix(X_proc %*% V)
  colnames(X_score) <- paste0("PC", axes)

  frame_out <- dp_frame(
    X       = X_score,
    eps_q   = budget$eps_q,
    delta_q = budget$delta_q,
    inflate = 0.20,
    q       = NULL,
    frame   = frame
  )

  xlim <- frame_out$xlim
  ylim <- frame_out$ylim

  if (bin_method == "none") {
    if (is.null(m_x) || is.null(m_y)) {
      stop("When bin_method = 'none', both m_x and m_y must be supplied.")
    }
  } else {
    if (is.null(m_x)) m_x <- recommend_bins(X_score, method = bin_method)
    if (is.null(m_y)) m_y <- recommend_bins(X_score, method = bin_method)
  }

  if (!is.numeric(m_x) || length(m_x) != 1 || m_x < 1) stop("'m_x' must be a positive integer.")
  if (!is.numeric(m_y) || length(m_y) != 1 || m_y < 1) stop("'m_y' must be a positive integer.")

  m_x <- as.integer(m_x)
  m_y <- as.integer(m_y)
  m   <- m_x * m_y

  x_breaks <- seq(xlim[1], xlim[2], length.out = m_x + 1)
  y_breaks <- seq(ylim[1], ylim[2], length.out = m_y + 1)

  base_coord <- do.call(
    rbind,
    lapply(seq_len(m_y), function(j) {
      data.frame(
        xmin = x_breaks[-length(x_breaks)],
        xmax = x_breaks[-1],
        ymin = y_breaks[j],
        ymax = y_breaks[j + 1]
      )
    })
  )

  groups_out <- list()

  eps_hist_method <- budget$eps_hist / length(method)
  delta_hist_method <- budget$delta_hist / length(method)

  for (g in g_levels) {
    idx_g <- which(as.character(group_vec) == g)
    Xg    <- X_score[idx_g, , drop = FALSE]
    n_g   <- nrow(Xg)

    bx <- cut(Xg[, 1], breaks = x_breaks, include.lowest = TRUE, labels = FALSE)
    by <- cut(Xg[, 2], breaks = y_breaks, include.lowest = TRUE, labels = FALSE)

    if (anyNA(bx)) {
      idx <- which(is.na(bx))
      bx[idx] <- findInterval(Xg[idx, 1], x_breaks, all.inside = TRUE)
    }
    if (anyNA(by)) {
      idx <- which(is.na(by))
      by[idx] <- findInterval(Xg[idx, 2], y_breaks, all.inside = TRUE)
    }

    bidx   <- (by - 1) * m_x + bx
    counts <- as.numeric(table(factor(bidx, levels = 1:m)))
    p_hat  <- counts / n_g

    hist_none <- base_coord
    hist_none$prob <- p_hat

    hist_add <- NULL
    hist_sparse <- NULL

    if ("add" %in% method) {
      sigma <- sqrt(2) * sqrt(2 * log(1.25 / delta_hist_method)) / eps_hist_method
      c_tilde <- pmax(counts + stats::rnorm(m, mean = 0, sd = sigma), 0)

      if (sum(c_tilde) <= 0) {
        stop("Group '", g, "': all privatized counts are zero after Gaussian noise.")
      }

      hist_add <- base_coord
      hist_add$prob <- c_tilde / sum(c_tilde)
    }

    if ("sparse" %in% method) {
      q_sparse <- numeric(m)
      scale_lap <- 2 / (eps_hist_method * n_g)
      thr <- (2 * log(2 / delta_hist_method)) / (eps_hist_method * n_g) + (1 / n_g)

      for (kk in seq_len(m)) {
        if (p_hat[kk] == 0) {
          q_sparse[kk] <- 0
        } else {
          z_k <- VGAM::rlaplace(1, location = 0, scale = scale_lap)
          q_k <- max(p_hat[kk] + z_k, 0)
          q_sparse[kk] <- if (q_k < thr) 0 else q_k
        }
      }

      if (sum(q_sparse) > 0) q_sparse <- q_sparse / sum(q_sparse)

      hist_sparse <- base_coord
      hist_sparse$prob <- q_sparse
    }

    groups_out[[g]] <- list(
      n      = n_g,
      none   = hist_none,
      add    = hist_add,
      sparse = hist_sparse
    )
  }

  list(
    score = X_score,
    frame = list(xlim = xlim, ylim = ylim),
    groups = groups_out,
    method = method
  )
}

#' Plot group-wise DP score histograms
#'
#' @description
#' Computes group-wise PCA score histograms through \code{dp_score_group()} and
#' returns both the calculation results and the corresponding plots.
#'
#' @param X Numeric matrix or data frame.
#' @param G Group labels, either a vector of length \code{nrow(X)} or a single
#'   column name in \code{X}.
#' @param center Logical; whether the variables should be centered before PCA.
#' @param standardize Logical; whether the variables should be scaled before PCA.
#' @param g_dppca Logical; whether to use a DP PCA direction matrix.
#' @param cpp.option Logical passed to \code{mech_tau_sph()} when private
#'   directions are computed.
#' @param axes Integer vector of length 2 specifying the principal components
#'   used for score construction.
#' @param eps_total Total privacy budget.
#' @param delta_total Total privacy parameter.
#' @param method Character vector specifying which DP histogram method(s) to
#'   plot. Options are \code{"add"} and \code{"sparse"}.
#' @param frame Optional user-specified frame.
#' @param m_x Optional number of bins along the x-axis.
#' @param m_y Optional number of bins along the y-axis.
#' @param bin_method Character string specifying the bin recommendation rule.
#'
#' @return A list with components \code{score}, \code{plot}, and
#'   \code{group_colors}.
#'
#' @export
dp_score_plot_group <- function(
    X,
    G,
    center = TRUE,
    standardize = FALSE,
    g_dppca = FALSE,
    cpp.option = FALSE,
    axes = c(1, 2),
    eps_total,
    delta_total,
    method = c("add", "sparse"),
    frame = NULL,
    m_x = NULL,
    m_y = NULL,
    bin_method = c("WZ", "Lei", "none")
) {
  method <- match.arg(method, choices = c("add", "sparse"), several.ok = TRUE)
  bin_method <- match.arg(bin_method)

  X_df <- as.data.frame(X)
  if (is.character(G) && length(G) == 1L) {
    if (!G %in% colnames(X_df)) stop("Column '", G, "' not found in X.")
    group_vec <- X_df[[G]]
    X_feat    <- X_df[, setdiff(colnames(X_df), G), drop = FALSE]
  } else {
    group_vec <- G
    X_feat    <- X_df
  }

  g_levels <- as.character(unique(group_vec))

  if (.is_color_vec(g_levels)) {
    col_map <- stats::setNames(as.character(g_levels), g_levels)
  } else {
    pal <- grDevices::hcl.colors(length(g_levels), "Dark3")
    col_map <- stats::setNames(pal, g_levels)
  }

  score_res <- dp_score_group(
    X = X_feat,
    G = group_vec,
    center = center,
    standardize = standardize,
    g_dppca = g_dppca,
    cpp.option = cpp.option,
    axes = axes,
    eps_total = eps_total,
    delta_total = delta_total,
    method = method,
    frame = frame,
    m_x = m_x,
    m_y = m_y,
    bin_method = bin_method
  )

  X_score <- as.data.frame(score_res$score)
  colnames(X_score) <- c("pc_x", "pc_y")
  X_score$group <- as.character(group_vec)

  xlim <- score_res$frame$xlim
  ylim <- score_res$frame$ylim
  pc_names <- colnames(score_res$score)

  p_scatter <- ggplot2::ggplot(
    X_score,
    ggplot2::aes(x = .data$pc_x, y = .data$pc_y, colour = .data$group)
  ) +
    ggplot2::geom_point(alpha = 0.6, size = 1.8) +
    ggplot2::scale_colour_manual(values = col_map) +
    ggplot2::coord_fixed(xlim = xlim, ylim = ylim) +
    ggplot2::scale_x_continuous(expand = c(0, 0), breaks = pretty(xlim, n = 5)) +
    ggplot2::scale_y_continuous(expand = c(0, 0), breaks = pretty(ylim, n = 5)) +
    theme_dp_base() +
    ggplot2::labs(x = pc_names[1], y = pc_names[2])

  p_scatter <- add_title_dp(p_scatter, "Original Scatter")

  coord_none_all <- dplyr::bind_rows(lapply(names(score_res$groups), function(g) {
    df <- score_res$groups[[g]]$none
    df$group <- g
    df
  }))

  p_none_all <- make_hist_all_dp(coord_none_all, xlim, ylim, col_map, "Original Hist")
  p_add_all <- NULL
  p_sparse_all <- NULL

  plot_panels <- list(p_scatter, p_none_all)

  if ("add" %in% method) {
    coord_add_all <- dplyr::bind_rows(lapply(names(score_res$groups), function(g) {
      df <- score_res$groups[[g]]$add
      df$group <- g
      df
    }))
    p_add_all <- make_hist_all_dp(coord_add_all, xlim, ylim, col_map, "Add DP Hist")
    plot_panels <- c(plot_panels, list(p_add_all))
  }

  if ("sparse" %in% method) {
    coord_sparse_all <- dplyr::bind_rows(lapply(names(score_res$groups), function(g) {
      df <- score_res$groups[[g]]$sparse
      df$group <- g
      df
    }))
    p_sparse_all <- make_hist_all_dp(coord_sparse_all, xlim, ylim, col_map, "Sparse DP Hist")
    plot_panels <- c(plot_panels, list(p_sparse_all))
  }

  p_all <- patchwork::wrap_plots(plot_panels, nrow = 1)

  make_group_layout <- function(which = c("none", "add", "sparse"), title_prefix) {
    which <- match.arg(which)

    group_plots <- lapply(g_levels, function(g) {
      df_g <- score_res$groups[[g]][[which]]
      title_g <- if (.is_color_vec(g_levels)) NULL else paste0(title_prefix, g)

      make_hist_single_dp(
        df_g,
        xlim = xlim,
        ylim = ylim,
        col  = col_map[[g]],
        title = title_g
      )
    })
    names(group_plots) <- g_levels

    list(
      all = patchwork::wrap_plots(group_plots, ncol = 3),
      groups = group_plots
    )
  }

  none_layout <- make_group_layout("none", "Original Hist: ")
  add_layout <- if ("add" %in% method) make_group_layout("add", "Add DP Hist: ") else NULL
  sparse_layout <- if ("sparse" %in% method) make_group_layout("sparse", "Sparse DP Hist: ") else NULL

  list(
    score = score_res,
    plot = list(
      scatter = p_scatter,
      none = none_layout,
      add = add_layout,
      sparse = sparse_layout,
      all = p_all
    ),
    group_colors = col_map
  )
}
