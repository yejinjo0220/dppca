# ============================================================
# score_helpers.R
# Internal helper functions for score-related routines
# ============================================================

#' Differentially private quantile by smooth sensitivity
#'
#' Internal helper for estimating a univariate quantile with a smooth-sensitivity
#' based Laplace mechanism.
#'
#' @param x Numeric vector of observations.
#' @param q Quantile level in `(0, 1)`.
#' @param epsilon Positive `epsilon` privacy parameter.
#' @param delta Number in `(0, 1)` defining the `delta` privacy parameter.
#'
#' @return A numeric value giving the noisy quantile estimate.
#'
#' @noRd
dp_quantile_ss <- function(x, q, epsilon, delta) {
  x <- as.numeric(x)

  if (length(x) < 1L || anyNA(x) || any(!is.finite(x))) {
    stop("`x` must contain finite numeric values.", call. = FALSE)
  }
  if (!is.numeric(q) || length(q) != 1L || !is.finite(q) || q <= 0 || q >= 1) {
    stop("`q` must be a number in `(0, 1)`.", call. = FALSE)
  }
  if (
    !is.numeric(epsilon) || length(epsilon) != 1L ||
    !is.finite(epsilon) || epsilon <= 0
  ) {
    stop("`epsilon` must be a positive number.", call. = FALSE)
  }
  if (
    !is.numeric(delta) || length(delta) != 1L ||
    !is.finite(delta) || delta <= 0 || delta >= 1
  ) {
    stop("`delta` must be a number in `(0, 1)`.", call. = FALSE)
  }

  xs <- sort(x)
  n <- length(xs)
  r <- max(1L, min(n, ceiling(q * n)))
  beta <- epsilon / (2 * log(1 / delta))

  ss_beta <- 0
  for (k in 0:(n - 1L)) {
    t_vals <- 0:(k + 1L)
    left <- r + t_vals - (k + 1L)
    right <- r + t_vals
    valid <- (left >= 1L) & (right <= n)

    if (any(valid)) {
      a_k <- max(xs[right[valid]] - xs[left[valid]])
      val <- exp(-beta * k) * a_k
      if (val > ss_beta) {
        ss_beta <- val
      }
    }
  }

  scale <- (2 * ss_beta) / epsilon
  xs[r] + VGAM::rlaplace(1, location = 0, scale = scale)
}

#' Construct a private plotting frame for two-dimensional scores
#'
#' The plotting frame is constructed using coordinate-wise private medians as
#' the center and the private 0.99 quantile of Euclidean distances from that
#' center as the radius. The same radius is used on both axes, producing a
#' square plotting frame.
#'
#' @param X Numeric matrix with two columns.
#' @param eps_frame Positive `epsilon` privacy parameter for private frame
#'   construction.
#' @param delta_frame Number in `(0, 1)` defining the `delta` privacy
#'   parameter for private frame construction.
#' @param inflate Nonnegative inflation factor applied to the private radius.
#'
#' @return A list with components `xlim` and `ylim`.
#'
#' @noRd
dp_frame <- function(
    X,
    eps_frame,
    delta_frame,
    inflate = 0.20
) {
  X <- as.matrix(X)

  if (!is.numeric(X) || ncol(X) != 2L) {
    stop("`X` must be a numeric matrix with exactly two columns.", call. = FALSE)
  }
  if (nrow(X) < 2L) {
    stop("`X` must have at least two rows.", call. = FALSE)
  }
  if (anyNA(X) || any(!is.finite(X))) {
    stop("`X` must contain only finite values.", call. = FALSE)
  }
  if (
    !is.numeric(eps_frame) || length(eps_frame) != 1L ||
    !is.finite(eps_frame) || eps_frame <= 0
  ) {
    stop("`eps_frame` must be a positive number.", call. = FALSE)
  }
  if (
    !is.numeric(delta_frame) || length(delta_frame) != 1L ||
    !is.finite(delta_frame) || delta_frame <= 0 || delta_frame >= 1
  ) {
    stop("`delta_frame` must be a number in `(0, 1)`.", call. = FALSE)
  }
  if (
    !is.numeric(inflate) || length(inflate) != 1L ||
    !is.finite(inflate) || inflate < 0
  ) {
    stop("`inflate` must be a nonnegative number.", call. = FALSE)
  }

  eps_each <- eps_frame / 3
  delta_each <- delta_frame / 3
  q_radius <- 0.99

  center_x <- dp_quantile_ss(
    X[, 1], q = 0.5, epsilon = eps_each, delta = delta_each
  )
  center_y <- dp_quantile_ss(
    X[, 2], q = 0.5, epsilon = eps_each, delta = delta_each
  )

  radius_values <- sqrt((X[, 1] - center_x)^2 + (X[, 2] - center_y)^2)
  radius <- dp_quantile_ss(
    radius_values, q = q_radius, epsilon = eps_each, delta = delta_each
  )

  if (!is.finite(radius) || radius <= 0) {
    stop(
      "The private frame radius is not positive. ",
      "Try a larger privacy budget.",
      call. = FALSE
    )
  }

  inflated_radius <- (1 + inflate) * radius

  list(
    xlim = c(center_x - inflated_radius, center_x + inflated_radius),
    ylim = c(center_y - inflated_radius, center_y + inflated_radius)
  )
}

#' Add a centered title to a ggplot object
#'
#' @param plot A `ggplot` object.
#' @param title_text Plot title.
#'
#' @return A `ggplot` object.
#'
#' @noRd
add_title_dp <- function(plot, title_text) {
  if (is.null(plot) || is.null(title_text)) {
    return(plot)
  }

  plot +
    ggplot2::ggtitle(title_text) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        face = "bold",
        size = 14
      )
    )
}

#' Base theme for score plots
#'
#' @return A `ggplot2` theme object.
#'
#' @noRd
theme_dp_base <- function() {
  ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      legend.position = "none",
      panel.border = ggplot2::element_rect(
        color = "black",
        fill = NA,
        linewidth = 0.5
      ),
      plot.margin = ggplot2::margin(2, 2, 2, 2, unit = "pt")
    )
}

#' Plot a score histogram panel
#'
#' @param hist_df Histogram data frame with bin coordinates and probabilities.
#' @param xlim,ylim Plotting limits.
#' @param color Fill color.
#' @param title Optional plot title.
#' @param xlab,ylab Axis labels.
#'
#' @return A `ggplot` object or `NULL`.
#'
#' @noRd
make_hist_plot_dp <- function(hist_df, xlim, ylim, color, title = NULL,
                              xlab = "PC1", ylab = "PC2") {
  if (is.null(hist_df)) {
    return(NULL)
  }

  p <- ggplot2::ggplot(hist_df) +
    ggplot2::geom_rect(
      ggplot2::aes(
        xmin = .data$xmin,
        xmax = .data$xmax,
        ymin = .data$ymin,
        ymax = .data$ymax,
        alpha = .data$prob
      ),
      fill = color,
      linewidth = 0
    ) +
    ggplot2::scale_alpha(range = c(0, 1), guide = "none") +
    ggplot2::coord_fixed(xlim = xlim, ylim = ylim) +
    ggplot2::scale_x_continuous(expand = c(0, 0), breaks = pretty(xlim, n = 5)) +
    ggplot2::scale_y_continuous(expand = c(0, 0), breaks = pretty(ylim, n = 5)) +
    theme_dp_base() +
    ggplot2::labs(x = xlab, y = ylab)

  add_title_dp(p, title)
}

#' Plot a pooled group-wise histogram panel
#'
#' @param df Histogram data frame containing a `group` column.
#' @param xlim,ylim Plotting limits.
#' @param col_map Named color vector.
#' @param title Optional plot title.
#' @param xlab,ylab Axis labels.
#'
#' @return A `ggplot` object or a patchwork spacer.
#'
#' @noRd
make_hist_all_dp <- function(df, xlim, ylim, col_map, title = NULL,
                             xlab = "PC1", ylab = "PC2") {
  if (is.null(df)) {
    return(patchwork::plot_spacer())
  }

  df <- as.data.frame(df)
  required_cols <- c("xmin", "xmax", "ymin", "ymax", "prob")
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0L) {
    stop(
      "Histogram data frame is missing required column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  df$prob <- as.numeric(df$prob)
  df$prob[!is.finite(df$prob)] <- 0

  if ("group" %in% names(df)) {
    df$group <- as.factor(df$group)
    group_levels <- levels(df$group)
    col_map <- complete_color_map(col_map, group_levels)

    p <- ggplot2::ggplot(df) +
      ggplot2::geom_rect(
        ggplot2::aes(
          xmin = .data$xmin,
          xmax = .data$xmax,
          ymin = .data$ymin,
          ymax = .data$ymax,
          fill = .data$group,
          alpha = .data$prob
        ),
        linewidth = 0
      ) +
      ggplot2::scale_fill_manual(values = col_map, drop = FALSE) +
      ggplot2::scale_alpha_continuous(range = c(0, 0.85), guide = "none")
  } else {
    p <- ggplot2::ggplot(df) +
      ggplot2::geom_rect(
        ggplot2::aes(
          xmin = .data$xmin,
          xmax = .data$xmax,
          ymin = .data$ymin,
          ymax = .data$ymax,
          fill = .data$prob
        ),
        linewidth = 0
      ) +
      ggplot2::scale_fill_gradient(low = "white", high = "steelblue")
  }

  p <- p +
    ggplot2::coord_fixed(xlim = xlim, ylim = ylim) +
    ggplot2::scale_x_continuous(expand = c(0, 0), breaks = pretty(xlim, n = 5)) +
    ggplot2::scale_y_continuous(expand = c(0, 0), breaks = pretty(ylim, n = 5)) +
    theme_dp_base() +
    ggplot2::labs(x = xlab, y = ylab)

  add_title_dp(p, title)
}

#' Plot a single group histogram panel
#'
#' @param df Histogram data frame with bin coordinates and probabilities.
#' @param xlim,ylim Plotting limits.
#' @param col Fill color.
#' @param title Optional plot title.
#' @param xlab,ylab Axis labels.
#'
#' @return A `ggplot` object or a patchwork spacer.
#'
#' @noRd
make_hist_single_dp <- function(df, xlim, ylim, col, title = NULL,
                                xlab = "PC1", ylab = "PC2") {
  if (is.null(df)) {
    return(patchwork::plot_spacer())
  }

  p <- ggplot2::ggplot(df) +
    ggplot2::geom_rect(
      ggplot2::aes(
        xmin = .data$xmin,
        xmax = .data$xmax,
        ymin = .data$ymin,
        ymax = .data$ymax,
        alpha = .data$prob
      ),
      fill = col,
      linewidth = 0
    ) +
    ggplot2::scale_alpha(range = c(0, 1), guide = "none") +
    ggplot2::coord_fixed(xlim = xlim, ylim = ylim) +
    ggplot2::scale_x_continuous(expand = c(0, 0), breaks = pretty(xlim, n = 5)) +
    ggplot2::scale_y_continuous(expand = c(0, 0), breaks = pretty(ylim, n = 5)) +
    theme_dp_base() +
    ggplot2::labs(x = xlab, y = ylab)

  add_title_dp(p, title)
}

#' Complete a named color map
#'
#' @param col_map Optional color map.
#' @param group_levels Character vector of group levels.
#'
#' @return A named character vector of colors.
#'
#' @noRd
complete_color_map <- function(col_map, group_levels) {
  if (is.null(col_map) || length(col_map) == 0L) {
    return(stats::setNames(grDevices::hcl.colors(length(group_levels), "Dark 3"), group_levels))
  }

  col_map <- as.character(col_map)
  if (is.null(names(col_map))) {
    return(stats::setNames(rep_len(col_map, length(group_levels)), group_levels))
  }

  missing_groups <- setdiff(group_levels, names(col_map))
  if (length(missing_groups) > 0L) {
    extra_cols <- grDevices::hcl.colors(length(missing_groups), "Dark 3")
    col_map <- c(col_map, stats::setNames(extra_cols, missing_groups))
  }

  col_map[group_levels]
}

#' Check whether labels are valid colors
#'
#' @param x Character vector.
#'
#' @return A logical value.
#'
#' @noRd
.is_color_vec <- function(x) {
  x <- as.character(x)
  all(vapply(x, function(z) {
    out <- try(grDevices::col2rgb(z), silent = TRUE)
    !inherits(out, "try-error")
  }, logical(1)))
}

#' @importFrom rlang .data
NULL

# Internal helpers ------------------------------------------------------------

split_score_privacy_budget <- function(eps, delta, g_dppca) {
  if (isTRUE(g_dppca)) {
    list(
      eps_pc = eps / 3,
      eps_frame = eps / 3,
      eps_hist = eps / 3,
      delta_pc = delta / 3,
      delta_frame = delta / 3,
      delta_hist = delta / 3
    )
  } else {
    list(
      eps_pc = NULL,
      eps_frame = eps / 2,
      eps_hist = eps / 2,
      delta_pc = NULL,
      delta_frame = delta / 2,
      delta_hist = delta / 2
    )
  }
}

validate_score_matrix <- function(X) {
  X <- as.matrix(X)

  if (!is.numeric(X)) {
    stop("`X` must be numeric or coercible to a numeric matrix.", call. = FALSE)
  }
  if (nrow(X) < 2L) {
    stop("`X` must have at least two rows.", call. = FALSE)
  }
  if (ncol(X) < 2L) {
    stop("`X` must have at least two columns.", call. = FALSE)
  }
  if (anyNA(X) || any(!is.finite(X))) {
    stop("`X` must contain only finite values.", call. = FALSE)
  }

  X
}

validate_score_common <- function(
    X,
    eps,
    delta,
    bins,
    center,
    standardize,
    g_dppca,
    cpp.option,
    axes
) {
  validate_logical_value(center, "center")
  validate_logical_value(standardize, "standardize")
  validate_logical_value(g_dppca, "g_dppca")
  validate_logical_value(cpp.option, "cpp.option")

  if (!is.numeric(eps) || length(eps) != 1L || !is.finite(eps) || eps <= 0) {
    stop("`eps` must be a positive number.", call. = FALSE)
  }
  if (
    !is.numeric(delta) || length(delta) != 1L ||
    !is.finite(delta) || delta <= 0 || delta >= 1
  ) {
    stop("`delta` must be a number in `(0, 1)`.", call. = FALSE)
  }

  validate_bins(bins)

  if (length(axes) != 2L || !is.numeric(axes)) {
    stop("`axes` must be an integer vector of length 2.", call. = FALSE)
  }
  if (anyNA(axes) || any(!is.finite(axes)) || any(axes <= 0) ||
      any(axes != as.integer(axes))) {
    stop("`axes` must contain positive integers.", call. = FALSE)
  }
  if (max(axes) > ncol(X)) {
    stop("The largest value in `axes` must not exceed `ncol(X)`.", call. = FALSE)
  }

  invisible(TRUE)
}

validate_logical_value <- function(x, arg) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop("`", arg, "` must be `TRUE` or `FALSE`.", call. = FALSE)
  }

  invisible(TRUE)
}

validate_positive_integer <- function(x, arg) {
  if (
    !is.numeric(x) || length(x) != 1L || !is.finite(x) ||
    x < 1 || x != as.integer(x)
  ) {
    stop("`", arg, "` must be a positive integer.", call. = FALSE)
  }

  invisible(TRUE)
}

validate_bins <- function(bins) {
  if (
    !is.numeric(bins) || length(bins) != 2L || anyNA(bins) ||
    any(!is.finite(bins)) || any(bins < 1) || any(bins != as.integer(bins))
  ) {
    stop("`bins` must be an integer vector of length 2 with positive values.", call. = FALSE)
  }

  invisible(TRUE)
}

compute_score_coordinates <- function(
    X,
    axes,
    center,
    standardize,
    g_dppca,
    cpp.option,
    eps_pc,
    delta_pc
) {
  k_max <- max(axes)

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
    eps = eps_pc,
    delta = delta_pc,
    cpp.option = cpp.option
  )

  V <- V_all[, axes, drop = FALSE]
  X_score <- as.matrix(X_proc %*% V)
  colnames(X_score) <- paste0("PC", axes)

  list(score = X_score, directions = V)
}

score_histogram_grid <- function(xlim, ylim, m_x, m_y) {
  x_breaks <- seq(xlim[1], xlim[2], length.out = m_x + 1L)
  y_breaks <- seq(ylim[1], ylim[2], length.out = m_y + 1L)

  base_coord <- do.call(
    rbind,
    lapply(seq_len(m_y), function(j) {
      data.frame(
        xmin = x_breaks[-length(x_breaks)],
        xmax = x_breaks[-1],
        ymin = y_breaks[j],
        ymax = y_breaks[j + 1L]
      )
    })
  )

  list(
    x_breaks = x_breaks,
    y_breaks = y_breaks,
    base_coord = base_coord,
    m_x = m_x,
    m_y = m_y,
    m = m_x * m_y
  )
}

score_histograms <- function(
    X_score,
    xlim,
    ylim,
    bins,
    eps_hist,
    delta_hist,
    method
) {
  validate_bins(bins)
  bins <- as.integer(bins)
  grid <- score_histogram_grid(
    xlim = xlim,
    ylim = ylim,
    m_x = bins[1],
    m_y = bins[2]
  )
  eps_hist_method <- eps_hist / length(method)
  delta_hist_method <- delta_hist / length(method)

  score_histograms_from_grid(
    X_score = X_score,
    grid = grid,
    eps_hist_method = eps_hist_method,
    delta_hist_method = delta_hist_method,
    method = method,
    group_name = NULL
  )
}

score_histograms_from_grid <- function(
    X_score,
    grid,
    eps_hist_method,
    delta_hist_method,
    method,
    group_name = NULL
) {
  n <- nrow(X_score)
  if (n < 1L) {
    stop("Each histogram must contain at least one observation.", call. = FALSE)
  }

  bx <- cut(
    X_score[, 1],
    breaks = grid$x_breaks,
    include.lowest = TRUE,
    labels = FALSE
  )
  by <- cut(
    X_score[, 2],
    breaks = grid$y_breaks,
    include.lowest = TRUE,
    labels = FALSE
  )

  if (anyNA(bx)) {
    idx <- which(is.na(bx))
    bx[idx] <- findInterval(X_score[idx, 1], grid$x_breaks, all.inside = TRUE)
  }
  if (anyNA(by)) {
    idx <- which(is.na(by))
    by[idx] <- findInterval(X_score[idx, 2], grid$y_breaks, all.inside = TRUE)
  }

  bidx <- (by - 1L) * grid$m_x + bx
  counts <- as.numeric(table(factor(bidx, levels = seq_len(grid$m))))
  p_hat <- counts / n

  hist_none <- grid$base_coord
  hist_none$prob <- p_hat

  hist_add <- NULL
  hist_sparse <- NULL

  if ("add" %in% method) {
    sigma <- sqrt(2) * sqrt(2 * log(1.25 / delta_hist_method)) / eps_hist_method
    c_tilde <- pmax(counts + stats::rnorm(grid$m, mean = 0, sd = sigma), 0)

    if (sum(c_tilde) <= 0) {
      prefix <- if (is.null(group_name)) "" else paste0("Group `", group_name, "`: ")
      stop(
        prefix,
        "all privatized bin counts are zero after Gaussian noise and clipping. ",
        "Try a larger `eps` or fewer bins.",
        call. = FALSE
      )
    }

    hist_add <- grid$base_coord
    hist_add$prob <- c_tilde / sum(c_tilde)
  }

  if ("sparse" %in% method) {
    q_sparse <- numeric(grid$m)
    scale_lap <- 2 / (eps_hist_method * n)
    threshold <- (2 * log(2 / delta_hist_method)) / (eps_hist_method * n) + 1 / n

    for (kk in seq_len(grid$m)) {
      if (p_hat[kk] == 0) {
        q_sparse[kk] <- 0
      } else {
        z_k <- VGAM::rlaplace(1, location = 0, scale = scale_lap)
        q_k <- max(p_hat[kk] + z_k, 0)
        q_sparse[kk] <- if (q_k < threshold) 0 else q_k
      }
    }

    if (sum(q_sparse) > 0) {
      q_sparse <- q_sparse / sum(q_sparse)
    }

    hist_sparse <- grid$base_coord
    hist_sparse$prob <- q_sparse
  }

  list(
    none = hist_none,
    add = hist_add,
    sparse = hist_sparse
  )
}

#' Split features and group labels
#'
#' Separates the feature data from group labels for grouped score
#' visualizations. If `group` is a column name, that column is removed from `X`
#' and used as the group vector. Otherwise, `group` is treated as an external
#' group vector.
#'
#' @param X Input data frame or matrix.
#' @param group Group column name or group vector.
#'
#' @return A list with feature data `X` and group vector `group`.
#'
#' @noRd
split_group_input <- function(X, group) {
  X_df <- as.data.frame(X)

  if (is.character(group) && length(group) == 1L) {
    if (!group %in% colnames(X_df)) {
      stop("Column `", group, "` was not found in `X`.", call. = FALSE)
    }
    group_vec <- X_df[[group]]
    X_feat <- X_df[, setdiff(colnames(X_df), group), drop = FALSE]
  } else {
    group_vec <- group
    X_feat <- X_df
  }

  list(X = X_feat, group = group_vec)
}
