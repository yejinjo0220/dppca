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


#' Construct a pure-DP plotting frame for two-dimensional scores
#'
#' The frame uses coordinate-wise pure-DP median estimates from
#' `unbounded_quantile()` as its center. For each coordinate separately,
#' `unbounded_quantile_upper()` estimates the 0.995 quantile of absolute
#' deviations from the corresponding private center. Each center and radius
#' mechanism receives `eps_frame / 4`, so the four mechanisms compose to a pure
#' `eps_frame`-DP frame under fixed-size replacement adjacency.
#'
#' Each private radius is multiplied by `1 + inflate`. Separate horizontal and
#' vertical radii produce a rectangular frame centered at the two private
#' median estimates.
#'
#' @param X Numeric matrix with exactly two columns and at least two rows. All
#'   entries must be finite.
#' @param eps_frame Positive total `epsilon` privacy parameter for the four
#'   private mechanisms used to construct the frame.
#' @param inflate Nonnegative inflation factor. Each private radius is
#'   multiplied by `1 + inflate`.
#'
#' @return A list with components `xlim` and `ylim`, each a length-two numeric
#'   vector giving the corresponding plotting limits.
#'
#' @noRd
dp_frame <- function(
    X,
    eps_frame,
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
    !is.numeric(inflate) || length(inflate) != 1L ||
    !is.finite(inflate) || inflate < 0
  ) {
    stop("`inflate` must be a nonnegative number.", call. = FALSE)
  }

  eps_each <- eps_frame / 4

  center_x <- unbounded_quantile(
    x = X[, 1], q = 0.5, epsilon = eps_each
  )
  center_y <- unbounded_quantile(
    x = X[, 2], q = 0.5, epsilon = eps_each
  )

  q_radius <- 0.995

  r_x_values <- abs( X[, 1] - center_x )
  r_x_max <- unbounded_quantile_upper(
    x = r_x_values, q = q_radius, epsilon = eps_each
  )

  r_y_values <- abs( X[, 2] - center_y )
  r_y_max <- unbounded_quantile_upper(
    x = r_y_values, q = q_radius, epsilon = eps_each
  )

  if (!is.finite(r_x_max) || r_x_max <= 0) {
    stop(
      "The private frame radius for x-axis is not positive. ",
      "Try a larger privacy budget.",
      call. = FALSE
    )
  }
  if (!is.finite(r_y_max) || r_y_max <= 0) {
    stop(
      "The private frame radius for y-axis is not positive. ",
      "Try a larger privacy budget.",
      call. = FALSE
    )
  }

  r_x_max <- (1 + inflate) * r_x_max
  r_y_max <- (1 + inflate) * r_y_max

  list(
    xlim = c(center_x - r_x_max, center_x + r_x_max),
    ylim = c(center_y - r_y_max, center_y + r_y_max)
  )
}

#' Add a centered title to a ggplot object
#'
#' @param plot A `ggplot` object.
#' @param title_text Plot title.
#' @param title_size Positive plot-title font size.
#'
#' @return A `ggplot` object.
#'
#' @noRd
add_title_dp <- function(plot, title_text, title_size = 14) {
  if (is.null(plot) || is.null(title_text)) {
    return(plot)
  }

  plot +
    ggplot2::ggtitle(title_text) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        face = "bold",
        size = title_size
      )
    )
}


#' Base theme for score plots
#'
#' @param base_size Positive base font size.
#' @param legend_position Legend position passed to `ggplot2::theme()`.
#'
#' @return A `ggplot2` theme object.
#'
#' @noRd
theme_dp_base <- function(base_size = 12, legend_position = "none") {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      legend.position = legend_position,
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
#' @param alpha_range Numeric alpha range for histogram cells.
#' @param base_size Positive base font size.
#' @param title_size Positive title font size.
#'
#' @return A `ggplot` object or `NULL`.
#'
#' @noRd
make_hist_plot_dp <- function(
    hist_df,
    xlim,
    ylim,
    color,
    title = NULL,
    xlab = "PC1",
    ylab = "PC2",
    alpha_range = c(0, 1),
    base_size = 12,
    title_size = 14
) {
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
    ggplot2::scale_alpha(range = alpha_range, guide = "none") +
    ggplot2::coord_fixed(xlim = xlim, ylim = ylim) +
    ggplot2::scale_x_continuous(
      expand = c(0, 0),
      breaks = pretty(xlim, n = 5)
    ) +
    ggplot2::scale_y_continuous(
      expand = c(0, 0),
      breaks = pretty(ylim, n = 5)
    ) +
    theme_dp_base(base_size = base_size) +
    ggplot2::labs(x = xlab, y = ylab)

  add_title_dp(p, title, title_size = title_size)
}


#' Plot synthetic score points sampled from a private histogram
#'
#' @param sample_df Data frame with columns `pc_x` and `pc_y`.
#' @param xlim,ylim Plotting limits.
#' @param color Point color.
#' @param title Optional plot title.
#' @param xlab,ylab Axis labels.
#' @param point_alpha Point transparency.
#' @param point_size Point size.
#' @param base_size Positive base font size.
#' @param title_size Positive title font size.
#'
#' @return A `ggplot` object.
#'
#' @noRd
make_sample_plot_dp <- function(
    sample_df,
    xlim,
    ylim,
    color,
    title = NULL,
    xlab = "PC1",
    ylab = "PC2",
    point_alpha = 0.6,
    point_size = 1.8,
    base_size = 12,
    title_size = 14
) {
  sample_df <- as.data.frame(sample_df)

  if (!all(c("pc_x", "pc_y") %in% names(sample_df))) {
    stop(
      "`sample_df` must contain columns `pc_x` and `pc_y`.",
      call. = FALSE
    )
  }

  p <- ggplot2::ggplot(
    sample_df,
    ggplot2::aes(x = .data$pc_x, y = .data$pc_y)
  ) +
    ggplot2::geom_point(
      alpha = point_alpha,
      size = point_size,
      color = color
    ) +
    ggplot2::coord_fixed(xlim = xlim, ylim = ylim) +
    ggplot2::scale_x_continuous(
      expand = c(0, 0),
      breaks = pretty(xlim, n = 5)
    ) +
    ggplot2::scale_y_continuous(
      expand = c(0, 0),
      breaks = pretty(ylim, n = 5)
    ) +
    theme_dp_base(base_size = base_size) +
    ggplot2::labs(x = xlab, y = ylab)

  add_title_dp(p, title, title_size = title_size)
}


#' Plot a pooled group-wise histogram panel
#'
#' @param df Histogram data frame containing a `group` column.
#' @param xlim,ylim Plotting limits.
#' @param col_map Named color vector.
#' @param title Optional plot title.
#' @param xlab,ylab Axis labels.
#' @param alpha_range Numeric alpha range for histogram cells.
#' @param base_size Positive base font size.
#' @param title_size Positive title font size.
#' @param legend_position Legend position.
#'
#' @return A `ggplot` object or a patchwork spacer.
#'
#' @noRd
make_hist_all_dp <- function(
    df,
    xlim,
    ylim,
    col_map,
    title = NULL,
    xlab = "PC1",
    ylab = "PC2",
    alpha_range = c(0, 0.85),
    base_size = 12,
    title_size = 14,
    legend_position = "none"
) {
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
      ggplot2::scale_alpha_continuous(
        range = alpha_range,
        guide = "none"
      )
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
    ggplot2::scale_x_continuous(
      expand = c(0, 0),
      breaks = pretty(xlim, n = 5)
    ) +
    ggplot2::scale_y_continuous(
      expand = c(0, 0),
      breaks = pretty(ylim, n = 5)
    ) +
    theme_dp_base(
      base_size = base_size,
      legend_position = legend_position
    ) +
    ggplot2::labs(x = xlab, y = ylab)

  add_title_dp(p, title, title_size = title_size)
}


#' Plot a single group histogram panel
#'
#' @param df Histogram data frame with bin coordinates and probabilities.
#' @param xlim,ylim Plotting limits.
#' @param col Fill color.
#' @param title Optional plot title.
#' @param xlab,ylab Axis labels.
#' @param alpha_range Numeric alpha range for histogram cells.
#' @param base_size Positive base font size.
#' @param title_size Positive title font size.
#'
#' @return A `ggplot` object or a patchwork spacer.
#'
#' @noRd
make_hist_single_dp <- function(
    df,
    xlim,
    ylim,
    col,
    title = NULL,
    xlab = "PC1",
    ylab = "PC2",
    alpha_range = c(0, 1),
    base_size = 12,
    title_size = 14
) {
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
    ggplot2::scale_alpha(range = alpha_range, guide = "none") +
    ggplot2::coord_fixed(xlim = xlim, ylim = ylim) +
    ggplot2::scale_x_continuous(
      expand = c(0, 0),
      breaks = pretty(xlim, n = 5)
    ) +
    ggplot2::scale_y_continuous(
      expand = c(0, 0),
      breaks = pretty(ylim, n = 5)
    ) +
    theme_dp_base(base_size = base_size) +
    ggplot2::labs(x = xlab, y = ylab)

  add_title_dp(p, title, title_size = title_size)
}



#' Plot grouped synthetic score points
#'
#' @param sample_df Data frame with columns `pc_x`, `pc_y`, and `group`.
#' @param xlim,ylim Plotting limits.
#' @param col_map Named vector mapping groups to colors.
#' @param title Optional plot title.
#' @param xlab,ylab Axis labels.
#' @param point_alpha Point transparency.
#' @param point_size Point size.
#' @param base_size Positive base font size.
#' @param title_size Positive title font size.
#' @param legend_position Legend position passed to `ggplot2::theme()`.
#'
#' @return A `ggplot` object.
#'
#' @noRd
make_group_sample_plot_dp <- function(
    sample_df,
    xlim,
    ylim,
    col_map,
    title = NULL,
    xlab = "PC1",
    ylab = "PC2",
    point_alpha = 0.6,
    point_size = 1.8,
    base_size = 12,
    title_size = 14,
    legend_position = "none"
) {
  sample_df <- as.data.frame(sample_df)

  required <- c("pc_x", "pc_y", "group")
  missing <- setdiff(required, names(sample_df))
  if (length(missing) > 0L) {
    stop(
      "`sample_df` is missing required column(s): ",
      paste(missing, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  sample_df$group <- factor(
    as.character(sample_df$group),
    levels = names(col_map)
  )

  p <- ggplot2::ggplot(
    sample_df,
    ggplot2::aes(
      x = .data$pc_x,
      y = .data$pc_y,
      colour = .data$group
    )
  ) +
    ggplot2::geom_point(
      alpha = point_alpha,
      size = point_size
    ) +
    ggplot2::scale_colour_manual(
      values = col_map,
      drop = FALSE
    ) +
    ggplot2::coord_fixed(xlim = xlim, ylim = ylim) +
    ggplot2::scale_x_continuous(
      expand = c(0, 0),
      breaks = pretty(xlim, n = 5)
    ) +
    ggplot2::scale_y_continuous(
      expand = c(0, 0),
      breaks = pretty(ylim, n = 5)
    ) +
    theme_dp_base(
      base_size = base_size,
      legend_position = legend_position
    ) +
    ggplot2::labs(
      x = xlab,
      y = ylab,
      colour = "Group"
    )

  add_title_dp(
    p,
    title,
    title_size = title_size
  )
}


#' Allocate a total synthetic sample size across groups
#'
#' Allocates a requested total number of synthetic observations in proportion
#' to observed group sizes using the largest-remainder rule. The returned
#' integer counts sum exactly to `total_size`.
#'
#' @param group_sizes Named numeric vector of nonnegative group sizes.
#' @param total_size Positive integer total synthetic sample size.
#'
#' @return Named nonnegative integer vector with the same names as
#'   `group_sizes`.
#'
#' @noRd
allocate_group_sample_sizes <- function(group_sizes, total_size) {
  if (
    !is.numeric(group_sizes) ||
    length(group_sizes) < 1L ||
    anyNA(group_sizes) ||
    any(!is.finite(group_sizes)) ||
    any(group_sizes < 0)
  ) {
    stop(
      "`group_sizes` must contain finite nonnegative values.",
      call. = FALSE
    )
  }

  if (is.null(names(group_sizes)) || any(names(group_sizes) == "")) {
    stop("`group_sizes` must be named.", call. = FALSE)
  }

  if (sum(group_sizes) <= 0) {
    stop("The total group size must be positive.", call. = FALSE)
  }

  if (
    !is.numeric(total_size) ||
    length(total_size) != 1L ||
    !is.finite(total_size) ||
    total_size < 1 ||
    total_size != as.integer(total_size)
  ) {
    stop("`total_size` must be a positive integer.", call. = FALSE)
  }

  total_size <- as.integer(total_size)

  raw <- total_size * group_sizes / sum(group_sizes)
  out <- floor(raw)
  remainder <- total_size - sum(out)

  if (remainder > 0L) {
    frac <- raw - out

    # Stable tie-breaking by original group order.
    ord <- order(
      -frac,
      seq_along(frac)
    )

    take <- ord[seq_len(remainder)]
    out[take] <- out[take] + 1L
  }

  out <- as.integer(out)
  names(out) <- names(group_sizes)

  if (any(out == 0L)) {
    warning(
      "The requested total `sample_size` is small relative to the number of ",
      "groups; at least one group receives zero synthetic observations.",
      call. = FALSE
    )
  }

  out
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
      eps_pc = 0.2 * eps,
      eps_frame = 0.2 * eps,
      eps_hist = 0.6 * eps,
      delta_pc = 0.2 * delta,
      delta_hist = 0.8 * delta
    )
  } else {
    list(
      eps_pc = NULL,
      eps_frame = 0.2 * eps,
      eps_hist = 0.8 * eps,
      delta_pc = NULL,
      delta_hist = delta
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

#' Compute score histograms on a privately constructed frame
#'
#' Internal helper that constructs the histogram grid and computes the
#' non-private histogram together with the requested private histogram methods.
#'
#' When multiple private methods are requested, the same `eps_hist` and
#' `delta_hist` values are supplied separately to each method for comparison.
#' The histogram privacy budget is not divided by the number of methods.
#' Therefore, if outputs from multiple private histogram methods are released
#' together, the privacy cost of the joint release must be handled by
#' composition.
#'
#' @noRd
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

  method <- unique(
    match.arg(
      method,
      choices = c("add", "sparse"),
      several.ok = TRUE
    )
  )

  grid <- score_histogram_grid(
    xlim = xlim,
    ylim = ylim,
    m_x = bins[1],
    m_y = bins[2]
  )

  # The same histogram privacy budget is supplied separately to each
  # requested method for method comparison. It is not divided across methods.
  eps_hist_method <- eps_hist
  delta_hist_method <- delta_hist

  score_histograms_from_grid(
    X_score = X_score,
    grid = grid,
    eps_hist_method = eps_hist_method,
    delta_hist_method = delta_hist_method,
    method = method,
    group_name = NULL
  )
}

#' Compute score histograms on a fixed grid
#'
#' Internal helper that computes the empirical non-private histogram and the
#' requested private histogram estimates on a common grid.
#'
#' The values `eps_hist_method` and `delta_hist_method` are applied separately
#' to each requested private histogram method. They are not divided across
#' `"add"` and `"sparse"`.
#'
#' @noRd
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

  method <- unique(
    match.arg(
      method,
      choices = c("add", "sparse"),
      several.ok = TRUE
    )
  )

  if (
    !is.numeric(eps_hist_method) ||
    length(eps_hist_method) != 1L ||
    !is.finite(eps_hist_method) ||
    eps_hist_method <= 0
  ) {
    stop("`eps_hist_method` must be a positive number.", call. = FALSE)
  }

  if (
    !is.numeric(delta_hist_method) ||
    length(delta_hist_method) != 1L ||
    !is.finite(delta_hist_method) ||
    delta_hist_method <= 0 ||
    delta_hist_method >= 1
  ) {
    stop("`delta_hist_method` must be a number in `(0, 1)`.", call. = FALSE)
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
    bx[idx] <- findInterval(
      X_score[idx, 1],
      grid$x_breaks,
      all.inside = TRUE
    )
  }

  if (anyNA(by)) {
    idx <- which(is.na(by))
    by[idx] <- findInterval(
      X_score[idx, 2],
      grid$y_breaks,
      all.inside = TRUE
    )
  }

  bidx <- (by - 1L) * grid$m_x + bx
  counts <- as.numeric(
    table(factor(bidx, levels = seq_len(grid$m)))
  )
  p_hat <- counts / n

  hist_nonprivate <- grid$base_coord
  hist_nonprivate$prob <- p_hat

  out <- list(
    nonprivate = hist_nonprivate
  )

  if ("add" %in% method) {
    sigma <- sqrt(2) *
      sqrt(2 * log(1.25 / delta_hist_method)) /
      eps_hist_method

    c_tilde <- pmax(
      counts + stats::rnorm(grid$m, mean = 0, sd = sigma),
      0
    )

    if (sum(c_tilde) <= 0) {
      prefix <- if (is.null(group_name)) {
        ""
      } else {
        paste0("Group `", group_name, "`: ")
      }

      stop(
        prefix,
        "all privatized bin counts are zero after Gaussian noise and clipping. ",
        "Try a larger `eps` or fewer bins.",
        call. = FALSE
      )
    }

    hist_add <- grid$base_coord
    hist_add$prob <- c_tilde / sum(c_tilde)
    out$add <- hist_add
  }

  if ("sparse" %in% method) {
    q_sparse <- numeric(grid$m)
    scale_lap <- 2 / (eps_hist_method * n)
    threshold <- (
      2 * log(2 / delta_hist_method)
    ) / (eps_hist_method * n) + 1 / n

    for (kk in seq_len(grid$m)) {
      if (p_hat[kk] == 0) {
        q_sparse[kk] <- 0
      } else {
        z_k <- VGAM::rlaplace(
          1,
          location = 0,
          scale = scale_lap
        )
        q_k <- max(p_hat[kk] + z_k, 0)
        q_sparse[kk] <- if (q_k < threshold) 0 else q_k
      }
    }

    if (sum(q_sparse) > 0) {
      q_sparse <- q_sparse / sum(q_sparse)
    }

    hist_sparse <- grid$base_coord
    hist_sparse$prob <- q_sparse
    out$sparse <- hist_sparse
  }

  out
}

#' Sample synthetic score points from a private histogram
#'
#' Internal post-processing helper for converting a private histogram into
#' synthetic score coordinates.
#'
#' @param hist_df Private histogram data frame containing bin boundaries and
#'   probabilities.
#' @param sample_size Positive integer number of synthetic points.
#' @param sample_method Character vector containing `"center"` and/or
#'   `"uniform"`.
#' @param bandwidth_scale Named nonnegative numeric vector with one value for
#'   each requested sampling method.
#'
#' @return A named list of data frames. Each data frame has columns `pc_x` and
#'   `pc_y`.
#'
#' @noRd
sample_private_score_histogram <- function(
    hist_df,
    sample_size,
    sample_method,
    bandwidth_scale
) {
  hist_df <- as.data.frame(hist_df)

  required <- c("xmin", "xmax", "ymin", "ymax", "prob")
  missing <- setdiff(required, names(hist_df))
  if (length(missing) > 0L) {
    stop(
      "`hist_df` is missing required column(s): ",
      paste(missing, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  if (
    !is.numeric(sample_size) ||
    length(sample_size) != 1L ||
    !is.finite(sample_size) ||
    sample_size < 1 ||
    sample_size != as.integer(sample_size)
  ) {
    stop("`sample_size` must be a positive integer.", call. = FALSE)
  }
  sample_size <- as.integer(sample_size)

  sample_method <- unique(
    match.arg(
      sample_method,
      choices = c("center", "uniform"),
      several.ok = TRUE
    )
  )

  if (
    !is.numeric(bandwidth_scale) ||
    anyNA(bandwidth_scale) ||
    any(!is.finite(bandwidth_scale)) ||
    any(bandwidth_scale < 0) ||
    is.null(names(bandwidth_scale))
  ) {
    stop(
      "`bandwidth_scale` must be a named vector of finite nonnegative values.",
      call. = FALSE
    )
  }

  missing_scale <- setdiff(sample_method, names(bandwidth_scale))
  if (length(missing_scale) > 0L) {
    stop(
      "`bandwidth_scale` is missing method(s): ",
      paste(missing_scale, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  probs <- as.numeric(hist_df$prob)
  probs[!is.finite(probs) | probs < 0] <- 0

  total_mass <- sum(probs)
  if (total_mass <= 0) {
    warning(
      "The private histogram contains no positive probability mass; ",
      "no synthetic score points were generated.",
      call. = FALSE
    )

    empty <- data.frame(
      pc_x = numeric(0),
      pc_y = numeric(0)
    )

    return(
      stats::setNames(
        rep(list(empty), length(sample_method)),
        sample_method
      )
    )
  }

  probs <- probs / total_mass

  dx <- hist_df$xmax - hist_df$xmin
  dy <- hist_df$ymax - hist_df$ymin

  if (
    any(!is.finite(dx)) ||
    any(!is.finite(dy)) ||
    any(dx <= 0) ||
    any(dy <= 0)
  ) {
    stop(
      "Histogram bin widths must be finite and positive.",
      call. = FALSE
    )
  }

  # The score histogram uses a regular rectangular grid. Median widths make
  # this helper robust to harmless floating-point differences.
  delta_x <- stats::median(dx)
  delta_y <- stats::median(dy)

  out <- list()

  for (sm in sample_method) {
    bin_id <- sample.int(
      n = nrow(hist_df),
      size = sample_size,
      replace = TRUE,
      prob = probs
    )

    if (sm == "center") {
      base_x <- (
        hist_df$xmin[bin_id] +
          hist_df$xmax[bin_id]
      ) / 2
      base_y <- (
        hist_df$ymin[bin_id] +
          hist_df$ymax[bin_id]
      ) / 2
    } else {
      base_x <- stats::runif(
        sample_size,
        min = hist_df$xmin[bin_id],
        max = hist_df$xmax[bin_id]
      )
      base_y <- stats::runif(
        sample_size,
        min = hist_df$ymin[bin_id],
        max = hist_df$ymax[bin_id]
      )
    }

    scale <- as.numeric(bandwidth_scale[[sm]])
    h_x <- scale * delta_x
    h_y <- scale * delta_y

    pc_x <- base_x + stats::rnorm(
      sample_size,
      mean = 0,
      sd = h_x
    )
    pc_y <- base_y + stats::rnorm(
      sample_size,
      mean = 0,
      sd = h_y
    )

    out[[sm]] <- data.frame(
      pc_x = pc_x,
      pc_y = pc_y
    )
  }

  out
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
