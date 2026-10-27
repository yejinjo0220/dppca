# ============================================================
# score_helpers.R
# Internal helper functions for score-related routines
# ============================================================

#' Construct a pure-DP plotting frame for two-dimensional scores
#'
#' The frame uses coordinate-wise pure-DP median estimates from
#' `unbounded_quantile()` as its center. `unbounded_quantile_upper()` then
#' estimates the 0.99 quantile of the Euclidean distances from that private
#' center. Each of the two center mechanisms and the radial-quantile mechanism
#' receives `eps_frame / 3`, so the three mechanisms compose to a pure
#' `eps_frame`-DP frame under fixed-size replacement adjacency of the score
#' rows. For a guarantee with respect to the original observations, directions
#' and preprocessing must satisfy the conditions in [dp_score()], and earlier
#' private releases must be composed with the frame mechanism.
#'
#' The common private radius is multiplied by `1 + inflate` and applied to both
#' coordinates. This produces an axis-aligned square centered at the two private
#' median estimates.
#'
#' @param X Numeric matrix with exactly two columns and at least two rows. All
#'   entries must be finite.
#' @param eps_frame Positive total `epsilon` privacy parameter for the three
#'   private mechanisms used to construct the frame.
#' @param inflate Nonnegative inflation factor. The common private radius is
#'   multiplied by `1 + inflate`.
#'
#' @return A list with components `xlim` and `ylim`, each a length-two numeric
#'   vector giving the corresponding plotting limits.
#'
#' @noRd
dp_frame <- function(
    X,
    eps_frame,
    inflate = 0
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

  eps_center_x <- eps_frame / 3
  eps_center_y <- eps_frame / 3
  eps_radius <- eps_frame / 3

  center_x <- unbounded_quantile(
    x = X[, 1], q = 0.5, epsilon = eps_center_x
  )
  center_y <- unbounded_quantile(
    x = X[, 2], q = 0.5, epsilon = eps_center_y
  )

  q_radius <- 0.99

  radius_values <- sqrt(
    (X[, 1] - center_x)^2 + (X[, 2] - center_y)^2
  )

  radius <- unbounded_quantile_upper(
    x = radius_values,
    q = q_radius,
    epsilon = eps_radius
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
#' @param plot A `ggplot` object or `NULL`.
#' @param title_text Optional plot title. If `NULL`, `plot` is returned
#'   unchanged.
#' @param title_size Positive plot-title font size.
#'
#' @return The titled `ggplot` object, or the original `plot` when either
#'   `plot` or `title_text` is `NULL`.
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
#' @param df Histogram data frame with bin coordinates and probabilities. An
#'   optional `group` column enables group-specific fills.
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
    # Show only groups with at least one nonzero histogram bin in the legend.
    # The stored histogram and displayed rectangles remain unchanged.
    active_groups <- unique(as.character(df$group[df$prob > 0]))

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
      ggplot2::scale_fill_manual(
        values = col_map, breaks = group_levels[group_levels %in% active_groups],
        drop = TRUE
      ) +
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
    # Omit groups with no sampled points from the sample legend.
    ggplot2::scale_colour_manual(
      values = col_map,
      drop = TRUE
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


#' Assign colors to groups
#'
#' Named colors are matched to group labels and missing groups receive an HCL
#' color. Unnamed colors are recycled in `group_levels` order.
#'
#' @param col_map Optional color vector, either named by group or in group order.
#'   If `NULL` or empty, an HCL palette is generated.
#' @param group_levels Character vector of group levels.
#'
#' @return A named character vector of colors.
#'
#' @noRd
complete_color_map <- function(col_map, group_levels) {
  if (is.null(col_map) || length(col_map) == 0L) {
    return(stats::setNames(grDevices::hcl.colors(length(group_levels), "Dark 3"), group_levels))
  }

  color_names <- names(col_map)
  col_map <- as.character(col_map)
  names(col_map) <- color_names
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

#' Extract loading and score-estimation budgets
#'
#' Use the loading allocation for directions, then split the score `eps`
#' 0.35/0.65 across frame and histogram estimation. All score `delta` is
#' assigned to the histogram. Frame construction is pure DP. Each requested
#' histogram method receives the full histogram portion.
#'
#' @param privacy Output of [privacy_control()] with `loading` and `score`
#'   components, plus an optional `scree` component for integrated PCA.
#'
#' @return A list with components `eps_pc`, `eps_frame`, `eps_hist`, `delta_pc`,
#'   and `delta_hist`.
#' @noRd
split_score_privacy_budget <- function(privacy) {
  components <- c("loading", "score")
  if (is.list(privacy) && "scree" %in% names(privacy$budgets)) {
    components <- c("loading", "scree", "score")
  }
  privacy <- .validate_privacy_control(privacy, components)
  budget <- privacy$budgets

  list(
    eps_pc = budget$loading[["eps"]],
    eps_frame = 0.35 * budget$score[["eps"]],
    eps_hist = 0.65 * budget$score[["eps"]],
    delta_pc = budget$loading[["delta"]],
    delta_hist = budget$score[["delta"]]
  )
}

#' Validate common score-estimation arguments
#'
#' @param X Numeric matrix used to check the available component indices.
#' @param bins Positive integer vector of length 2.
#' @param cpp.option Logical scalar option.
#' @param axes Positive integer vector of length 2 whose largest value does not
#'   exceed `ncol(X)`.
#'
#' @return Invisibly returns `TRUE`.
#' @noRd
validate_score_common <- function(
    X,
    bins,
    cpp.option,
    axes
) {
  validate_logical_value(cpp.option, "cpp.option")

  validate_bins(bins)

  if (length(axes) != 2L || !is.numeric(axes)) {
    stop("`axes` must be an integer vector of length 2.", call. = FALSE)
  }
  if (anyNA(axes) || any(!is.finite(axes)) || any(axes <= 0) ||
      any(axes != as.integer(axes))) {
    stop("`axes` must contain positive integers.", call. = FALSE)
  }
  if (anyDuplicated(axes)) {
    stop("`axes` must contain two distinct component indices.", call. = FALSE)
  }
  if (max(axes) > ncol(X)) {
    stop("The largest value in `axes` must not exceed `ncol(X)`.", call. = FALSE)
  }

  invisible(TRUE)
}

#' Validate two-dimensional histogram bin counts
#'
#' @param bins Object to validate as a positive integer vector of length 2.
#'
#' @return Invisibly returns `TRUE`.
#' @noRd
validate_bins <- function(bins) {
  if (
    !is.numeric(bins) || length(bins) != 2L || anyNA(bins) ||
    any(!is.finite(bins)) || any(bins < 1) || any(bins != as.integer(bins))
  ) {
    stop("`bins` must be an integer vector of length 2 with positive values.", call. = FALSE)
  }

  invisible(TRUE)
}

#' Compute selected principal component score coordinates
#'
#' @param X Numeric data matrix with observations in rows.
#' @param axes Positive integer vector of length 2 selecting score axes.
#' @param cpp.option Whether to use the Rcpp spherical-Kendall implementation.
#' @param eps_pc,delta_pc Privacy parameters used if directions are estimated;
#'   with supplied `V_dp`, they record its earlier loading allocation.
#' @param V_dp Optional full private direction matrix with orthonormal columns,
#'   matching feature order and preprocessing. See [dp_score()]. Supplying it
#'   skips direction estimation; the caller must establish its DP provenance.
#'
#' @return A list containing the observed projected `score` matrix and selected
#'   `directions`. These raw scores are needed internally and are not themselves
#'   a private release; public callers control whether to return them.
#' @noRd
compute_score_coordinates <- function(
    X,
    axes,
    cpp.option,
    eps_pc,
    delta_pc,
    V_dp = NULL
) {
  X_proc <- prep_matrix_for_pca(
    X = X
  )

  V_all <- .get_dp_directions(
    X_proc,
    c(eps = eps_pc, delta = delta_pc),
    cpp.option,
    V_dp
  )

  V <- V_all[, axes, drop = FALSE]

  X_score <- as.matrix(X_proc %*% V)
  colnames(X_score) <- paste0("PC", axes)

  list(
    score = X_score,
    directions = V
  )
}

#' Construct a regular two-dimensional histogram grid
#'
#' @param xlim,ylim Length-two numeric vectors giving grid limits.
#' @param m_x,m_y Positive integers giving the numbers of horizontal and
#'   vertical bins.
#'
#' @return A list containing `x_breaks`, `y_breaks`, the bin-boundary data frame
#'   `base_coord`, bin counts `m_x` and `m_y`, and total bin count `m`.
#' @noRd
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
#' Internal helper that constructs the histogram grid and computes requested
#' private histograms, with an optional empirical reference on the same grid.
#'
#' When multiple private methods are requested, the same `eps_hist` and
#' `delta_hist` values are supplied separately to each method for comparison.
#' The histogram privacy budget is not divided by the number of methods.
#' Therefore, if outputs from multiple private histogram methods are released
#' together, the privacy cost of the joint release must be handled by
#' composition.
#'
#' @param X_score Numeric matrix with two score columns.
#' @param xlim,ylim Length-two numeric vectors defining the histogram frame.
#' @param bins Positive integer vector of length 2 giving bin counts.
#' @param eps_hist Positive histogram `epsilon` privacy parameter.
#' @param delta_hist Number in `(0, 1)` defining the histogram `delta` privacy
#'   parameter.
#' @param method Character vector containing `"add"` and/or `"sparse"`.
#' @param non_private Whether to include the empirical histogram as `nonprivate`.
#'   The default is `TRUE`. `FALSE` omits it; internal counts are still computed.
#'
#' @return A named list containing requested `add` and/or `sparse` histogram
#'   data frames. The `nonprivate` data frame is included only when
#'   `non_private = TRUE` and describes the supplied score coordinates.
#' @noRd
score_histograms <- function(
    X_score,
    xlim,
    ylim,
    bins,
    eps_hist,
    delta_hist,
    method,
    non_private = TRUE
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
    group_name = NULL,
    non_private = non_private
  )
}

#' Compute score histograms on a fixed grid
#'
#' Internal helper that computes requested private histogram estimates and an
#' optional empirical reference on a common grid. Count-sensitivity arguments
#' require fixed directions, a fixed grid and suitable preprocessing; see
#' [dp_score()]. Omitting the reference does not skip count computation.
#'
#' The values `eps_hist_method` and `delta_hist_method` are applied separately
#' to each requested private histogram method. They are not divided across
#' `"add"` and `"sparse"`.
#'
#' @param X_score Numeric matrix with two score columns.
#' @param grid Histogram-grid list produced by `score_histogram_grid()`.
#' @param eps_hist_method Positive `epsilon` privacy parameter supplied to each
#'   requested histogram method.
#' @param delta_hist_method Number in `(0, 1)` defining the `delta` privacy
#'   parameter supplied to each requested histogram method.
#' @param method Character vector containing `"add"` and/or `"sparse"`.
#' @param group_name Optional group label used to prefix an error message.
#' @param non_private Whether to include an empirical `nonprivate` histogram.
#'   The default is `TRUE`; this reference is not a private release.
#'
#' @return A named list containing requested `add` and/or `sparse` histogram
#'   data frames, plus `nonprivate` only when `non_private = TRUE`. The empirical
#'   reference uses the supplied coordinates and does not compute ordinary PCA.
#'   Private histograms also retain a `dp_mass` attribute of unnormalized
#'   nonnegative noised counts for pooling already released group histograms.
#' @noRd
score_histograms_from_grid <- function(
    X_score,
    grid,
    eps_hist_method,
    delta_hist_method,
    method,
    group_name = NULL,
    non_private = TRUE
) {
  validate_logical_value(non_private, "non_private")

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

  out <- list()

  if (non_private) {
    hist_nonprivate <- grid$base_coord
    hist_nonprivate$prob <- p_hat
    out$nonprivate <- hist_nonprivate
  }

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
    attr(hist_add, "dp_mass") <- c_tilde
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

    # Save the noised count scale before probabilities are normalized.
    # This permits pooled displays to reuse the group releases.
    private_mass <- q_sparse * n
    if (sum(q_sparse) > 0) {
      q_sparse <- q_sparse / sum(q_sparse)
    }

    hist_sparse <- grid$base_coord
    hist_sparse$prob <- q_sparse
    attr(hist_sparse, "dp_mass") <- private_mass
    out$sparse <- hist_sparse
  }

  out
}

#' Sample synthetic score points from a histogram
#'
#' Internal helper for converting a histogram into synthetic score coordinates.
#' With a private histogram and public sampling choices, this is private
#' post-processing. It is also used for empirical comparison histograms by
#' [plot.dp_pca()]; sampling those does not make them private.
#'
#' @param hist_df Histogram data frame with columns `xmin`, `xmax`, `ymin`,
#'   `ymax` and `prob`. Nonfinite or negative probabilities are set to zero,
#'   and positive remaining mass is normalized before sampling. Bin widths
#'   must be positive and finite; the median widths determine jitter scales.
#' @param sample_size Positive integer number of synthetic points.
#' @param sample_method Character vector containing `"center"` and/or
#'   `"uniform"`.
#' @param bandwidth_scale Named nonnegative numeric vector with one value for
#'   each requested sampling method.
#'
#' @return A named list of data frames. Each data frame has columns `pc_x` and
#'   `pc_y`. If no positive probability mass remains, a warning is issued and
#'   each data frame has zero rows. Gaussian jitter can place points outside
#'   the histogram frame.
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


# Compute each group on one shared grid. The directions and frame have
# already been estimated by dp_score() or dp_pca().
.score_group_histograms <- function(
    scores,
    group,
    frame,
    bins,
    eps_hist,
    delta_hist,
    method,
    non_private
) {
  group <- as.character(group)
  labels <- unique(group)
  grid <- score_histogram_grid(
    xlim = frame$xlim, ylim = frame$ylim,
    m_x = bins[1L], m_y = bins[2L]
  )
  out <- stats::setNames(vector("list", length(labels)), labels)
  for (g in labels) {
    group_scores <- scores[group == g, , drop = FALSE]
    hist <- score_histograms_from_grid(
      X_score = group_scores, grid = grid,
      eps_hist_method = eps_hist, delta_hist_method = delta_hist,
      method = method, group_name = g, non_private = non_private
    )
    out[[g]] <- c(if (non_private) list(n = nrow(group_scores)), hist)
  }
  out
}

# Sum released noised counts rather than estimating a second pooled histogram.
# Empirical references, when retained, may use their observed group sizes.
.score_pool_histograms <- function(groups, method) {
  method <- unique(match.arg(method, c("add", "sparse"), several.ok = TRUE))
  out <- list()
  if (all(vapply(groups, function(g) !is.null(g$nonprivate), logical(1)))) {
    sizes <- vapply(groups, function(g) g$n, numeric(1))
    hist <- groups[[1L]]$nonprivate
    mass <- Reduce(`+`, lapply(seq_along(groups), function(i) {
      groups[[i]]$nonprivate$prob * sizes[i]
    }))
    hist$prob <- mass / sum(sizes)
    out$nonprivate <- hist
  }
  for (m in method) {
    hist <- groups[[1L]][[m]]
    masses <- lapply(groups, function(g) attr(g[[m]], "dp_mass"))
    if (any(vapply(masses, is.null, logical(1)))) {
      stop("Stored group histograms are missing private bin masses.", call. = FALSE)
    }
    mass <- Reduce(`+`, masses)
    hist$prob <- if (sum(mass) > 0) mass / sum(mass) else numeric(length(mass))
    attr(hist, "dp_mass") <- mass
    out[[m]] <- hist
  }
  out
}

# Render only stored histograms, frames and optional reference scores.
# The pooled and grouped S3 plot methods share this code.
.render_dp_score <- function(
    x,
    private_plot,
    sampling_control,
    plot_control,
    non_private,
    center = FALSE
) {
  validate_logical_value(center, "center")
  private_plot <- unique(match.arg(
    private_plot, c("histogram", "sample"), several.ok = TRUE
  ))
  validate_logical_value(non_private, "non_private")
  grouped <- inherits(x, "dp_score_group")
  metadata <- attr(x, "metadata")
  pctrl <- .merge_score_plot_control(plot_control)
  available_methods <- if (grouped) names(x$groups[[1L]]) else names(x)
  methods <- intersect(metadata$method, available_methods)
  if (!length(methods)) {
    stop("The stored score result contains no private histogram.", call. = FALSE)
  }

  has_reference <- if (grouped) {
    !is.null(x$groups[[1L]]$nonprivate)
  } else {
    !is.null(x$nonprivate)
  }
  if (non_private && !has_reference) {
    stop(
      "Non-private references were not retained; use `non_private = FALSE`.",
      call. = FALSE
    )
  }

  # Display-only centering: keep stored frames, bins and samples unchanged.
  private_shift <- if (center) {
    c(mean(x$frame$xlim), mean(x$frame$ylim))
  } else {
    c(0, 0)
  }
  xlim <- x$frame$xlim - private_shift[1L]
  ylim <- x$frame$ylim - private_shift[2L]
  reference_frame <- if (non_private) x$nonprivate_frame else NULL
  if (non_private && is.null(reference_frame)) {
    stop("The non-private frame is missing from the stored result.", call. = FALSE)
  }
  pc_names <- paste0("PC", metadata$axes)
  xlab <- if (is.null(pctrl$xlab)) pc_names[1L] else pctrl$xlab
  ylab <- if (is.null(pctrl$ylab)) pc_names[2L] else pctrl$ylab
  col_map <- NULL
  legend_position <- "none"

  if (grouped) {
    group_levels <- names(x$groups)
    col_map <- if (!is.null(pctrl$group_colors)) {
      complete_color_map(pctrl$group_colors, group_levels)
    } else if (.is_color_vec(group_levels)) {
      stats::setNames(group_levels, group_levels)
    } else {
      complete_color_map(NULL, group_levels)
    }
    legend_position <- if (pctrl$show_group_legend) "right" else "none"
  }

  # Keep each group legend as a single vertical column.
  legend_cols <- 1L
  format_group_plot <- function(pl) {
    if (!grouped || !pctrl$show_group_legend) return(pl)
    pl + ggplot2::theme(
      legend.position = "right",
      legend.text = ggplot2::element_text(size = max(8, 0.8 * pctrl$base_size)),
      legend.key.height = grid::unit(3, "mm"),
      legend.key.width = grid::unit(3, "mm")
    ) + ggplot2::guides(
      colour = ggplot2::guide_legend(ncol = legend_cols, byrow = FALSE),
      fill = ggplot2::guide_legend(ncol = legend_cols, byrow = FALSE)
    )
  }

  scatter <- NULL
  if (non_private && !is.null(x$score)) {
    observed <- as.data.frame(x$score)
    names(observed) <- c("pc_x", "pc_y")
    if (grouped) {
      observed$group <- factor(metadata$group, levels = group_levels)
      scatter <- make_group_sample_plot_dp(
        observed, reference_frame$xlim, reference_frame$ylim, col_map,
        title = pctrl$scatter_title, xlab = xlab, ylab = ylab,
        point_alpha = pctrl$scatter_alpha, point_size = pctrl$scatter_size,
        base_size = pctrl$base_size, title_size = pctrl$title_size,
        legend_position = legend_position
      ) + ggplot2::labs(colour = NULL)
      scatter <- format_group_plot(scatter)
    } else {
      scatter <- make_sample_plot_dp(
        observed, reference_frame$xlim, reference_frame$ylim, pctrl$color,
        title = pctrl$scatter_title, xlab = xlab, ylab = ylab,
        point_alpha = pctrl$scatter_alpha, point_size = pctrl$scatter_size,
        base_size = pctrl$base_size, title_size = pctrl$title_size
      )
    }
  }

  histograms <- list()
  hist_methods <- c(if (non_private) "nonprivate", methods)
  for (m in hist_methods) {
    title <- pctrl[[paste0(m, "_title")]]
    this_frame <- if (m == "nonprivate") {
      reference_frame
    } else {
      list(xlim = xlim, ylim = ylim)
    }
    if (grouped) {
      pieces <- lapply(group_levels, function(g) {
        df <- x$groups[[g]][[m]]
        if (m != "nonprivate") {
          df$xmin <- df$xmin - private_shift[1L]
          df$xmax <- df$xmax - private_shift[1L]
          df$ymin <- df$ymin - private_shift[2L]
          df$ymax <- df$ymax - private_shift[2L]
        }
        df$group <- rep(g, nrow(df))
        df
      })
      histograms[[m]] <- make_hist_all_dp(
        do.call(rbind, pieces), this_frame$xlim, this_frame$ylim, col_map,
        title = title, xlab = xlab, ylab = ylab,
        alpha_range = pctrl$hist_alpha_range, base_size = pctrl$base_size,
        title_size = pctrl$title_size, legend_position = legend_position
      ) + ggplot2::labs(fill = NULL)
      histograms[[m]] <- format_group_plot(histograms[[m]])
    } else {
      display_hist <- x[[m]]
      if (m != "nonprivate") {
        display_hist$xmin <- display_hist$xmin - private_shift[1L]
        display_hist$xmax <- display_hist$xmax - private_shift[1L]
        display_hist$ymin <- display_hist$ymin - private_shift[2L]
        display_hist$ymax <- display_hist$ymax - private_shift[2L]
      }
      histograms[[m]] <- make_hist_plot_dp(
        display_hist, this_frame$xlim, this_frame$ylim, pctrl$color,
        title = title, xlab = xlab, ylab = ylab,
        alpha_range = pctrl$hist_alpha_range, base_size = pctrl$base_size,
        title_size = pctrl$title_size
      )
    }
  }

  samples <- NULL
  sample_plots <- NULL
  reference_sample_plots <- NULL
  sampling_methods <- character(0)
  if ("sample" %in% private_plot) {
    sctrl <- .resolve_sampling_control(sampling_control, n = metadata$n)
    sampling_methods <- sctrl$method
    samples <- list()
    sample_plots <- list()

    # An explicit size requests sampling on *both* sides. With NULL, the
    # reference retains the actual centered PCA scores shown in `scatter`.
    sample_reference <- non_private && is.list(sampling_control) &&
      !is.null(sampling_control$sample_size)
    if (non_private) {
      reference_sample_plots <- stats::setNames(
        rep(list(scatter), length(sctrl$method)), sctrl$method
      )
    }
    if (sample_reference) {
      if (!grouped) {
        reference_samples <- sample_private_score_histogram(
          x$nonprivate, sctrl$sample_size, sctrl$method, sctrl$bandwidth_scale
        )
      }
      for (sm in sctrl$method) {
        if (grouped) {
          # Same sampling rule as dp_pca(): sample each reference group
          # proportional to its observed size, on the ordinary grid.
          one_method <- sctrl
          one_method$method <- sm
          sizes <- vapply(x$groups, function(g) g$n, numeric(1))
          ref <- .dppca_group_sample(
            x$groups, "nonprivate", one_method, group_sizes = sizes
          )
          names(ref)[1:2] <- c("pc_x", "pc_y")
          reference_sample_plots[[sm]] <- make_group_sample_plot_dp(
            ref, reference_frame$xlim, reference_frame$ylim, col_map,
            title = "Non-private Sample", xlab = xlab, ylab = ylab,
            point_alpha = pctrl$scatter_alpha, point_size = pctrl$scatter_size,
            base_size = pctrl$base_size, title_size = pctrl$title_size,
            legend_position = legend_position
          ) + ggplot2::labs(colour = NULL)
          reference_sample_plots[[sm]] <- format_group_plot(reference_sample_plots[[sm]])
        } else {
          reference_sample_plots[[sm]] <- make_sample_plot_dp(
            reference_samples[[sm]], reference_frame$xlim, reference_frame$ylim,
            pctrl$color, title = "Non-private Sample", xlab = xlab, ylab = ylab,
            point_alpha = pctrl$scatter_alpha, point_size = pctrl$scatter_size,
            base_size = pctrl$base_size, title_size = pctrl$title_size
          )
        }
      }
    }

    if (grouped) {
      # Use stored reference counts even when the reference plot is hidden.
      # Private-only results omit counts and retain equal allocation.
      group_sizes <- stats::setNames(rep(1, length(group_levels)), group_levels)
      if (has_reference) {
        group_sizes <- vapply(x$groups, function(g) g$n, numeric(1))
      }
      sample_sizes <- allocate_group_sample_sizes(group_sizes, sctrl$sample_size)
    }

    for (m in methods) {
      if (grouped) {
        pieces <- stats::setNames(vector("list", length(group_levels)), group_levels)
        for (g in group_levels) {
          if (sample_sizes[[g]] == 0L) {
            pieces[[g]] <- stats::setNames(
              lapply(sctrl$method, function(sm) {
                data.frame(pc_x = numeric(0), pc_y = numeric(0))
              }), sctrl$method
            )
          } else {
            pieces[[g]] <- sample_private_score_histogram(
              x$groups[[g]][[m]], sample_size = sample_sizes[[g]],
              sample_method = sctrl$method, bandwidth_scale = sctrl$bandwidth_scale
            )
          }
        }
        samples[[m]] <- stats::setNames(lapply(sctrl$method, function(sm) {
          combined <- do.call(rbind, lapply(group_levels, function(g) {
            df <- pieces[[g]][[sm]]
            df$group <- rep(g, nrow(df))
            df
          }))
          combined$group <- factor(combined$group, levels = group_levels)
          combined
        }), sctrl$method)
      } else {
        samples[[m]] <- sample_private_score_histogram(
          x[[m]], sample_size = sctrl$sample_size,
          sample_method = sctrl$method, bandwidth_scale = sctrl$bandwidth_scale
        )
      }

      sample_plots[[m]] <- list()
      for (sm in sctrl$method) {
        hist_label <- if (m == "add") "Add" else "Sparse"
        sample_label <- if (sm == "center") "Center" else "Uniform"
        title <- paste0(hist_label, " Sample (", sample_label, ")")
        display_sample <- samples[[m]][[sm]]
        display_sample$pc_x <- display_sample$pc_x - private_shift[1L]
        display_sample$pc_y <- display_sample$pc_y - private_shift[2L]
        if (grouped) {
          sample_plots[[m]][[sm]] <- make_group_sample_plot_dp(
            display_sample, xlim, ylim, col_map,
            title = title, xlab = xlab, ylab = ylab,
            point_alpha = pctrl$scatter_alpha, point_size = pctrl$scatter_size,
            base_size = pctrl$base_size, title_size = pctrl$title_size,
            legend_position = legend_position
          ) + ggplot2::labs(colour = NULL)
          sample_plots[[m]][[sm]] <- format_group_plot(sample_plots[[m]][[sm]])
        } else {
          sample_plots[[m]][[sm]] <- make_sample_plot_dp(
            display_sample, xlim, ylim, pctrl$color,
            title = title, xlab = xlab, ylab = ylab,
            point_alpha = pctrl$scatter_alpha, point_size = pctrl$scatter_size,
            base_size = pctrl$base_size, title_size = pctrl$title_size
          )
        }
      }
    }
  }

  plot_out <- c(if (!is.null(scatter)) list(scatter = scatter), histograms)
  if (!is.null(reference_sample_plots)) {
    plot_out$reference_samples <- reference_sample_plots
  }
  if (!is.null(sample_plots)) {
    plot_out$sample <- sample_plots
  }
  plot_out$all <- .score_plot_layout(
    scatter, histograms, sample_plots, methods, private_plot,
    sampling_methods, reference_sample_plots, grouped = grouped
  )
  out <- list(score = x)
  if (!is.null(samples)) {
    out$sample <- samples
  }
  out$plot <- plot_out
  if (grouped) {
    out$group_colors <- col_map
  }
  class(out) <- c(if (grouped) "dp_score_plot_group" else "dp_score_plot", "list")
  out
}

# Keep the same reference-column layout for both pooled and grouped displays.
.score_plot_layout <- function(
    scatter, histograms, samples, methods, private_plot, sampling_methods,
    reference_sample_plots = NULL, grouped = FALSE
) {
  use_samples <- "sample" %in% private_plot
  has_reference <- !is.null(scatter) || !is.null(histograms$nonprivate)
  keys <- c(if (has_reference) "reference", methods)
  columns <- list()

  # One method, histogram + one sampled panel, and no ordinary reference:
  # display the two plots side by side. Otherwise, keep one vertical stack
  # per method (ordinary reference, Add, Sparse).
  single_pair <- length(keys) == 1L &&
    length(private_plot) == 2L &&
    all(c("histogram", "sample") %in% private_plot) &&
    length(sampling_methods) == 1L

  for (key in keys) {
    method_panels <- list()
    if ("histogram" %in% private_plot) {
      if (key == "reference") {
        if (!use_samples && !is.null(scatter)) {
          method_panels[[length(method_panels) + 1L]] <- scatter
        }
        if (!is.null(histograms$nonprivate)) {
          method_panels[[length(method_panels) + 1L]] <- histograms$nonprivate
        }
      } else {
        method_panels[[length(method_panels) + 1L]] <- histograms[[key]]
      }
    }
    if (use_samples) {
      for (sm in sampling_methods) {
        panel <- if (key == "reference") {
          if (!is.null(reference_sample_plots)) reference_sample_plots[[sm]]
          else scatter
        } else {
          samples[[key]][[sm]]
        }
        if (!is.null(panel)) {
          method_panels[[length(method_panels) + 1L]] <- panel
        }
      }
    }
    if (!length(method_panels)) next

    # Display the group legend only alongside the last panel, i.e. the
    # sampled points when sampling was requested. The original per-panel
    # ggplot objects in the returned result are never modified.
    displayed <- lapply(seq_along(method_panels), function(i) {
      panel <- method_panels[[i]]
      if (grouped && i < length(method_panels) && inherits(panel, "ggplot")) {
        panel <- panel + ggplot2::theme(legend.position = "none")
      }
      panel
    })
    columns[[length(columns) + 1L]] <- patchwork::wrap_plots(
      displayed, ncol = if (single_pair) 2L else 1L
    )
  }
  if (!length(columns)) return(patchwork::plot_spacer())
  patchwork::wrap_plots(columns, ncol = length(columns))
}

# A short table rather than printing every bin and observed coordinate.
.score_summary <- function(object, grouped = FALSE) {
  metadata <- attr(object, "metadata")
  histograms <- if (grouped) object$groups else list(All = object)
  rows <- list()
  for (g in names(histograms)) {
    hist <- histograms[[g]]
    for (m in intersect(c("nonprivate", "add", "sparse"), names(hist))) {
      prob <- hist[[m]]$prob
      rows[[length(rows) + 1L]] <- data.frame(
        group = g, method = m, bins = length(prob),
        positive_bins = sum(prob > 0), probability_sum = sum(prob)
      )
    }
  }
  table <- do.call(rbind, rows)
  if (!grouped) {
    table$group <- NULL
  }
  out <- list(
    n = metadata$n, axes = metadata$axes, bins = metadata$bins,
    frame = object$frame, histograms = table, privacy = metadata$privacy,
    non_private = !is.null(object$score), grouped = grouped,
    group_name = metadata$group_name, group_levels = metadata$group_levels,
    excluded_columns = metadata$excluded_columns
  )
  class(out) <- "summary.dp_score"
  out
}

.score_plot_summary <- function(object) {
  out <- summary(object$score)
  counts <- list()
  for (m in names(object$sample)) {
    for (sm in names(object$sample[[m]])) {
      counts[[length(counts) + 1L]] <- data.frame(
        method = m, sampling = sm, points = nrow(object$sample[[m]][[sm]])
      )
    }
  }
  out$samples <- if (length(counts)) do.call(rbind, counts) else NULL
  out
}

# Number of columns used in a right-side group legend. Keep each column
# short enough to fit alongside a square score/biplot panel.
.dppca_group_legend_cols <- function(n) {
  n <- as.integer(n)
  if (n <= 7L) return(1L)
  if (n <= 14L) return(2L)
  max(3L, as.integer(ceiling(n / 7)))
}

# Infer the total number of known groups from an individual plot.
.dppca_plot_group_count <- function(plot) {
  if (!inherits(plot, "ggplot") || is.null(plot$data$group)) return(0L)
  group <- plot$data$group
  if (is.factor(group)) nlevels(group) else length(unique(group))
}

# Extract a right-positioned group guide from a plot without changing the
# original ggplot returned for individual extraction. A plot with no legend
# (for example an all-zero sparse sample) returns NULL.
.dppca_right_group_guide <- function(plot, group_count, title = NULL) {
  if (!inherits(plot, "ggplot")) return(NULL)
  if (identical(plot$theme$legend.position, "none")) return(NULL)
  guide_plot <- plot + ggplot2::theme(
    legend.position = "right",
    legend.text = ggplot2::element_text(size = 9),
    legend.title = ggplot2::element_text(size = 10, face = "bold"),
    legend.key.height = grid::unit(3, "mm"),
    legend.key.width = grid::unit(3, "mm")
  ) + ggplot2::guides(
    colour = ggplot2::guide_legend(ncol = 1L),
    fill = ggplot2::guide_legend(ncol = 1L)
  )
  if (!is.null(title)) {
    guide_plot <- guide_plot + ggplot2::labs(colour = title, fill = title)
  }
  g <- ggplot2::ggplotGrob(guide_plot)
  guides <- g$grobs[grepl("^guide-box", g$layout$name)]
  guides <- Filter(function(x) !inherits(x, "zeroGrob"), guides)
  if (length(guides)) guides[[1L]] else NULL
}
