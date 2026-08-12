# ============================================================
# plot_controls.R
# plotting controls for scree/PVE and score plots
# ============================================================


# ============================================================
# Scree / PVE plot controls
# ============================================================

#' Control options for scree and PVE plots
#'
#' Creates a control list for the appearance of plots produced by
#' [dp_scree_plot()]. The plotting function uses base R graphics.
#'
#' @param title Optional plot title. If `NULL`, [dp_scree_plot()] uses
#'   `"DP PVE Plot"` when `type = "pve"` and `"DP Scree Plot"` when
#'   `type = "scree"`.
#' @param xlab Optional x-axis label. If `NULL`, `"Component"` is used.
#' @param ylab Optional y-axis label. If `NULL`, the label is chosen according
#'   to the plotted quantity.
#' @param legend_position Legend position passed to [graphics::legend()].
#'   The default is `"topright"`.
#' @param legend_labels Optional named character vector giving legend labels.
#'   Recognized names are `"nonprivate"`, `"clipped"`, `"pmwm"`, and `"huber"`.
#' @param col Optional named vector of plotting colors.
#' @param lty Optional named vector of line types.
#' @param pch Optional named vector of point symbols.
#' @param lwd Positive line width. The default is `1`.
#' @param point_cex Positive point-size multiplier. The default is `1`.
#' @param cex_main Positive title-size multiplier. The default is `1.5`.
#' @param cex_lab Positive axis-label-size multiplier. The default is `1.2`.
#' @param cex_axis Positive axis-text-size multiplier. The default is `1.1`.
#' @param cex_legend Positive legend-text-size multiplier. The default is `1.1`.
#' @param legend_bty Box type passed to [graphics::legend()]. The default
#'   `"n"` suppresses the legend border.
#'
#' @return A list of plotting options for [dp_scree_plot()].
#'
#' @examples
#' scree_plot_control()
#'
#' scree_plot_control(
#'   title = "PVE comparison",
#'   legend_position = "topright",
#'   cex_main = 1.4
#' )
#'
#' @export
scree_plot_control <- function(
    title = NULL,
    xlab = NULL,
    ylab = NULL,
    legend_position = "topright",
    legend_labels = NULL,
    col = NULL,
    lty = NULL,
    pch = NULL,
    lwd = 1,
    point_cex = 1,
    cex_main = 1.5,
    cex_lab = 1.2,
    cex_axis = 1.1,
    cex_legend = 1.1,
    legend_bty = "n"
) {
  .validate_plot_optional_string(title, "title")
  .validate_plot_optional_string(xlab, "xlab")
  .validate_plot_optional_string(ylab, "ylab")
  .validate_plot_string(legend_position, "legend_position")
  .validate_plot_string(legend_bty, "legend_bty")

  .validate_plot_positive_number(lwd, "lwd")
  .validate_plot_positive_number(point_cex, "point_cex")
  .validate_plot_positive_number(cex_main, "cex_main")
  .validate_plot_positive_number(cex_lab, "cex_lab")
  .validate_plot_positive_number(cex_axis, "cex_axis")
  .validate_plot_positive_number(cex_legend, "cex_legend")

  valid_series <- c("nonprivate", "clipped", "pmwm", "huber")

  if (!is.null(legend_labels)) {
    .validate_named_plot_vector(
      legend_labels,
      "legend_labels",
      valid_series,
      mode = "character"
    )
  }

  if (!is.null(col)) {
    .validate_named_plot_vector(
      col,
      "col",
      valid_series,
      mode = "character"
    )
  }

  if (!is.null(lty)) {
    .validate_named_plot_vector(
      lty,
      "lty",
      valid_series,
      mode = "any"
    )
  }

  if (!is.null(pch)) {
    .validate_named_plot_vector(
      pch,
      "pch",
      valid_series,
      mode = "any"
    )
  }

  list(
    title = title,
    xlab = xlab,
    ylab = ylab,
    legend_position = legend_position,
    legend_labels = legend_labels,
    col = col,
    lty = lty,
    pch = pch,
    lwd = lwd,
    point_cex = point_cex,
    cex_main = cex_main,
    cex_lab = cex_lab,
    cex_axis = cex_axis,
    cex_legend = cex_legend,
    legend_bty = legend_bty
  )
}


#' Default plotting options for scree and PVE plots
#'
#' @param type Plot type, either `"pve"` or `"scree"`.
#'
#' @return A complete default plotting-control list.
#' @noRd
.default_scree_plot_control <- function(type) {
  type <- match.arg(type, c("pve", "scree"))

  list(
    title = if (type == "pve") "DP PVE Plot" else "DP Scree Plot",
    xlab = "Component",
    ylab = if (type == "pve") {
      "Proportion of Variance Explained"
    } else {
      "Scree Value"
    },
    legend_position = "topright",
    legend_labels = c(
      nonprivate = "Non-private",
      clipped = "Clipped",
      pmwm = "PMWM",
      huber = "Huber"
    ),
    col = c(
      nonprivate = "black",
      clipped = "red",
      pmwm = "forestgreen",
      huber = "blue"
    ),
    lty = c(
      nonprivate = 1,
      clipped = 1,
      pmwm = 1,
      huber = 1
    ),
    pch = c(
      nonprivate = 16,
      clipped = 17,
      pmwm = 15,
      huber = 18
    ),
    lwd = 1,
    point_cex = 1,
    cex_main = 1.5,
    cex_lab = 1.2,
    cex_axis = 1.1,
    cex_legend = 1.1,
    legend_bty = "n"
  )
}


#' Merge user-supplied scree plot controls with defaults
#'
#' @param plot_control Optional output of `scree_plot_control()`.
#' @param type Plot type.
#'
#' @return A complete plotting-control list.
#' @noRd
.merge_scree_plot_control <- function(plot_control, type) {
  default <- .default_scree_plot_control(type)

  if (is.null(plot_control)) {
    return(default)
  }

  if (!is.list(plot_control)) {
    stop(
      "`plot_control` must be a list created by `scree_plot_control()`.",
      call. = FALSE
    )
  }

  out <- default

  scalar_names <- c(
    "title", "xlab", "ylab", "legend_position",
    "lwd", "point_cex", "cex_main", "cex_lab",
    "cex_axis", "cex_legend", "legend_bty"
  )

  for (nm in scalar_names) {
    if (!is.null(plot_control[[nm]])) {
      out[[nm]] <- plot_control[[nm]]
    }
  }

  named_vector_names <- c("legend_labels", "col", "lty", "pch")

  for (nm in named_vector_names) {
    if (!is.null(plot_control[[nm]])) {
      replacement <- plot_control[[nm]]
      out[[nm]][names(replacement)] <- replacement
    }
  }

  out
}


# ============================================================
# Score plot controls
# ============================================================

#' Control options for score plots
#'
#' Creates a control list for the appearance of plots produced by
#' [dp_score_plot()] and [dp_score_plot_group()]. Score plots are constructed
#' with `ggplot2` and combined with `patchwork`.
#'
#' @param color Single color used for non-group score plots. The default is
#'   `"#6A5ACD"`.
#' @param scatter_alpha Number in `[0, 1]` controlling point transparency in
#'   non-private and sampled scatter panels. The default is `0.6`.
#' @param scatter_size Positive number controlling point size in non-private
#'   and sampled scatter panels. The default is `1.8`.
#' @param hist_alpha_range Numeric vector of length 2 giving the minimum and
#'   maximum alpha values used for histogram cells. Values must satisfy
#'   `0 <= hist_alpha_range[1] <= hist_alpha_range[2] <= 1`. The default is
#'   `c(0, 1)`.
#' @param base_size Positive base font size for score plots. The default is `12`.
#' @param title_size Positive plot-title font size. The default is `14`.
#' @param scatter_title Optional title for the non-private scatter panel.
#'   The default is `"Non-private Scatter"`.
#' @param nonprivate_title Optional title for the non-private histogram panel.
#'   The default is `"Non-private Histogram"`.
#' @param add_title Optional title for the additive private histogram panel.
#'   The default is `"Add DP Histogram"`.
#' @param sparse_title Optional title for the sparse private histogram panel.
#'   The default is `"Sparse DP Histogram"`.
#' @param xlab Optional x-axis label. If `NULL`, the selected PC name is used.
#' @param ylab Optional y-axis label. If `NULL`, the selected PC name is used.
#' @param group_colors Optional named character vector of colors for grouped
#'   score plots. If `NULL`, group colors are determined automatically.
#' @param group_ncol Positive integer retained for backward compatibility.
#'   Grouped score plots are currently overlaid by group color rather than
#'   arranged as separate group-specific panels, so this option is not used by
#'   [dp_score_plot_group()]. The default is `3`.
#' @param show_group_legend Logical value indicating whether to display the
#'   group legend in grouped score plots. The default is `FALSE`.
#'
#' @return A list of plotting options for [dp_score_plot()] and
#'   [dp_score_plot_group()].
#'
#' @examples
#' score_plot_control()
#'
#' score_plot_control(
#'   scatter_size = 2.2,
#'   title_size = 13,
#'   nonprivate_title = "Non-private score histogram"
#' )
#'
#' @export
score_plot_control <- function(
    color = "#6A5ACD",
    scatter_alpha = 0.6,
    scatter_size = 1.8,
    hist_alpha_range = c(0, 1),
    base_size = 12,
    title_size = 14,
    scatter_title = "Non-private Scatter",
    nonprivate_title = "Non-private Histogram",
    add_title = "Add DP Histogram",
    sparse_title = "Sparse DP Histogram",
    xlab = NULL,
    ylab = NULL,
    group_colors = NULL,
    group_ncol = 3L,
    show_group_legend = FALSE
) {
  .validate_plot_color(color, "color")
  .validate_plot_probability(scatter_alpha, "scatter_alpha")
  .validate_plot_positive_number(scatter_size, "scatter_size")
  .validate_plot_alpha_range(hist_alpha_range, "hist_alpha_range")
  .validate_plot_positive_number(base_size, "base_size")
  .validate_plot_positive_number(title_size, "title_size")

  .validate_plot_optional_string(scatter_title, "scatter_title")
  .validate_plot_optional_string(nonprivate_title, "nonprivate_title")
  .validate_plot_optional_string(add_title, "add_title")
  .validate_plot_optional_string(sparse_title, "sparse_title")
  .validate_plot_optional_string(xlab, "xlab")
  .validate_plot_optional_string(ylab, "ylab")

  if (!is.null(group_colors)) {
    if (!is.character(group_colors) || length(group_colors) < 1L) {
      stop(
        "`group_colors` must be a non-empty character vector of colors.",
        call. = FALSE
      )
    }

    invalid <- vapply(
      group_colors,
      function(z) {
        inherits(try(grDevices::col2rgb(z), silent = TRUE), "try-error")
      },
      logical(1)
    )

    if (any(invalid)) {
      stop(
        "`group_colors` contains invalid color value(s).",
        call. = FALSE
      )
    }
  }

  .validate_plot_positive_integer(group_ncol, "group_ncol")
  .validate_plot_logical(show_group_legend, "show_group_legend")

  list(
    color = color,
    scatter_alpha = scatter_alpha,
    scatter_size = scatter_size,
    hist_alpha_range = as.numeric(hist_alpha_range),
    base_size = base_size,
    title_size = title_size,
    scatter_title = scatter_title,
    nonprivate_title = nonprivate_title,
    add_title = add_title,
    sparse_title = sparse_title,
    xlab = xlab,
    ylab = ylab,
    group_colors = group_colors,
    group_ncol = as.integer(group_ncol),
    show_group_legend = show_group_legend
  )
}


#' Default plotting options for score plots
#'
#' @return A complete default score plotting-control list.
#' @noRd
.default_score_plot_control <- function() {
  list(
    color = "#6A5ACD",
    scatter_alpha = 0.6,
    scatter_size = 1.8,
    hist_alpha_range = c(0, 1),
    base_size = 12,
    title_size = 14,
    scatter_title = "Non-private Scatter",
    nonprivate_title = "Non-private Histogram",
    add_title = "Add DP Histogram",
    sparse_title = "Sparse DP Histogram",
    xlab = NULL,
    ylab = NULL,
    group_colors = NULL,
    group_ncol = 3L,
    show_group_legend = FALSE
  )
}


#' Merge user-supplied score plot controls with defaults
#'
#' @param plot_control Optional output of `score_plot_control()`.
#'
#' @return A complete score plotting-control list.
#' @noRd
.merge_score_plot_control <- function(plot_control) {
  default <- .default_score_plot_control()

  if (is.null(plot_control)) {
    return(default)
  }

  if (!is.list(plot_control)) {
    stop(
      "`plot_control` must be a list created by `score_plot_control()`.",
      call. = FALSE
    )
  }

  out <- default
  valid_names <- names(default)
  supplied_names <- names(plot_control)

  if (is.null(supplied_names) || any(supplied_names == "")) {
    stop(
      "`plot_control` must be a named list created by `score_plot_control()`.",
      call. = FALSE
    )
  }

  unknown_names <- setdiff(supplied_names, valid_names)
  if (length(unknown_names) > 0L) {
    stop(
      "Unknown score plot control option(s): ",
      paste(unknown_names, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  for (nm in supplied_names) {
    # NULL is meaningful for optional labels and group_colors, so preserve it.
    out[[nm]] <- plot_control[[nm]]
  }

  out
}


# ============================================================
# Shared plot-control validation helpers
# ============================================================

#' Validate an optional character scalar
#'
#' @param x Object to validate.
#' @param arg Argument name.
#'
#' @return Invisibly returns `TRUE`.
#' @noRd
.validate_plot_optional_string <- function(x, arg) {
  if (is.null(x)) {
    return(invisible(TRUE))
  }

  .validate_plot_string(x, arg)
}


#' Validate a character scalar
#'
#' @param x Object to validate.
#' @param arg Argument name.
#'
#' @return Invisibly returns `TRUE`.
#' @noRd
.validate_plot_string <- function(x, arg) {
  if (!is.character(x) || length(x) != 1L || is.na(x)) {
    stop(
      "`", arg, "` must be a single character string.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}


#' Validate a positive numeric scalar used for plotting
#'
#' @param x Object to validate.
#' @param arg Argument name.
#'
#' @return Invisibly returns `TRUE`.
#' @noRd
.validate_plot_positive_number <- function(x, arg) {
  if (
    !is.numeric(x) ||
    length(x) != 1L ||
    !is.finite(x) ||
    x <= 0
  ) {
    stop(
      "`", arg, "` must be a single positive number.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}


#' Validate a probability-like plotting value
#'
#' @param x Object to validate.
#' @param arg Argument name.
#'
#' @return Invisibly returns `TRUE`.
#' @noRd
.validate_plot_probability <- function(x, arg) {
  if (
    !is.numeric(x) ||
    length(x) != 1L ||
    !is.finite(x) ||
    x < 0 ||
    x > 1
  ) {
    stop(
      "`", arg, "` must be a number in `[0, 1]`.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}


#' Validate a histogram alpha range
#'
#' @param x Numeric vector of length 2.
#' @param arg Argument name.
#'
#' @return Invisibly returns `TRUE`.
#' @noRd
.validate_plot_alpha_range <- function(x, arg) {
  if (
    !is.numeric(x) ||
    length(x) != 2L ||
    anyNA(x) ||
    any(!is.finite(x)) ||
    x[1] < 0 ||
    x[2] > 1 ||
    x[1] > x[2]
  ) {
    stop(
      "`", arg,
      "` must be a numeric vector `c(min, max)` with ",
      "`0 <= min <= max <= 1`.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}


#' Validate a color scalar
#'
#' @param x Color value.
#' @param arg Argument name.
#'
#' @return Invisibly returns `TRUE`.
#' @noRd
.validate_plot_color <- function(x, arg) {
  if (!is.character(x) || length(x) != 1L || is.na(x)) {
    stop(
      "`", arg, "` must be a single valid color string.",
      call. = FALSE
    )
  }

  ok <- !inherits(
    try(grDevices::col2rgb(x), silent = TRUE),
    "try-error"
  )

  if (!ok) {
    stop(
      "`", arg, "` must be a valid R color.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}


#' Validate a positive integer used for plotting
#'
#' @param x Object to validate.
#' @param arg Argument name.
#'
#' @return Invisibly returns `TRUE`.
#' @noRd
.validate_plot_positive_integer <- function(x, arg) {
  if (
    !is.numeric(x) ||
    length(x) != 1L ||
    !is.finite(x) ||
    x < 1 ||
    x != as.integer(x)
  ) {
    stop(
      "`", arg, "` must be a positive integer.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}


#' Validate a logical plotting option
#'
#' @param x Object to validate.
#' @param arg Argument name.
#'
#' @return Invisibly returns `TRUE`.
#' @noRd
.validate_plot_logical <- function(x, arg) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop(
      "`", arg, "` must be `TRUE` or `FALSE`.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}


#' Validate a named plotting vector
#'
#' @param x Vector to validate.
#' @param arg Argument name.
#' @param valid_names Allowed series names.
#' @param mode Required mode: `"character"` or `"any"`.
#'
#' @return Invisibly returns `TRUE`.
#' @noRd
.validate_named_plot_vector <- function(
    x,
    arg,
    valid_names,
    mode = c("character", "any")
) {
  mode <- match.arg(mode)

  if (is.null(names(x)) || any(names(x) == "")) {
    stop(
      "`", arg, "` must be a named vector with names among: ",
      paste(valid_names, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  if (any(!names(x) %in% valid_names)) {
    stop(
      "`", arg, "` contains unsupported name(s). Valid names are: ",
      paste(valid_names, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  if (mode == "character" && !is.character(x)) {
    stop(
      "`", arg, "` must be a named character vector.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}
