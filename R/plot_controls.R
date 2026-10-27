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
#' [plot.dp_scree()]. The plotting function uses base R graphics.
#'
#' @param title Optional plot title. If `NULL`, [plot.dp_scree()] uses
#'   `"DP PVE Plot"` when `type = "pve"` and `"DP Scree Plot"` when
#'   `type = "scree"`.
#' @param xlab Optional x-axis label. If `NULL`, `"Component"` is used.
#' @param ylab Optional y-axis label. If `NULL`, `"PVE"` is used for PVE
#'   and `"Scree Value"` for scree values.
#' @param legend_position Legend position passed to [graphics::legend()].
#'   The default is `"topright"`.
#' @param legend_labels Optional named character vector giving legend labels.
#'   Recognized names are `"nonprivate"`, `"clipped"`, `"pmwm"`, and `"huber"`.
#' @param col Optional named character vector of plotting colors, using the
#'   same series names as `legend_labels`.
#' @param lty Optional named vector of line types, using the same series names
#'   as `legend_labels`.
#' @param pch Optional named vector of point symbols, using the same series
#'   names as `legend_labels`.
#' @param lwd Positive line width. The default is `1`.
#' @param point_cex Positive point-size multiplier. The default is `1`.
#' @param cex_main Positive title-size multiplier. The default is `1.5`.
#' @param cex_lab Positive axis-label-size multiplier. The default is `1.2`.
#' @param cex_axis Positive axis-text-size multiplier. The default is `1.1`.
#' @param cex_legend Positive legend-text-size multiplier. The default is `1.1`.
#' @param legend_bty Box type passed to [graphics::legend()]. The default
#'   `"n"` suppresses the legend border.
#'
#' @return A list of plotting options for [plot.dp_scree()].
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
      "PVE"
    } else {
      "Scree Value"
    },
    legend_position = "topright",
    legend_labels = c(
      nonprivate = "Non-Private",
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
#' [plot.dp_score()], including grouped score plots. Score plots are
#' constructed with `ggplot2` and combined with `patchwork`.
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
#'   The default is `"Non-Private Scatter"`.
#' @param nonprivate_title Optional title for the non-private histogram panel.
#'   The default is `"Non-Private Histogram"`.
#' @param add_title Optional title for the additive private histogram panel.
#'   The default is `"Add DP Histogram"`.
#' @param sparse_title Optional title for the sparse private histogram panel.
#'   The default is `"Sparse DP Histogram"`.
#' @param xlab Optional x-axis label. If `NULL`, the selected PC name is used.
#' @param ylab Optional y-axis label. If `NULL`, the selected PC name is used.
#' @param group_colors Optional nonempty character vector of colors for grouped
#'   score plots. Named colors are matched to group labels; missing groups
#'   receive generated colors. Unnamed colors are recycled in group order.
#'   If `NULL`, valid color labels are used directly, otherwise a palette is
#'   generated automatically.
#' @param group_ncol Positive integer retained for backward compatibility.
#'   Grouped score plots are currently overlaid by group color rather than
#'   arranged as separate group-specific panels, so this option is not used by
#'   [plot.dp_score()]. The default is `3`.
#' @param show_group_legend Logical value indicating whether to display the
#'   group legend in grouped score plots. The default is `TRUE`. In combined
#'   grouped displays, each histogram method uses one legend; standalone panels
#'   retain their own legends. Groups with no sampled points are omitted from
#'   sampled-point legends.
#'
#' @return A list of plotting options for [plot.dp_score()].
#'
#' @examples
#' score_plot_control()
#'
#' score_plot_control(
#'   scatter_size = 2.2,
#'   title_size = 13,
#'   nonprivate_title = "Non-Private score histogram"
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
    scatter_title = "Non-Private Scatter",
    nonprivate_title = "Non-Private Histogram",
    add_title = "Add DP Histogram",
    sparse_title = "Sparse DP Histogram",
    xlab = NULL,
    ylab = NULL,
    group_colors = NULL,
    group_ncol = 3L,
    show_group_legend = TRUE
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

  validate_positive_integer(group_ncol, "group_ncol")
  validate_logical_value(show_group_legend, "show_group_legend")

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
    scatter_title = "Non-Private Scatter",
    nonprivate_title = "Non-Private Histogram",
    add_title = "Add DP Histogram",
    sparse_title = "Sparse DP Histogram",
    xlab = NULL,
    ylab = NULL,
    group_colors = NULL,
    group_ncol = 3L,
    show_group_legend = TRUE
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

# Synthetic score sampling controls ----------------------------------------

#' Control options for synthetic score sampling
#'
#' Creates a control list for generating synthetic score points from a score
#' histogram. Sampling from a released private histogram is post-processing
#' under the conditions below; it does not alter the histogram mechanism.
#' [plot.dp_pca()] can also use these controls to sample an empirical reference
#' histogram, which does not make that reference private.
#'
#' @param method Character vector specifying the sampling construction.
#'   `"center"` selects a histogram bin according to its histogram probability
#'   and samples around that bin center using Gaussian jitter.
#'   `"uniform"` selects a histogram bin according to its histogram probability,
#'   samples uniformly inside the selected bin, and then adds Gaussian jitter.
#'   The default is `"uniform"`. Both methods can be requested with
#'   `c("center", "uniform")` in [plot.dp_score()];
#'   [plot.dp_pca()] requires exactly one sampling method.
#' @param sample_size Optional positive integer giving the synthetic sample
#'   size for each requested sampling construction. For [plot.dp_score()], this
#'   is the number of synthetic score points, or the total across all groups
#'   for a grouped fit. For [plot.dp_pca()], it is the size shared by score and biplot
#'   samples within each version; an explicit value also samples non-private
#'   reference histograms when comparison panels are requested. With `NULL`,
#'   private samples use the number of input observations and integrated
#'   non-private panels use observed scores. The default private sample size
#'   assumes that the input sample size is public. See [dp_score()] and
#'   [dp_pca()] for grouped sampling weights and the applicable privacy
#'   conditions.
#' @param bandwidth_scale Gaussian smoothing scale. If `NULL`, method-specific
#'   defaults are used: `0.5` for `"center"` and `0.25` for `"uniform"`.
#'   A single nonnegative number applies the same scale to all requested
#'   sampling methods. A named numeric vector can specify method-specific
#'   values, for example `c(center = 0.5, uniform = 0.2)`; it must include every
#'   requested method. A scale of zero disables Gaussian jitter.
#'
#' @details
#' For released histogram probabilities `p_k` and cells `B_k`, first draw
#' \deqn{J \sim \mathrm{Categorical}(p_1,\ldots,p_m).}
#' With `method = "center"`, sample `S^* = c_J + Z`, where `c_J` is the
#' selected cell center. With `method = "uniform"`, sample `U` uniformly in
#' `B_J` and use `S^* = U + Z`. In both cases Gaussian jitter has covariance
#' \deqn{\mathrm{Cov}(Z)=c^2\,\mathrm{diag}(\Delta_x^2,\Delta_y^2),}
#' where `c` is `bandwidth_scale` and the deltas are bin widths.
#'
#' Sampling uses only released histogram probabilities and bin geometry.
#' With publicly chosen settings it is post-processing and requires no new
#' privacy budget. Jitter can place points outside the frame; an all-zero
#' histogram yields no sampled points and a warning.
#'
#' @return A `dppca_sampling_control` list with components `method`,
#'   `sample_size`, and `bandwidth_scale`. The bandwidth values are resolved
#'   to a named numeric vector, while `sample_size = NULL` is resolved later
#'   by the plotting function.
#'
#' @examples
#' sampling_control()
#'
#' sampling_control(
#'   method = c("center", "uniform"),
#'   sample_size = 1000,
#'   bandwidth_scale = c(center = 0.5, uniform = 0.2)
#' )
#'
#' @export
sampling_control <- function(
    method = "uniform",
    sample_size = NULL,
    bandwidth_scale = NULL
) {
  method <- unique(
    match.arg(
      method,
      choices = c("center", "uniform"),
      several.ok = TRUE
    )
  )

  .validate_sampling_size(sample_size)
  scales <- .resolve_sampling_bandwidth(
    method = method,
    bandwidth_scale = bandwidth_scale
  )

  out <- list(
    method = method,
    sample_size = sample_size,
    bandwidth_scale = scales
  )

  class(out) <- c("dppca_sampling_control", "list")
  out
}


#' Default sampling control
#'
#' @return A sampling-control list using the package defaults.
#' @noRd
.default_sampling_control <- function() {
  sampling_control(
    method = "uniform",
    sample_size = NULL,
    bandwidth_scale = NULL
  )
}


#' Resolve and validate a sampling control
#'
#' @param control Optional sampling-control list.
#' @param n Number of input observations.
#'
#' @return A validated list with a resolved positive integer `sample_size`.
#' @noRd
.resolve_sampling_control <- function(control, n) {
  if (is.null(control)) {
    control <- .default_sampling_control()
  }

  if (!is.list(control)) {
    stop(
      "`sampling_control` must be a list created by `sampling_control()`.",
      call. = FALSE
    )
  }

  required <- c("method", "sample_size", "bandwidth_scale")
  missing <- setdiff(required, names(control))
  if (length(missing) > 0L) {
    stop(
      "`sampling_control` is missing component(s): ",
      paste(missing, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  method <- unique(
    match.arg(
      control$method,
      choices = c("center", "uniform"),
      several.ok = TRUE
    )
  )

  .validate_sampling_size(control$sample_size)

  if (
    !is.numeric(n) ||
    length(n) != 1L ||
    !is.finite(n) ||
    n < 1 ||
    n != as.integer(n)
  ) {
    stop("`n` must be a positive integer.", call. = FALSE)
  }

  sample_size <- if (is.null(control$sample_size)) {
    as.integer(n)
  } else {
    as.integer(control$sample_size)
  }

  scales <- .resolve_sampling_bandwidth(
    method = method,
    bandwidth_scale = control$bandwidth_scale
  )

  list(
    method = method,
    sample_size = sample_size,
    bandwidth_scale = scales
  )
}


#' Resolve method-specific bandwidth scales
#'
#' @param method Requested sampling methods.
#' @param bandwidth_scale User-supplied scale specification.
#'
#' @return Named numeric vector with one scale per requested method.
#' @noRd
.resolve_sampling_bandwidth <- function(method, bandwidth_scale) {
  defaults <- c(
    center = 0.5,
    uniform = 0.25
  )

  method <- unique(
    match.arg(
      method,
      choices = names(defaults),
      several.ok = TRUE
    )
  )

  if (is.null(bandwidth_scale)) {
    return(defaults[method])
  }

  if (
    !is.numeric(bandwidth_scale) ||
    length(bandwidth_scale) < 1L ||
    anyNA(bandwidth_scale) ||
    any(!is.finite(bandwidth_scale)) ||
    any(bandwidth_scale < 0)
  ) {
    stop(
      "`bandwidth_scale` must be `NULL` or contain finite nonnegative numbers.",
      call. = FALSE
    )
  }

  if (length(bandwidth_scale) == 1L) {
    return(
      stats::setNames(
        rep(as.numeric(bandwidth_scale), length(method)),
        method
      )
    )
  }

  nm <- names(bandwidth_scale)
  if (is.null(nm) || any(nm == "")) {
    stop(
      "When `bandwidth_scale` has length greater than one, it must be a named ",
      "vector with names `center` and/or `uniform`.",
      call. = FALSE
    )
  }

  if (anyDuplicated(nm)) {
    stop("`bandwidth_scale` must not contain duplicated names.", call. = FALSE)
  }

  invalid <- setdiff(nm, names(defaults))
  if (length(invalid) > 0L) {
    stop(
      "Unsupported `bandwidth_scale` name(s): ",
      paste(invalid, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  missing_method <- setdiff(method, nm)
  if (length(missing_method) > 0L) {
    stop(
      "`bandwidth_scale` must provide a value for each requested sampling ",
      "method. Missing: ",
      paste(missing_method, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  as.numeric(bandwidth_scale[method]) |>
    stats::setNames(method)
}


#' Validate synthetic sample size
#'
#' @param sample_size Sample-size specification.
#'
#' @return Invisibly returns `TRUE`.
#' @noRd
.validate_sampling_size <- function(sample_size) {
  if (is.null(sample_size)) {
    return(invisible(TRUE))
  }

  if (
    !is.numeric(sample_size) ||
    length(sample_size) != 1L ||
    !is.finite(sample_size) ||
    sample_size < 1 ||
    sample_size != as.integer(sample_size)
  ) {
    stop(
      "`sample_size` must be `NULL` or a positive integer.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}
