# ============================================================
# scree_plot_controls.R
# plotting controls for DP scree and PVE plots
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
#' @param plot_control Optional output of [scree_plot_control()].
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
    stop("`", arg, "` must be a single character string.", call. = FALSE)
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
    !is.numeric(x) || length(x) != 1L ||
    !is.finite(x) || x <= 0
  ) {
    stop("`", arg, "` must be a single positive number.", call. = FALSE)
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
    stop("`", arg, "` must be a named character vector.", call. = FALSE)
  }

  invisible(TRUE)
}
