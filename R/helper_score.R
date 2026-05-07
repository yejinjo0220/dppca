# ============================================================
# helper_score_final.R
# Internal helper functions for score-related routines
# ============================================================

#' Differentially private quantile via smooth sensitivity
#'
#' @description
#' Internal helper that computes a differentially private estimate of a quantile
#' using the smooth sensitivity framework. This function is primarily used by
#' \code{dp_frame()} to construct a differentially private plotting frame for
#' score-based visualizations.
#'
#' @param x Numeric vector of observations.
#' @param q Target quantile level in \eqn{(0, 1)}.
#' @param epsilon Privacy budget \eqn{\varepsilon}.
#' @param delta Privacy parameter \eqn{\delta}.
#'
#' @return A numeric scalar giving the noisy quantile estimate.
#'
#' @keywords internal
dp_quantile_ss <- function(x, q, epsilon, delta) {
  x <- as.numeric(x)

  stopifnot(
    length(x) >= 1,
    is.finite(epsilon), epsilon > 0,
    is.finite(delta), delta > 0, delta < 1,
    is.finite(q), q > 0, q < 1
  )

  xs <- sort(x)
  n  <- length(xs)
  r  <- max(1L, min(n, ceiling(q * n)))

  beta <- epsilon / (2 * log(1 / delta))

  ss_beta <- 0
  for (k in 0:(n - 1)) {
    t_vals <- 0:(k + 1)
    left   <- r + t_vals - (k + 1)
    right  <- r + t_vals
    valid  <- (left >= 1) & (right <= n)

    if (any(valid)) {
      A_k <- max(xs[right[valid]] - xs[left[valid]])
      val <- exp(-beta * k) * A_k
      if (val > ss_beta) {
        ss_beta <- val
      }
    }
  }

  b <- (2 * ss_beta) / epsilon
  xs[r] + VGAM::rlaplace(1, location = 0, scale = b)
}

#' Plotting frame for 2D scores
#'
#' @description
#' Internal helper that constructs a plotting frame for two-dimensional score
#' data. If \code{frame} is supplied, it is used directly. Otherwise a
#' differentially private frame is constructed from privatized lower and upper
#' quantiles computed from the pooled score coordinates.
#'
#' @param X Numeric matrix with exactly two columns.
#' @param eps_q Privacy budget for DP quantile estimation. Required only when
#'   \code{frame = NULL}.
#' @param delta_q Privacy parameter for DP quantile estimation. Required only
#'   when \code{frame = NULL}.
#' @param inflate Non-negative numeric value controlling frame expansion.
#' @param q Optional symmetric quantile level in \eqn{(0, 1)}. If \code{NULL},
#'   extreme order statistics are used.
#' @param frame Optional user-specified frame. This can be either a numeric
#'   vector of length 2, interpreted as a common range for both axes, or a list
#'   with components \code{xlim} and \code{ylim}.
#'
#' @return A list with components \code{xlim} and \code{ylim}.
#'
#' @keywords internal
dp_frame <- function(X, eps_q = NULL, delta_q = NULL,
                     inflate = 0.10, q = NULL, frame = NULL) {
  stopifnot(is.matrix(X), ncol(X) == 2)

  if (!is.null(frame)) {
    if (is.numeric(frame) && length(frame) == 2) {
      if (!all(is.finite(frame)) || frame[1] >= frame[2]) {
        stop("If 'frame' is a numeric vector, it must satisfy frame[1] < frame[2].")
      }

      return(list(xlim = frame, ylim = frame))
    }

    if (is.list(frame) && !is.null(frame$xlim) && !is.null(frame$ylim)) {
      xlim <- frame$xlim
      ylim <- frame$ylim

      if (!is.numeric(xlim) || length(xlim) != 2 || !all(is.finite(xlim)) ||
          xlim[1] >= xlim[2]) {
        stop("frame$xlim must be a numeric vector of length 2 with xlim[1] < xlim[2].")
      }
      if (!is.numeric(ylim) || length(ylim) != 2 || !all(is.finite(ylim)) ||
          ylim[1] >= ylim[2]) {
        stop("frame$ylim must be a numeric vector of length 2 with ylim[1] < ylim[2].")
      }

      return(list(xlim = xlim, ylim = ylim))
    }

    stop("frame must be either a numeric vector of length 2 or a list with xlim and ylim.")
  }

  if (is.null(eps_q) || is.null(delta_q)) {
    stop("When 'frame = NULL', both 'eps_q' and 'delta_q' must be supplied.")
  }

  z <- as.numeric(X)
  n <- length(z)

  if (n < 2) stop("X must contain at least two values.")

  if (is.null(q)) {
    q_min <- 1 / n
    q_max <- (n - 1) / n
  } else {
    if (!is.numeric(q) || length(q) != 1 || !is.finite(q) || q <= 0 || q >= 1) {
      stop("'q' must be one number in (0, 1).")
    }
    q_min <- min(q, 1 - q)
    q_max <- max(q, 1 - q)
  }

  eps_each   <- eps_q / 2
  delta_each <- delta_q / 2

  dp_min <- dp_quantile_ss(z, q = q_min, epsilon = eps_each, delta = delta_each)
  dp_max <- dp_quantile_ss(z, q = q_max, epsilon = eps_each, delta = delta_each)

  if (dp_min > dp_max) {
    stop(
      "The privatized lower quantile exceeded the privatized upper quantile. ",
      "Please rerun dp_frame()."
    )
  }

  center   <- (dp_min + dp_max) / 2
  half_len <- (dp_max - dp_min) / 2
  inflated_half_len <- (1 + inflate) * half_len

  frame_min <- center - inflated_half_len
  frame_max <- center + inflated_half_len

  list(
    xlim = c(frame_min, frame_max),
    ylim = c(frame_min, frame_max)
  )
}

#' Recommend the number of histogram bins
#'
#' @description
#' Internal helper that recommends the number of bins per axis for a
#' two-dimensional histogram based on the sample size.
#'
#' @param X Matrix or data frame. Only the number of rows is used.
#' @param method Character string specifying the binning rule. Supported options
#'   are \code{"WZ"} and \code{"Lei"}.
#'
#' @return A positive integer giving the recommended number of bins per axis.
#'
#' @keywords internal
recommend_bins <- function(X, method = c("WZ", "Lei")) {
  X <- as.data.frame(X)
  n <- nrow(X)

  if (n < 2) stop("X must contain at least two rows.")

  method <- match.arg(method)

  m_axis <- switch(
    method,
    "WZ" = n^(1/4),
    "Lei" = (n / log(n))^(1/3)
  )

  max(1L, as.integer(round(m_axis)))
}

#' Add a centered title to a ggplot object
#'
#' @description
#' Internal plotting helper used in score-related plotting functions.
#'
#' @param plot A ggplot object.
#' @param title_text Plot title.
#'
#' @return A ggplot object.
#'
#' @keywords internal
add_title_dp <- function(plot, title_text) {
  if (is.null(plot) || is.null(title_text)) return(plot)

  plot +
    ggplot2::ggtitle(title_text) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        hjust = 0.5, face = "bold", size = 14
      )
    )
}

#' Base theme for score plots
#'
#' @description
#' Internal plotting helper that returns a minimal ggplot theme used across
#' score-related plots.
#'
#' @return A ggplot theme object.
#'
#' @keywords internal
theme_dp_base <- function() {
  ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      legend.position = "none",
      panel.border    = ggplot2::element_rect(
        color = "black", fill = NA, linewidth = 0.5
      ),
      plot.margin     = ggplot2::margin(2, 2, 2, 2, unit = "pt")
    )
}

#' Plot a single histogram panel
#'
#' @description
#' Internal plotting helper used by \code{dp_score_plot()} for a single
#' histogram panel.
#'
#' @param hist_df Histogram data frame with bin coordinates and probabilities.
#' @param xlim,ylim Plotting limits.
#' @param color Fill color.
#' @param title Optional plot title.
#'
#' @return A ggplot object.
#'
#' @keywords internal
make_hist_plot_dp <- function(hist_df, xlim, ylim, color, title = NULL) {
  if (is.null(hist_df)) return(NULL)

  p <- ggplot2::ggplot(hist_df) +
    ggplot2::geom_rect(
      ggplot2::aes(
        xmin = .data$xmin, xmax = .data$xmax,
        ymin = .data$ymin, ymax = .data$ymax,
        alpha = .data$prob
      ),
      fill = color, linewidth = 0
    ) +
    ggplot2::scale_alpha(range = c(0, 1), guide = "none") +
    ggplot2::coord_fixed(xlim = xlim, ylim = ylim) +
    ggplot2::scale_x_continuous(expand = c(0, 0), breaks = pretty(xlim, n = 5)) +
    ggplot2::scale_y_continuous(expand = c(0, 0), breaks = pretty(ylim, n = 5)) +
    theme_dp_base() +
    ggplot2::theme(axis.title = ggplot2::element_blank())

  add_title_dp(p, title)
}

#' Plot all-group histogram panel
#'
#' @description
#' Internal plotting helper used by \code{dp_score_plot_group()} for pooled
#' group-wise histogram displays.
#'
#' @param df Histogram data frame containing a \code{group} column.
#' @param xlim,ylim Plotting limits.
#' @param col_map Named color vector.
#' @param title Optional plot title.
#'
#' @return A ggplot object or a patchwork spacer.
#'
#' @keywords internal
make_hist_all_dp <- function(df, xlim, ylim, col_map, title = NULL) {
  if (is.null(df)) return(patchwork::plot_spacer())

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
    ggplot2::scale_fill_manual(values = col_map) +
    ggplot2::scale_alpha(range = c(0, 1), guide = "none") +
    ggplot2::coord_fixed(xlim = xlim, ylim = ylim) +
    ggplot2::scale_x_continuous(expand = c(0, 0), breaks = pretty(xlim, 5)) +
    ggplot2::scale_y_continuous(expand = c(0, 0), breaks = pretty(ylim, 5)) +
    theme_dp_base() +
    ggplot2::theme(axis.title = ggplot2::element_blank())

  add_title_dp(p, title)
}

#' Plot a single-group histogram panel
#'
#' @description
#' Internal plotting helper used by \code{dp_score_plot_group()} for one
#' group's histogram display.
#'
#' @param df Histogram data frame with bin coordinates and probabilities.
#' @param xlim,ylim Plotting limits.
#' @param col Fill color.
#' @param title Optional plot title.
#'
#' @return A ggplot object or a patchwork spacer.
#'
#' @keywords internal
make_hist_single_dp <- function(df, xlim, ylim, col, title = NULL) {
  if (is.null(df)) return(patchwork::plot_spacer())

  p <- ggplot2::ggplot(df) +
    ggplot2::geom_rect(
      ggplot2::aes(
        xmin = .data$xmin, xmax = .data$xmax,
        ymin = .data$ymin, ymax = .data$ymax,
        alpha = .data$prob
      ),
      fill = col, linewidth = 0
    ) +
    ggplot2::scale_alpha(range = c(0, 1), guide = "none") +
    ggplot2::coord_fixed(xlim = xlim, ylim = ylim) +
    ggplot2::scale_x_continuous(expand = c(0, 0), breaks = pretty(xlim, 5)) +
    ggplot2::scale_y_continuous(expand = c(0, 0), breaks = pretty(ylim, 5)) +
    theme_dp_base() +
    ggplot2::theme(axis.title = ggplot2::element_blank())

  add_title_dp(p, title)
}

#' Check whether labels are valid colors
#'
#' @description
#' Internal helper used by group plotting functions to determine whether group
#' labels can themselves be interpreted as colors.
#'
#' @param x Character vector.
#'
#' @return Logical scalar.
#'
#' @keywords internal
.is_color_vec <- function(x) {
  x <- as.character(x)
  all(vapply(x, function(z) {
    out <- try(grDevices::col2rgb(z), silent = TRUE)
    !inherits(out, "try-error")
  }, logical(1)))
}

#' @importFrom rlang .data
NULL
