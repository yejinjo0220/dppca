# helper functions for plot -------------------------------------------------------

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

theme_dp_base_dp <- function() {
  ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      legend.position = "none",
      panel.border    = ggplot2::element_rect(
        color = "black", fill = NA, linewidth = 0.5
      ),
      plot.margin     = ggplot2::margin(2, 2, 2, 2, unit = "pt")
    )
}

make_hist_plot_dp <- function(hist_df, xlim, ylim, color, title = NULL) {
  if (is.null(hist_df)) return(NULL)

  p <- ggplot2::ggplot(hist_df) +
    ggplot2::geom_rect(
      ggplot2::aes(
        xmin = xmin, xmax = xmax,
        ymin = ymin, ymax = ymax,
        alpha = prob
      ),
      fill = color, linewidth = 0
    ) +
    ggplot2::scale_alpha(range = c(0, 1), guide = "none") +
    ggplot2::coord_fixed(xlim = xlim, ylim = ylim) +
    ggplot2::scale_x_continuous(
      expand = c(0, 0),
      breaks = pretty(xlim, n = 5)
    ) +
    ggplot2::scale_y_continuous(
      expand = c(0, 0),
      breaks = pretty(ylim, n = 5)
    ) +
    theme_dp_base_dp() +
    ggplot2::theme(axis.title = ggplot2::element_blank())

  add_title_dp(p, title)
}

make_sample_plot_dp <- function(df, xlim, ylim, color, title = NULL) {
  if (is.null(df)) return(NULL)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_point(alpha = 0.6, size = 1.8, color = color) +
    ggplot2::coord_fixed(xlim = xlim, ylim = ylim) +
    ggplot2::scale_x_continuous(
      expand = c(0, 0),
      breaks = pretty(xlim, n = 5)
    ) +
    ggplot2::scale_y_continuous(
      expand = c(0, 0),
      breaks = pretty(ylim, n = 5)
    ) +
    theme_dp_base_dp() +
    ggplot2::theme(axis.title = ggplot2::element_blank())

  add_title_dp(p, title)
}

make_hist_all_dp <- function(df, xlim, ylim, col_map, title = NULL) {
  if (is.null(df)) return(patchwork::plot_spacer())

  p <- ggplot2::ggplot(df) +
    ggplot2::geom_rect(
      ggplot2::aes(
        xmin = xmin, xmax = xmax,
        ymin = ymin, ymax = ymax,
        fill = group, alpha = prob
      ),
      linewidth = 0
    ) +
    ggplot2::scale_fill_manual(values = col_map) +
    ggplot2::scale_alpha(range = c(0, 1), guide = "none") +
    ggplot2::coord_fixed(xlim = xlim, ylim = ylim) +
    ggplot2::scale_x_continuous(expand = c(0, 0), breaks = pretty(xlim, 5)) +
    ggplot2::scale_y_continuous(expand = c(0, 0), breaks = pretty(ylim, 5)) +
    theme_dp_base_dp() +
    ggplot2::theme(axis.title = ggplot2::element_blank())

  add_title_dp(p, title)
}

make_hist_single_dp <- function(df, xlim, ylim, col, title = NULL) {
  if (is.null(df)) return(patchwork::plot_spacer())

  p <- ggplot2::ggplot(df) +
    ggplot2::geom_rect(
      ggplot2::aes(
        xmin = xmin, xmax = xmax,
        ymin = ymin, ymax = ymax,
        alpha = prob
      ),
      fill = col, linewidth = 0
    ) +
    ggplot2::scale_alpha(range = c(0, 1), guide = "none") +
    ggplot2::coord_fixed(xlim = xlim, ylim = ylim) +
    ggplot2::scale_x_continuous(expand = c(0, 0), breaks = pretty(xlim, 5)) +
    ggplot2::scale_y_continuous(expand = c(0, 0), breaks = pretty(ylim, 5)) +
    theme_dp_base_dp() +
    ggplot2::theme(axis.title = ggplot2::element_blank())

  add_title_dp(p, title)
}

make_sample_all_dp <- function(df, xlim, ylim, col_map, title = NULL) {
  if (is.null(df)) return(patchwork::plot_spacer())

  p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, colour = group)) +
    ggplot2::geom_point(alpha = 0.6, size = 1.6) +
    ggplot2::scale_colour_manual(values = col_map) +
    ggplot2::coord_fixed(xlim = xlim, ylim = ylim) +
    ggplot2::scale_x_continuous(expand = c(0, 0), breaks = pretty(xlim, 5)) +
    ggplot2::scale_y_continuous(expand = c(0, 0), breaks = pretty(ylim, 5)) +
    theme_dp_base_dp() +
    ggplot2::theme(axis.title = ggplot2::element_blank())

  add_title_dp(p, title)
}

make_sample_single_dp <- function(df, xlim, ylim, col, title = NULL) {
  if (is.null(df)) return(patchwork::plot_spacer())

  p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_point(alpha = 0.6, size = 1.6, colour = col) +
    ggplot2::coord_fixed(xlim = xlim, ylim = ylim) +
    ggplot2::scale_x_continuous(expand = c(0, 0), breaks = pretty(xlim, 5)) +
    ggplot2::scale_y_continuous(expand = c(0, 0), breaks = pretty(ylim, 5)) +
    theme_dp_base_dp() +
    ggplot2::theme(axis.title = ggplot2::element_blank())

  add_title_dp(p, title)
}

.is_color_vec <- function(x) {
  x <- as.character(x)
  all(vapply(x, function(z) {
    out <- try(grDevices::col2rgb(z), silent = TRUE)
    !inherits(out, "try-error")
  }, logical(1)))
}


