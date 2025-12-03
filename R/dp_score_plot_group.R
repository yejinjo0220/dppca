#' Group-wise DP PCA Score Plots with 2D Histograms
#'
#' @description
#' Creates group-wise PCA score plots combined with (optionally differentially
#' private) 2D histograms based on the first two principal components. Each
#' group shares a common DP frame and bin grid but has its own histogram (and
#' optional sampling) panels. Colors are assigned per group, either from the
#' group labels (if they are valid color codes) or from a qualitative palette.
#'
#' @description
#' Draw PCA scatter + DP hist panels for multiple groups.
#'
#' @param X A numeric matrix or data frame with \eqn{n} rows (observations) and
#'   \eqn{p} columns (variables and possibly a group label column). PCA is
#'   performed on the feature columns, and the first two PCs are used for
#'   plotting.
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
#'
#' @param eps_total Total privacy budget \eqn{\varepsilon_{\text{total}}} to be
#'   split among DP-PCA (optional), DP frame (\code{dp_frame()}), and group-wise
#'   DP histograms (inside \code{dp_hist_group()}).
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
#'   \code{dp_frame()}. A value of \code{inflate = 0.20} enlarges the interval
#'   by 20\%.
#' @param q_frame Optional symmetric quantile level in \eqn{(0,1)} passed to
#'   \code{dp_frame()} via \code{dp_hist_group()}. If \code{NULL}, extreme order
#'   statistics \eqn{1/n} and \eqn{(n-1)/n} are used.
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
#'   compute and display for each group:
#'   \itemize{
#'     \item \code{"none"}: (not used for group DP; empirical group histograms
#'           are derived from counts).
#'     \item \code{"add"}: Gaussian additive-noise DP histograms per group.
#'     \item \code{"sparse"}: sparse Laplace-thresholded DP histograms per group.
#'     \item \code{"all"}: compute and, if possible, display both
#'           \code{"add"} and \code{"sparse"}.
#'   }
#'
#' @param sampling Logical; if \code{TRUE}, additionally generate and display
#'   synthetic samples from each group's DP histogram(s) using
#'   \code{sample_from_hist()} with \code{k = n_g}, where \code{n_g} is the
#'   group size.
#'
#' @details
#' This function is a group-wise visualization wrapper around
#' \code{dp_hist_group()}. Internally, it:
#' \enumerate{
#'   \item Extracts feature columns and group labels from \code{X}, depending on
#'         \code{G}.
#'   \item Applies (optionally DP) PCA to the pooled data via \code{dp_pca()},
#'         and constructs a common DP frame and bin grid via
#'         \code{dp_hist_group()}.
#'   \item Builds a group-colored PCA scatter plot on the shared frame.
#'   \item Aggregates per-group DP histograms (and optionally samples) into
#'         panel layouts showing:
#'         \itemize{
#'           \item all-group histograms/samples (combined across groups),
#'           \item per-group histograms/samples arranged via patchwork.
#'         }
#'   \item Uses helper functions such as \code{theme_dp_base_dp()},
#'         \code{add_title_dp()}, \code{make_hist_all_dp()},
#'         \code{make_hist_single_dp()}, \code{make_sample_all_dp()},
#'         and \code{make_sample_single_dp()} to construct consistent ggplot2
#'         panels.
#' }
#'
#' Group colors are chosen as follows: if the group labels themselves are valid
#' color codes (checked via an internal helper \code{.is_color_vec()}), they are
#' used directly. Otherwise, colors are drawn from \code{hcl.colors()} with a
#' qualitative palette (e.g., \code{"Dark3"}).
#'
#' @return
#' A list equal to the output of \code{dp_hist_group()}, augmented with:
#' \itemize{
#'   \item \code{plot}: a list of ggplot2/patchwork objects, including
#'     \itemize{
#'       \item \code{all}: the main combined layout showing the group-colored
#'             scatter plot and group-wise histograms (and, optionally, samples).
#'       \item \code{none_hist}, \code{add_hist}, \code{sparse_hist}:
#'             each is a list whose first element \code{all} is a patchwork
#'             layout, followed by per-group histogram panels.
#'       \item \code{none_sampling}, \code{add_sampling},
#'             \code{sparse_sampling}: analogous lists for sampling-based
#'             visualizations (if \code{sampling = TRUE}).
#'     }
#'   \item \code{group_colors}: a named character vector mapping group levels
#'     to the colors used in the plots.
#' }
#'
#' @seealso
#'   \code{\link{dp_score_plot}},
#'   \code{\link{dp_hist_group}},
#'   \code{\link{dp_pca}},
#'   \code{\link{dp_frame}},
#'   \code{\link{number_bins}},
#'   \code{\link{sample_from_hist}}
#'
#' @importFrom ggplot2 ggplot aes geom_point coord_fixed scale_x_continuous
#'   scale_y_continuous scale_colour_manual labs
#' @importFrom dplyr bind_rows
#' @importFrom patchwork wrap_plots plot_spacer
#' @importFrom grDevices hcl.colors
#' @importFrom stats setNames
#' @export
dp_score_plot_group <- function(
    X,
    G,
    # PCA option
    center = TRUE,
    scale. = FALSE,
    dp_pca_flag = FALSE,
    cpp.option = FALSE,

    # Privacy budget
    eps_total,
    delta_total,
    eps_ratio   = NULL,
    delta_ratio = NULL,

    # Frame option
    inflate = 0.20,
    q_frame = NULL,

    # bin count
    m_x = NULL,
    m_y = NULL,
    bin_method = c("J", "W"),

    mechanism = c("all", "none", "add", "sparse"),

    # Sampling
    sampling  = FALSE
) {
  bin_method <- match.arg(bin_method)
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
  if (p < 2) stop("X (excluding group column) must have at least two columns for PCA.")

  groups   <- unique(group_vec)
  g_levels <- as.character(groups)

  # Color mapping
  if (.is_color_vec(g_levels)) {
    col_map <- stats::setNames(as.character(g_levels), g_levels)
  } else {
    pal <- grDevices::hcl.colors(length(g_levels), "Dark3")
    col_map <- stats::setNames(pal, g_levels)
  }

  # dp_hist_group
  hist_res <- dp_hist_group(
    X            = X_mat,
    G            = group_vec,
    center       = center,
    scale.       = scale.,
    dp_pca_flag  = dp_pca_flag,
    cpp.option   = cpp.option,
    eps_total    = eps_total,
    delta_total  = delta_total,
    eps_ratio    = eps_ratio,
    delta_ratio  = delta_ratio,
    inflate      = inflate,
    q_frame      = q_frame,
    m_x          = m_x,
    m_y          = m_y,
    bin_method   = bin_method,
    mechanism    = mechanism,
    sampling     = sampling
  )

  # PCA result
  X_pca_raw <- hist_res$pca$X_pca
  pc_names  <- colnames(X_pca_raw)

  df_pca <- data.frame(
    pc_x  = X_pca_raw[, 1],
    pc_y  = X_pca_raw[, 2],
    group = as.character(group_vec)
  )

  # frame / breaks / base bins
  xlim <- hist_res$frame$xlim
  ylim <- hist_res$frame$ylim

  x_breaks <- hist_res$breaks$x
  y_breaks <- hist_res$breaks$y
  m_x <- length(x_breaks) - 1
  m_y <- length(y_breaks) - 1
  m   <- m_x * m_y

  base_coord <- do.call(
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

  coord_none_list   <- list()
  coord_add_list    <- list()
  coord_sparse_list <- list()
  samp_none_list    <- list()
  samp_add_list     <- list()
  samp_sparse_list  <- list()

  # hist / sampling per group
  for (g in g_levels) {
    g_res <- hist_res$groups[[g]]
    if (is.null(g_res)) next

    n_g   <- g_res$n
    col_g <- col_map[[g]]

    # Original hist
    prob_none <- g_res$counts / n_g
    coord_none <- base_coord
    coord_none$prob  <- prob_none
    coord_none$group <- g
    coord_none_list[[g]] <- coord_none

    # Original sampling
    if (sampling) {
      coord_tmp <- coord_none
      coord_tmp$mechanism <- "none"
      samp <- sample_from_hist(coord_tmp, k = n_g)
      samp$group <- g
      samp_none_list[[g]] <- samp
    }

    # Add
    if (!is.null(g_res$add)) {
      ca <- g_res$add$coord
      ca$group <- g
      coord_add_list[[g]] <- ca

      if (sampling && !is.null(g_res$add$sample)) {
        sa <- g_res$add$sample
        sa$group <- g
        samp_add_list[[g]] <- sa
      }
    }

    # Sparse
    if (!is.null(g_res$sparse)) {
      cs <- g_res$sparse$coord
      cs$group <- g
      coord_sparse_list[[g]] <- cs

      if (sampling && !is.null(g_res$sparse$sample)) {
        ss <- g_res$sparse$sample
        ss$group <- g
        samp_sparse_list[[g]] <- ss
      }
    }
  }

  coord_none_all   <- if (length(coord_none_list))   dplyr::bind_rows(coord_none_list)   else NULL
  coord_add_all    <- if (length(coord_add_list))    dplyr::bind_rows(coord_add_list)    else NULL
  coord_sparse_all <- if (length(coord_sparse_list)) dplyr::bind_rows(coord_sparse_list) else NULL

  samp_none_all   <- if (length(samp_none_list))   dplyr::bind_rows(samp_none_list)   else NULL
  samp_add_all    <- if (length(samp_add_list))    dplyr::bind_rows(samp_add_list)    else NULL
  samp_sparse_all <- if (length(samp_sparse_list)) dplyr::bind_rows(samp_sparse_list) else NULL

  # Original scatter plot
  p_org_scatter <- ggplot2::ggplot(
    df_pca,
    ggplot2::aes(x = pc_x, y = pc_y, colour = as.character(group))
  ) +
    ggplot2::geom_point(alpha = 0.6, size = 1.8) +
    ggplot2::scale_colour_manual(values = col_map) +
    ggplot2::coord_fixed(xlim = xlim, ylim = ylim) +
    ggplot2::scale_x_continuous(expand = c(0, 0), breaks = pretty(xlim, 5)) +
    ggplot2::scale_y_continuous(expand = c(0, 0), breaks = pretty(ylim, 5)) +
    theme_dp_base_dp() +
    ggplot2::labs(x = pc_names[1], y = pc_names[2])

  p_org_scatter <- add_title_dp(p_org_scatter, "Original Scatter")

  # all-group layout (hist + sampling)
  p_none_hist_all   <- make_hist_all_dp(coord_none_all,   xlim, ylim, col_map, "Original Hist")
  p_add_hist_all    <- if (!is.null(coord_add_all))
    make_hist_all_dp(coord_add_all,    xlim, ylim, col_map, "Add DP Hist")
  else patchwork::plot_spacer()
  p_sparse_hist_all <- if (!is.null(coord_sparse_all))
    make_hist_all_dp(coord_sparse_all, xlim, ylim, col_map, "Sparse DP Hist")
  else patchwork::plot_spacer()

  if (sampling) {
    p_none_sample_all   <- make_sample_all_dp(samp_none_all,   xlim, ylim, col_map, "Original Sampling")
    p_add_sample_all    <- if (!is.null(samp_add_all))
      make_sample_all_dp(samp_add_all,    xlim, ylim, col_map, "Add DP Sampling")
    else patchwork::plot_spacer()
    p_sparse_sample_all <- if (!is.null(samp_sparse_all))
      make_sample_all_dp(samp_sparse_all, xlim, ylim, col_map, "Sparse DP Sampling")
    else patchwork::plot_spacer()

    row1 <- p_org_scatter | p_none_hist_all | p_add_hist_all | p_sparse_hist_all
    row2 <- patchwork::plot_spacer() | p_none_sample_all | p_add_sample_all | p_sparse_sample_all
    p_all <- row1 / row2
  } else {
    p_all <- p_org_scatter | p_none_hist_all | p_add_hist_all | p_sparse_hist_all
  }

  # per-group histogram layout helper
  make_hist_layout <- function(which = c("none", "add", "sparse")) {
    which <- match.arg(which)

    if (which == "none") {
      third_panel <- patchwork::plot_spacer()
      source_list <- coord_none_list
      base_title  <- "Original Hist: "
    } else if (which == "add") {
      third_panel <- p_add_hist_all
      source_list <- coord_add_list
      base_title  <- "Add DP Hist: "
    } else {  # "sparse"
      third_panel <- p_sparse_hist_all
      source_list <- coord_sparse_list
      base_title  <- "Sparse DP Hist: "
    }

    plots_all <- list(
      p_org_scatter,
      p_none_hist_all,
      third_panel
    )

    group_plots <- list()
    for (g in g_levels) {
      if (!g %in% names(source_list)) next

      title_g <- NULL

      if (!.is_color_vec(g_levels)) {
        title_g <- paste0(base_title, g)
      }

      p_g <- make_hist_single_dp(
        source_list[[g]],
        xlim = xlim,
        ylim = ylim,
        col  = col_map[[g]],
        title = title_g
      )
      group_plots[[g]] <- p_g
      plots_all[[length(plots_all) + 1L]] <- p_g
    }

    all_layout <- patchwork::wrap_plots(plots_all, ncol = 3)

    list(
      all    = all_layout,
      groups = group_plots
    )
  }

  none_hist_layout   <- make_hist_layout("none")
  add_hist_layout    <- make_hist_layout("add")
  sparse_hist_layout <- make_hist_layout("sparse")

  make_sampling_layout <- function(which = c("none", "add", "sparse")) {
    which <- match.arg(which)
    if (!sampling) return(list(all = NULL, groups = list()))

    if (which == "none") {
      third_panel <- patchwork::plot_spacer()
      source_list <- samp_none_list
      base_title  <- "Original Sampling: "
    } else if (which == "add") {
      third_panel <- make_sample_all_dp(samp_add_all, xlim, ylim, col_map, "Add DP Sampling")
      source_list <- samp_add_list
      base_title  <- "Add DP Sampling: "
    } else {  # "sparse"
      third_panel <- make_sample_all_dp(samp_sparse_all, xlim, ylim, col_map, "Sparse DP Sampling")
      source_list <- samp_sparse_list
      base_title  <- "Sparse DP Sampling: "
    }

    plots_all <- list(
      p_org_scatter,
      make_sample_all_dp(samp_none_all, xlim, ylim, col_map, "Original Sampling"),
      third_panel
    )

    g_plots <- list()
    for (g in g_levels) {
      if (!g %in% names(source_list)) next

      title_g <- NULL
      if (!.is_color_vec(g_levels)) {
        title_g <- paste0(base_title, g)
      }

      p_g <- make_sample_single_dp(
        source_list[[g]],
        xlim = xlim,
        ylim = ylim,
        col  = col_map[[g]],
        title = title_g
      )
      g_plots[[g]] <- p_g
      plots_all[[length(plots_all) + 1L]] <- p_g
    }

    all_layout <- patchwork::wrap_plots(plots_all, ncol = 3)

    list(
      all    = all_layout,
      groups = g_plots
    )
  }

  none_sampling_layout   <- make_sampling_layout("none")
  add_sampling_layout    <- make_sampling_layout("add")
  sparse_sampling_layout <- make_sampling_layout("sparse")

  # result
  hist_res$plot <- list(
    all = p_all,

    none_hist   = c(list(all = none_hist_layout$all),   none_hist_layout$groups),
    add_hist    = c(list(all = add_hist_layout$all),    add_hist_layout$groups),
    sparse_hist = c(list(all = sparse_hist_layout$all), sparse_hist_layout$groups),

    none_sampling   = c(list(all = none_sampling_layout$all),   none_sampling_layout$groups),
    add_sampling    = c(list(all = add_sampling_layout$all),    add_sampling_layout$groups),
    sparse_sampling = c(list(all = sparse_sampling_layout$all), sparse_sampling_layout$groups)
  )

  hist_res$group_colors <- col_map
  hist_res
}

