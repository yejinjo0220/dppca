#' DP PCA Score Plot with 2D Histograms
#'
#' @description
#' Creates a visualization that combines the original PCA score plot and
#' (optionally differentially private) 2D histograms based on the first two
#' principal components. Depending on the selected mechanism, the function
#' displays the empirical histogram, a Gaussian-additive DP histogram, and/or
#' a sparse Laplace-thresholded DP histogram, together with optional synthetic
#' samples drawn from each histogram.
#'
#' @param X A numeric matrix or an object coercible to a matrix with
#'   \eqn{n} rows (observations) and \eqn{p} columns (variables).
#'   PCA is performed on \code{X}, and the first two PCs are used for plotting.
#' @param center Logical; should the variables be centered before PCA?
#'   Passed to \code{scale()}.
#' @param scale. Logical; should the variables be scaled to unit variance?
#'   Passed to \code{scale()}.
#' @param dp_pca_flag Logical; if \code{TRUE}, apply differentially private PCA
#'   using \code{dp_pca()} with part of the total privacy budget
#'   \code{eps_total}, \code{delta_total}. If \code{FALSE}, use non-DP PCA.
#' @param cpp.option Logical; if \code{TRUE}, allow C++ implementations inside
#'   DP mechanisms (e.g., \code{mech_tau_sph()}) when available.
#' @param axes Integer vector indicating which principal components to return
#'   (e.g., \code{c(1, 2)} for the first two PCs). The maximum index in
#'   \code{axes} must not exceed the number of variables \code{p}.
#'
#' @param eps_total Total privacy budget \eqn{\varepsilon_{\text{total}}} to be
#'   split among DP-PCA (optional), DP frame (\code{dp_frame()}), and DP
#'   histogram mechanisms (inside \code{dp_hist()}).
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
#'   \code{dp_frame()} via \code{dp_hist()}. If \code{NULL}, extreme order
#'   statistics \eqn{1/n} and \eqn{(n-1)/n} are used.
#'
#' @param m_x Optional integer specifying the number of bins along the first PC
#'   axis. If \code{NULL}, it is chosen via \code{number_bins()} inside
#'   \code{dp_hist()}.
#' @param m_y Optional integer specifying the number of bins along the second PC
#'   axis. If \code{NULL}, it is chosen via \code{number_bins()}.
#' @param bin_method Character string specifying the heuristic used by
#'   \code{number_bins()} when \code{m_x} and/or \code{m_y} are not supplied.
#'   One of \code{"J"} (Jing) or \code{"W"} (Wasserman).
#'
#' @param mechanism Character string indicating which histogram mechanism(s) to
#'   compute and display:
#'   \itemize{
#'     \item \code{"none"}: non-DP empirical histogram.
#'     \item \code{"add"}: Gaussian additive-noise DP histogram.
#'     \item \code{"sparse"}: sparse Laplace-thresholded DP histogram.
#'     \item \code{"all"}: compute and, if possible, display all three versions.
#'   }
#'
#' @param sampling Logical; if \code{TRUE}, additionally generate and display
#'   synthetic samples from each available histogram version using
#'   \code{sample_from_hist()} with \code{k = n}.
#' @param color A character string specifying the base color used for points and
#'   histogram tiles (default is a purple tone, \code{"#6A5ACD"}).
#'
#' @details
#' This function is a visualization wrapper around \code{dp_hist()}. Internally,
#' it first calls \code{dp_hist()} to compute (DP) PCA scores, a differentially
#' private frame, and one or more histogram mechanisms on the 2D scores. It then
#' builds a panel of ggplot2 plots showing:
#' \itemize{
#'   \item the original PCA scatter plot;
#'   \item the empirical (non-DP) histogram on the 2D scores;
#'   \item the Gaussian-additive DP histogram (if requested);
#'   \item the sparse DP histogram (if requested);
#'   \item optionally, synthetic samples drawn from each histogram version
#'         (if \code{sampling = TRUE}).
#' }
#'
#' Internally, helper functions such as \code{theme_dp_base_dp()},
#' \code{add_title_dp()}, \code{make_hist_plot_dp()}, and
#' \code{make_sample_plot_dp()} are used to construct a consistent visual layout
#' based on ggplot2 and patchwork.
#'
#' @return
#' A list equal to the output of \code{dp_hist()}, augmented with an additional
#' component \code{plot}, which is itself a list containing:
#' \itemize{
#'   \item \code{all}: the full patchwork object combining all panels.
#'   \item \code{org_scatter}: the original PCA scatter plot.
#'   \item \code{org_hist}: the empirical histogram plot.
#'   \item \code{org_sampling}: the empirical sampling plot (if available).
#'   \item \code{add_hist}: the Gaussian-additive DP histogram plot (if any).
#'   \item \code{add_sampling}: the Gaussian-additive sampling plot (if any).
#'   \item \code{sparse_hist}: the sparse DP histogram plot (if any).
#'   \item \code{sparse_sampling}: the sparse sampling plot (if any).
#' }
#'
#' @seealso
#'   \code{\link{dp_hist}},
#'   \code{\link{dp_score_plot_group}},
#'   \code{\link{dp_pca}},
#'   \code{\link{dp_frame}},
#'   \code{\link{number_bins}},
#'   \code{\link{sample_from_hist}}
#'
#' @importFrom ggplot2 ggplot aes geom_point coord_fixed scale_x_continuous
#'   scale_y_continuous labs
#' @importFrom patchwork plot_spacer
#' @export
dp_score_plot <- function(
    X,
    # PCA option
    center = TRUE,
    scale. = FALSE,
    dp_pca_flag = FALSE,
    cpp.option = TRUE,
    axes = c(1, 2),

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
    bin_method = c("J", "W", "none"),

    mechanism = c("all", "none", "add", "sparse"),

    # Sampling
    sampling  = FALSE,

    # color option
    color = "#6A5ACD"
) {
  bin_method <- match.arg(bin_method)
  mechanism  <- match.arg(mechanism)

  # bin method
  if (bin_method == "none") {
    if (is.null(m_x) || is.null(m_y)) {
      stop("When bin_method = 'none', you must supply both m_x and m_y.")
    }
    bin_method_for_hist <- "J"
    m_x_for_hist <- m_x
    m_y_for_hist <- m_y
  } else {
    bin_method_for_hist <- bin_method
    m_x_for_hist <- m_x
    m_y_for_hist <- m_y
  }

  # dp_hist
  hist_res <- .dp_hist_retry(
    X            = X,
    center       = center,
    scale.       = scale.,
    dp_pca_flag  = dp_pca_flag,
    cpp.option   = cpp.option,
    axes         = axes,
    eps_total    = eps_total,
    delta_total  = delta_total,
    eps_ratio    = eps_ratio,
    delta_ratio  = delta_ratio,
    inflate      = inflate,
    q_frame      = q_frame,
    m_x          = m_x_for_hist,
    m_y          = m_y_for_hist,
    bin_method   = bin_method_for_hist,
    mechanism    = mechanism,
    sampling     = sampling,
    verbose      = TRUE # numbers of trying until not errors
  )

  # PCA result
  X_pca_raw <- hist_res$pca$X_pca
  pc_names  <- colnames(X_pca_raw)

  X_pca <- as.data.frame(X_pca_raw)
  colnames(X_pca) <- c("pc_x", "pc_y")

  xlim <- hist_res$frame$xlim
  ylim <- hist_res$frame$ylim

  hist_none   <- hist_res$none$coord
  hist_add    <- if (!is.null(hist_res$add))    hist_res$add$coord    else NULL
  hist_sparse <- if (!is.null(hist_res$sparse)) hist_res$sparse$coord else NULL

  # Original scatter plot
  p_org_scatter <- ggplot2::ggplot(X_pca, ggplot2::aes(x = pc_x, y = pc_y)) +
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
    ggplot2::labs(x = pc_names[1], y = pc_names[2])

  p_org_scatter <- add_title_dp(p_org_scatter, "Original Scatter")

  # Histogram plot
  p_none_hist   <- make_hist_plot_dp(hist_none,   xlim, ylim, color, "Original Hist")
  p_add_hist    <- make_hist_plot_dp(hist_add,    xlim, ylim, color, "Add DP Hist")
  p_sparse_hist <- make_hist_plot_dp(hist_sparse, xlim, ylim, color, "Sparse DP Hist")

  # Sampling plot
  p_none_sample <- if (sampling && !is.null(hist_res$none$sample)) {
    make_sample_plot_dp(hist_res$none$sample, xlim, ylim, color, "Original Sampling")
  } else NULL

  p_add_sample <- if (sampling && !is.null(hist_res$add) && !is.null(hist_res$add$sample)) {
    make_sample_plot_dp(hist_res$add$sample, xlim, ylim, color, "Add DP Sampling")
  } else NULL

  p_sparse_sample <- if (sampling && !is.null(hist_res$sparse) && !is.null(hist_res$sparse$sample)) {
    make_sample_plot_dp(hist_res$sparse$sample, xlim, ylim, color, "Sparse DP Sampling")
  } else NULL

  # Layout
  if (sampling) {
    p_add_hist_use    <- if (!is.null(p_add_hist))    p_add_hist    else patchwork::plot_spacer()
    p_sparse_hist_use <- if (!is.null(p_sparse_hist)) p_sparse_hist else patchwork::plot_spacer()

    row1 <- p_org_scatter | p_none_hist | p_add_hist_use | p_sparse_hist_use

    p_none_sample_use   <- if (!is.null(p_none_sample))   p_none_sample   else patchwork::plot_spacer()
    p_add_sample_use    <- if (!is.null(p_add_sample))    p_add_sample    else patchwork::plot_spacer()
    p_sparse_sample_use <- if (!is.null(p_sparse_sample)) p_sparse_sample else patchwork::plot_spacer()

    row2 <- patchwork::plot_spacer() | p_none_sample_use | p_add_sample_use | p_sparse_sample_use

    p_all <- row1 / row2

  } else {
    # sampling = FALSE
    p_add_hist_use    <- if (!is.null(p_add_hist))    p_add_hist    else patchwork::plot_spacer()
    p_sparse_hist_use <- if (!is.null(p_sparse_hist)) p_sparse_hist else patchwork::plot_spacer()

    p_all <- p_org_scatter | p_none_hist | p_add_hist_use | p_sparse_hist_use
  }

  # result
  hist_res$plot <- list(
    all             = p_all,
    org_scatter     = p_org_scatter,
    org_hist        = p_none_hist,
    org_sampling    = p_none_sample,
    add_hist        = p_add_hist,
    add_sampling    = p_add_sample,
    sparse_hist     = p_sparse_hist,
    sparse_sampling = p_sparse_sample
  )

  hist_res
}



