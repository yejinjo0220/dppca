# ============================================================
# score.R
# Functions for differentially private score histograms
# ============================================================

#' Differentially private score histograms
#'
#' This function computes two-dimensional principal component scores and returns
#' differentially private histogram estimates on the score space. It returns the
#' score coordinates, the plotting frame, the non-private histogram, and the
#' requested private histogram estimates.
#'
#' @param X A numeric matrix or data frame. Rows correspond to observations and
#'  columns correspond to variables.
#' @param eps Positive number defining the `epsilon` privacy parameter supplied
#'   to each requested score-histogram procedure. When multiple histogram
#'   methods are requested, the same value is used for each method for
#'   comparison; it is not divided across methods.
#' @param delta Number in `(0, 1)` defining the `delta` privacy parameter
#'   supplied to each requested score-histogram procedure. When multiple
#'   histogram methods are requested, the same value is used for each method
#'   for comparison; it is not divided across methods.
#' @param bins Integer vector of length 2 defining the number of histogram bins
#'   along the first and second score axes, respectively.
#' @param method Character vector specifying which private histogram methods to
#'   compute. Use `"add"` for the additive Gaussian histogram and `"sparse"` for
#'   the sparse thresholded histogram. The default is `c("add", "sparse")`.
#' @param center A logical value indicating whether to center the columns of `X`
#'   before computing principal component directions. The default is `TRUE`.
#' @param standardize A logical value indicating whether to scale the columns of
#'   `X` by their sample standard deviations after optional centering. The
#'   default is `FALSE`.
#' @param g_dppca A logical value indicating whether to use private principal
#'   component directions. The default is `FALSE`. See [dp_pc_dir()] for details.
#' @param cpp.option A logical value passed to [dp_pc_dir()] when
#'   `g_dppca = TRUE`. The default is `FALSE`.
#' @param axes Integer vector of length 2 specifying the principal components
#'   used to construct the score coordinates. The default is `c(1, 2)`.
#'
#' @details
#' Let \eqn{v_a} and \eqn{v_b} be the principal component directions selected
#' by `axes = c(a, b)` for some \eqn{1 \le a < b \le ncol(X)}.
#' After preprocessing, the score point for \eqn{i}th observation
#' is \eqn{s_i = (x_i^\top v_a, x_i^\top v_b)}. A non-private score
#' plot would display the points \eqn{s_1, \ldots, s_n} directly. This function
#' instead summarizes their empirical distribution by a two-dimensional histogram
#' and releases private versions of the histogram for the visualization.
#'
#' The plotting frame is constructed privately from the score coordinates using
#' the pure-DP unbounded quantile mechanisms of
#' \insertCite{durfee2023unbounded;textual}{dppca}. Its center is estimated by
#' coordinate-wise private medians. For each coordinate, the private 0.995
#' quantile of the absolute deviations from its private median is then estimated
#' and inflated by a fixed factor. The two independently estimated radii form a
#' rectangular plotting frame.
#'
#' The private histogram is computed on the rectangular grid defined by the
#' private frame and the bin counts in `bins`. Under
#' row-level adjacency, changing one observation can increase one bin count by
#' one and decrease another by one, giving \eqn{\ell_1} sensitivity at most
#' \eqn{2} and \eqn{\ell_2} sensitivity at most \eqn{\sqrt{2}} for the count
#' vector.
#'
#' Two private histogram mechanisms are supported:
#' \itemize{
#'   \item `"add"` constructs an additive differentially private histogram by
#'   adding Gaussian noise to all bin counts, clipping negative noisy counts to
#'   zero, and normalizing the result. This additive-noise approach is commonly
#'   used for private histograms; see
#'   \insertCite{wasserman2010statistical;textual}{dppca}.
#'
#'   \item `"sparse"` constructs a sparse differentially private histogram for
#'   settings where many bins are empty. It perturbs only nonzero empirical bin
#'   proportions and keeps bins whose noisy values exceed a stability threshold,
#'   following the stability-based private histogram idea of
#'   \insertCite{karwa2017finite;textual}{dppca}.
#' }
#'
#' The privacy parameters are allocated across the privacy-consuming steps. If
#' `g_dppca = FALSE`, 20 percent of `eps` is used for the private frame and 80
#' percent for the private histogram; all of `delta` is used by the histogram.
#' If `g_dppca = TRUE`, `eps` is allocated in proportions 0.2, 0.2, and 0.6 to
#' private direction estimation, private frame construction, and private
#' histogram release, respectively. The corresponding `delta` proportions are
#' 0.2 for private directions and 0.8 for the histogram. Frame construction is
#' pure DP and divides `eps_frame` equally among its two private medians and two
#' private radius estimates.
#'
#' When multiple histogram methods are requested, the histogram privacy budget
#' is not divided across `"add"` and `"sparse"`. Instead, each requested method
#' receives the same histogram portion of `eps` and `delta` so that the methods
#' can be compared under the same privacy setting. If outputs from multiple
#' private histogram methods are released together, the privacy cost of the
#' joint release must be accounted for by composition.
#'
#' The returned score coordinates and the `nonprivate` histogram are included as
#' non-private references and are not themselves differentially private releases.
#'
#' For a detailed procedure and mathematical formulations,
#' refer \url{https://yejinjo0220.github.io/dppca/articles/dp_score}.
#'
#' @return A named list with components `score` and `frame`, followed by
#' histogram results. The `nonprivate` component contains the non-private
#' empirical histogram. Each requested private method is returned as an
#' additional component (`add` and/or `sparse`). Methods that are not requested
#' are omitted from the returned list.
#'
#' @seealso
#' [dp_score_plot()] for plotting the output of this function.
#' [dp_score_group()] and [dp_score_plot_group()] for group-wise score
#' histograms.
#' [dp_pc_dir()] for private principal component direction estimation.
#'
#' @references
#' \insertRef{dwork2014algorithmic}{dppca}
#'
#' \insertRef{durfee2023unbounded}{dppca}
#'
#' \insertRef{wasserman2010statistical}{dppca}
#'
#' \insertRef{karwa2017finite}{dppca}
#'
#' \insertRef{kim2025robustdppca}{dppca}
#'
#' @examples
#' data(gau, package = "dppca")
#'
#' # Use a small subset to keep the example fast.
#' X <- gau[1:300, ]
#'
#' # Compute private two-dimensional PCA scores using the additive histogram method.
#' set.seed(123)
#' score_gau <- dp_score(
#'   X,
#'   eps = 2,
#'   delta = 1e-3,
#'   method = "add",
#'   bins = c(10, 10)
#' )
#'
#' head(score_gau$score)
#' head(score_gau$add)
#'
#' @importFrom VGAM rlaplace
#' @export
dp_score <- function(
    X,
    eps,
    delta,
    bins,
    method = c("add", "sparse"),
    center = TRUE,
    standardize = FALSE,
    g_dppca = FALSE,
    cpp.option = FALSE,
    axes = c(1, 2)
) {
  X <- validate_score_matrix(X)
  validate_score_common(
    X = X,
    eps = eps,
    delta = delta,
    bins = bins,
    center = center,
    standardize = standardize,
    g_dppca = g_dppca,
    cpp.option = cpp.option,
    axes = axes
  )

  axes <- as.integer(axes)
  bins <- as.integer(bins)
  method <- unique(
    match.arg(
      method,
      choices = c("add", "sparse"),
      several.ok = TRUE
    )
  )

  budget <- split_score_privacy_budget(
    eps = eps,
    delta = delta,
    g_dppca = g_dppca
  )

  score_res <- compute_score_coordinates(
    X = X,
    axes = axes,
    center = center,
    standardize = standardize,
    g_dppca = g_dppca,
    cpp.option = cpp.option,
    eps_pc = budget$eps_pc,
    delta_pc = budget$delta_pc
  )

  frame_out <- dp_frame(
    X = score_res$score,
    eps_frame = budget$eps_frame,
    inflate = 0.20
  )

  hist <- score_histograms(
    X_score = score_res$score,
    xlim = frame_out$xlim,
    ylim = frame_out$ylim,
    bins = bins,
    eps_hist = budget$eps_hist,
    delta_hist = budget$delta_hist,
    method = method
  )

  c(
    list(
      score = score_res$score,
      frame = frame_out
    ),
    hist
  )
}

#' Plot differentially private score histograms or sampled score points
#'
#' This function computes two-dimensional principal component score histograms
#' and visualizes the results. Private histogram releases can be shown directly
#' as histogram panels or post-processed into synthetic scatter plots sampled
#' from the released histograms.
#'
#' @inheritParams dp_score
#' @param private_plot Character vector controlling how requested histogram
#'   results are visualized. `"histogram"` displays the histogram panels,
#'   while `"sample"` displays synthetic score points sampled from the
#'   histograms. Both can be requested with
#'   `c("histogram", "sample")`. The default is
#'   `c("histogram", "sample")`.
#' @param sampling_control Optional control list created by
#'   `sampling_control()`. It is used only when `"sample"` is included in
#'   `private_plot`. If `NULL`, the default sampling settings are used.
#' @param plot_control Optional control list created by `score_plot_control()`.
#'   If `NULL`, default score-plot settings are used.
#'
#' @details
#' Histogram plots and sampled-score plots are arranged dynamically.
#' When only `"histogram"` is requested, the original scatter plot is shown
#' in the first column of the histogram row. When `"sample"` is requested,
#' the original scatter plot is instead shown in the first column of the
#' first sample row. If both plot types are requested, the histogram row uses
#' a blank first panel so that the sample rows line up beneath it.
#'
#' Synthetic score sampling is applied not only to the private histograms,
#' but also to the non-private histogram, so that non-private and private
#' sampled-score visualizations can be compared side by side.
#'
#' @return A named list containing:
#' \describe{
#'   \item{score}{The output of [dp_score()].}
#'   \item{sample}{Present only when `"sample"` is requested. A nested list
#'   containing sampled coordinates for `nonprivate` and for each requested
#'   private histogram method (`add` and/or `sparse`). Each leaf is a data
#'   frame with columns `pc_x` and `pc_y`, indexed by sampling method
#'   (`center` and/or `uniform`).}
#'   \item{plot}{A list containing the original scatter plot, the non-private
#'   histogram plot, requested private histogram plots, sampled-score plots,
#'   and the combined patchwork layout in `plot$all`. Sampled-score plots are
#'   stored under `plot$sample`.}
#' }
#'
#' @seealso
#' [dp_score()] for computing score histograms without plotting.
#' [dp_score_plot_group()] for group-wise score histogram plots.
#'
#' @examples
#' data(gau, package = "dppca")
#'
#' X <- gau[1:300, ]
#'
#' # By default, both histogram and sampled panels are displayed.
#' set.seed(123)
#' score_plot <- dp_score_plot(
#'   X,
#'   eps = 3,
#'   delta = 1e-3,
#'   bins = c(8, 8),
#'   method = "add"
#' )
#' score_plot$plot$all
#'
#' # Sample-only display with two sampling methods.
#' set.seed(123)
#' sample_plot <- dp_score_plot(
#'   X,
#'   eps = 3,
#'   delta = 1e-3,
#'   bins = c(8, 8),
#'   method = "add",
#'   private_plot = "sample",
#'   sampling_control = sampling_control(
#'     method = c("center", "uniform")
#'   )
#' )
#' sample_plot$plot$all
#'
#' @export
dp_score_plot <- function(
    X,
    eps,
    delta,
    bins,
    method = c("add", "sparse"),
    private_plot = c("histogram", "sample"),
    sampling_control = NULL,
    center = TRUE,
    standardize = FALSE,
    g_dppca = FALSE,
    cpp.option = FALSE,
    axes = c(1, 2),
    plot_control = NULL
) {
  method <- unique(
    match.arg(
      method,
      choices = c("add", "sparse"),
      several.ok = TRUE
    )
  )

  private_plot <- unique(
    match.arg(
      private_plot,
      choices = c("histogram", "sample"),
      several.ok = TRUE
    )
  )

  pctrl <- .merge_score_plot_control(plot_control)

  score_res <- dp_score(
    X = X,
    eps = eps,
    delta = delta,
    bins = bins,
    center = center,
    standardize = standardize,
    g_dppca = g_dppca,
    cpp.option = cpp.option,
    axes = axes,
    method = method
  )

  X_score <- as.data.frame(score_res$score)
  colnames(X_score) <- c("pc_x", "pc_y")

  xlim <- score_res$frame$xlim
  ylim <- score_res$frame$ylim
  pc_names <- colnames(score_res$score)

  xlab <- if (is.null(pctrl$xlab)) pc_names[1] else pctrl$xlab
  ylab <- if (is.null(pctrl$ylab)) pc_names[2] else pctrl$ylab

  p_scatter <- ggplot2::ggplot(
    X_score,
    ggplot2::aes(x = .data$pc_x, y = .data$pc_y)
  ) +
    ggplot2::geom_point(
      alpha = pctrl$scatter_alpha,
      size = pctrl$scatter_size,
      color = pctrl$color
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
    theme_dp_base(base_size = pctrl$base_size) +
    ggplot2::labs(x = xlab, y = ylab)

  p_scatter <- add_title_dp(
    p_scatter,
    pctrl$scatter_title,
    title_size = pctrl$title_size
  )

  p_nonprivate_hist <- make_hist_plot_dp(
    score_res$nonprivate,
    xlim = xlim,
    ylim = ylim,
    color = pctrl$color,
    title = pctrl$nonprivate_title,
    xlab = xlab,
    ylab = ylab,
    alpha_range = pctrl$hist_alpha_range,
    base_size = pctrl$base_size,
    title_size = pctrl$title_size
  )

  plot_out <- list(
    scatter = p_scatter,
    nonprivate = p_nonprivate_hist
  )

  hist_plot_out <- list(nonprivate = p_nonprivate_hist)

  if ("add" %in% method) {
    hist_plot_out$add <- make_hist_plot_dp(
      score_res$add,
      xlim = xlim,
      ylim = ylim,
      color = pctrl$color,
      title = pctrl$add_title,
      xlab = xlab,
      ylab = ylab,
      alpha_range = pctrl$hist_alpha_range,
      base_size = pctrl$base_size,
      title_size = pctrl$title_size
    )
    plot_out$add <- hist_plot_out$add
  }

  if ("sparse" %in% method) {
    hist_plot_out$sparse <- make_hist_plot_dp(
      score_res$sparse,
      xlim = xlim,
      ylim = ylim,
      color = pctrl$color,
      title = pctrl$sparse_title,
      xlab = xlab,
      ylab = ylab,
      alpha_range = pctrl$hist_alpha_range,
      base_size = pctrl$base_size,
      title_size = pctrl$title_size
    )
    plot_out$sparse <- hist_plot_out$sparse
  }

  sample_out <- NULL
  sample_plot_out <- NULL
  sctrl <- NULL

  if ("sample" %in% private_plot) {
    sctrl <- .resolve_sampling_control(
      control = sampling_control,
      n = nrow(score_res$score)
    )

    sample_out <- list()
    sample_plot_out <- list()

    sampled_np <- sample_private_score_histogram(
      hist_df = score_res$nonprivate,
      sample_size = sctrl$sample_size,
      sample_method = sctrl$method,
      bandwidth_scale = sctrl$bandwidth_scale
    )

    sample_out$nonprivate <- sampled_np
    sample_plot_out$nonprivate <- list()

    for (sm in names(sampled_np)) {
      sm_label <- if (sm == "center") "Center" else "Uniform"
      sample_plot_out$nonprivate[[sm]] <- make_sample_plot_dp(
        sample_df = sampled_np[[sm]],
        xlim = xlim,
        ylim = ylim,
        color = pctrl$color,
        title = paste0("Non-private Sample (", sm_label, ")"),
        xlab = xlab,
        ylab = ylab,
        point_alpha = pctrl$scatter_alpha,
        point_size = pctrl$scatter_size,
        base_size = pctrl$base_size,
        title_size = pctrl$title_size
      )
    }

    for (hist_method in method) {
      sampled <- sample_private_score_histogram(
        hist_df = score_res[[hist_method]],
        sample_size = sctrl$sample_size,
        sample_method = sctrl$method,
        bandwidth_scale = sctrl$bandwidth_scale
      )

      sample_out[[hist_method]] <- sampled
      sample_plot_out[[hist_method]] <- list()

      hist_label <- if (hist_method == "add") "Add" else "Sparse"

      for (sm in names(sampled)) {
        sm_label <- if (sm == "center") "Center" else "Uniform"
        sample_plot_out[[hist_method]][[sm]] <- make_sample_plot_dp(
          sample_df = sampled[[sm]],
          xlim = xlim,
          ylim = ylim,
          color = pctrl$color,
          title = paste0(hist_label, " Sample (", sm_label, ")"),
          xlab = xlab,
          ylab = ylab,
          point_alpha = pctrl$scatter_alpha,
          point_size = pctrl$scatter_size,
          base_size = pctrl$base_size,
          title_size = pctrl$title_size
        )
      }
    }

    plot_out$sample <- sample_plot_out
  }

  # Build one flat patchwork grid rather than nesting row-wise patchworks.
  # Using a single grid is important here: nested patchworks can assign
  # different effective widths to otherwise corresponding columns, especially
  # when fixed-aspect plots and spacers are mixed. A single row-major grid keeps
  # every column aligned across histogram and sample rows.
  n_cols <- 2 + length(method)
  layout_panels <- list()

  if ("histogram" %in% private_plot) {
    if ("sample" %in% private_plot) {
      layout_panels[[length(layout_panels) + 1L]] <-
        patchwork::plot_spacer()
    } else {
      layout_panels[[length(layout_panels) + 1L]] <- p_scatter
    }

    layout_panels[[length(layout_panels) + 1L]] <-
      hist_plot_out$nonprivate

    for (m in method) {
      layout_panels[[length(layout_panels) + 1L]] <-
        hist_plot_out[[m]]
    }
  }

  if ("sample" %in% private_plot) {
    smethods <- sctrl$method

    for (i in seq_along(smethods)) {
      sm <- smethods[i]

      if (i == 1L) {
        layout_panels[[length(layout_panels) + 1L]] <- p_scatter
      } else {
        layout_panels[[length(layout_panels) + 1L]] <-
          patchwork::plot_spacer()
      }

      layout_panels[[length(layout_panels) + 1L]] <-
        sample_plot_out$nonprivate[[sm]]

      for (m in method) {
        layout_panels[[length(layout_panels) + 1L]] <-
          sample_plot_out[[m]][[sm]]
      }
    }
  }

  plot_out$all <- patchwork::wrap_plots(
    layout_panels,
    ncol = n_cols,
    byrow = TRUE
  )

  out <- list(score = score_res)

  if (!is.null(sample_out)) {
    out$sample <- sample_out
  }

  out$plot <- plot_out
  out
}


#' Group-wise differentially private score histograms
#'
#' This function computes two-dimensional principal component scores and releases
#' group-wise differentially private histograms on a common score frame and grid.
#' It is useful when observations have group labels and the low-dimensional score
#' distribution should be compared across groups.
#'
#' @inheritParams dp_score
#' @param X A matrix or data frame where rows correspond to observations
#'   and columns correspond to variables.
#'   `X` can additionally include a named column representing the group label
#'   for each observation.
#' @param group Group labels. This can be a vector of length `nrow(X)` or a
#'   single column name in `X`. If a column name is supplied, that column is
#'   used as the group label and removed from the feature matrix.
#'
#' @details
#' The score directions, plotting frame, and histogram grid are shared across all
#' groups. For each group \eqn{g}, the group-specific count in bin \eqn{B_k} is
#' \eqn{c_k^{(g)} = \sum_i 1\{s_i \in B_k, g_i = g\}}. Private histograms are
#' then computed separately for each group on the common grid. Because the groups
#' form a partition of the rows, group-wise histograms for the same histogram
#' method use the same histogram privacy parameters across groups by parallel
#' composition.
#'
#' When both `"add"` and `"sparse"` are requested, the histogram privacy budget
#' is not divided across the two methods. Each method receives the same histogram
#' portion of `eps` and `delta` for method comparison. Releasing outputs from
#' multiple private histogram methods together requires composition across
#' methods.
#'
#' @return A named list with components `score`, `frame`, and `groups`. Each
#' group entry contains `n`, `nonprivate`, and the requested private histogram
#' method components (`add` and/or `sparse`). Unrequested methods are omitted.
#'
#' @seealso
#' [dp_score_plot_group()] for plotting group-wise score histograms.
#' [dp_score()] for pooled score histograms.
#'
#' @examples
#' data(gau_g, package = "dppca")
#'
#' # Compute private grouped PCA scores.
#' set.seed(123)
#' score_gau_g <- dp_score_group(
#'   gau_g,
#'   group = "group",
#'   eps = 3,
#'   delta = 1e-3,
#'   bins = c(8, 8)
#' )
#'
#' head(score_gau_g$score)
#' head(score_gau_g$groups$group1$add)
#'
#' @importFrom VGAM rlaplace
#' @export
dp_score_group <- function(
    X,
    group,
    eps,
    delta,
    bins,
    method = c("add", "sparse"),
    center = TRUE,
    standardize = FALSE,
    g_dppca = FALSE,
    cpp.option = FALSE,
    axes = c(1, 2)
) {
  method <- unique(
    match.arg(
      method,
      choices = c("add", "sparse"),
      several.ok = TRUE
    )
  )

  group_data <- split_group_input(X, group)
  X_mat <- validate_score_matrix(group_data$X)
  group_vec <- group_data$group

  if (length(group_vec) != nrow(X_mat)) {
    stop("`group` must have length equal to `nrow(X)`.", call. = FALSE)
  }
  if (anyNA(group_vec)) {
    stop("`group` must not contain missing values.", call. = FALSE)
  }

  validate_score_common(
    X = X_mat,
    eps = eps,
    delta = delta,
    bins = bins,
    center = center,
    standardize = standardize,
    g_dppca = g_dppca,
    cpp.option = cpp.option,
    axes = axes
  )

  axes <- as.integer(axes)
  bins <- as.integer(bins)
  m_x <- bins[1]
  m_y <- bins[2]
  g_levels <- as.character(unique(group_vec))

  budget <- split_score_privacy_budget(
    eps = eps,
    delta = delta,
    g_dppca = g_dppca
  )

  score_res <- compute_score_coordinates(
    X = X_mat,
    axes = axes,
    center = center,
    standardize = standardize,
    g_dppca = g_dppca,
    cpp.option = cpp.option,
    eps_pc = budget$eps_pc,
    delta_pc = budget$delta_pc
  )

  frame_out <- dp_frame(
    X = score_res$score,
    eps_frame = budget$eps_frame,
    inflate = 0.20
  )

  grid <- score_histogram_grid(
    xlim = frame_out$xlim,
    ylim = frame_out$ylim,
    m_x = m_x,
    m_y = m_y
  )

  # The same histogram privacy budget is supplied separately to each
  # requested method. It is not divided across add and sparse.
  eps_hist_method <- budget$eps_hist
  delta_hist_method <- budget$delta_hist

  groups_out <- list()

  for (g in g_levels) {
    idx_g <- which(as.character(group_vec) == g)
    Xg <- score_res$score[idx_g, , drop = FALSE]
    n_g <- nrow(Xg)

    hist <- score_histograms_from_grid(
      X_score = Xg,
      grid = grid,
      eps_hist_method = eps_hist_method,
      delta_hist_method = delta_hist_method,
      method = method,
      group_name = g
    )

    groups_out[[g]] <- c(
      list(n = n_g),
      hist
    )
  }

  list(
    score = score_res$score,
    frame = frame_out,
    groups = groups_out
  )
}

#' Plot group-wise differentially private score histograms or samples
#'
#' This function computes group-wise two-dimensional principal component score
#' histograms and visualizes the groups together using a common color mapping.
#' Private histogram releases can be displayed directly or post-processed into
#' grouped synthetic scatter plots.
#'
#' @inheritParams dp_score_group
#' @param private_plot Character vector controlling how requested histogram
#'   results are visualized. `"histogram"` displays the overlaid group-wise
#'   histogram panels, while `"sample"` displays grouped synthetic score points
#'   sampled from the group-wise histograms. Both can be requested with
#'   `c("histogram", "sample")`. The default is
#'   `c("histogram", "sample")`.
#' @param sampling_control Optional control list created by
#'   `sampling_control()`. It is used only when `"sample"` is included in
#'   `private_plot`. If `NULL`, the default sampling settings are used.
#' @param plot_control Optional control list created by `score_plot_control()`.
#'   If `NULL`, default grouped score-plot settings are used.
#'
#' @details
#' All groups are shown together in each panel using group-specific colors.
#' Separate per-group plot layouts are not created.
#'
#' For sampled-score panels, sampling is carried out separately from each
#' group's histogram and the sampled coordinates are then combined with the
#' corresponding group label. The same procedure is applied to the non-private
#' group histograms and to each requested private histogram method.
#'
#' `sampling_control(sample_size = NULL)` uses a total synthetic sample size
#' equal to the number of input observations. The total is allocated across
#' groups in proportion to their observed sizes, which reproduces the original
#' group sizes when the default total is used. If a positive integer
#' `sample_size` is supplied, it is interpreted as the total synthetic sample
#' size across all groups and is allocated proportionally using a
#' largest-remainder rule.
#'
#' The combined plot follows the same dynamic layout as [dp_score_plot()].
#' When only `"histogram"` is requested, the original grouped scatter plot
#' appears in the first column of the histogram row. When sampling is requested,
#' the original grouped scatter plot appears in the first column of the first
#' sample row. If histogram and sample panels are both requested, the first
#' position of the histogram row is left blank to keep columns aligned.
#'
#' Sampling uses only the already computed histogram output and fixed
#' post-processing choices, so it does not consume an additional privacy
#' budget. A private histogram, especially a sparse histogram, can contain no
#' positive probability mass after thresholding. In that case no synthetic
#' points are generated for that group and histogram method, and a warning is
#' issued.
#'
#' @return A list with components:
#' \describe{
#'   \item{score}{The output of [dp_score_group()].}
#'   \item{sample}{Present only when `"sample"` is requested. A nested list
#'   containing combined grouped synthetic coordinates for `nonprivate` and
#'   each requested private histogram method. Each leaf is indexed by sampling
#'   method (`center` and/or `uniform`) and contains columns `pc_x`, `pc_y`,
#'   and `group`.}
#'   \item{plot}{A list containing the grouped original scatter plot, the
#'   overlaid non-private histogram, requested private histogram panels,
#'   grouped sampled-score panels under `plot$sample`, and the combined
#'   patchwork layout in `plot$all`.}
#'   \item{group_colors}{Named vector of colors used for the groups.}
#' }
#'
#' @seealso
#' [dp_score_group()] for computing group-wise score histograms without
#' plotting.
#' [dp_score_plot()] for pooled score histogram and sampled-score plots.
#'
#' @examples
#' data(gau_g, package = "dppca")
#'
#' set.seed(123)
#' score_plot_gau_g <- dp_score_plot_group(
#'   gau_g,
#'   group = "group",
#'   eps = 3,
#'   delta = 1e-3,
#'   bins = c(8, 8),
#'   sampling_control = sampling_control(
#'     method = c("center", "uniform")
#'   )
#' )
#'
#' score_plot_gau_g$plot$all
#'
#' @export
dp_score_plot_group <- function(
    X,
    group,
    eps,
    delta,
    bins,
    center = TRUE,
    standardize = FALSE,
    g_dppca = FALSE,
    cpp.option = FALSE,
    axes = c(1, 2),
    method = c("add", "sparse"),
    private_plot = c("histogram", "sample"),
    sampling_control = NULL,
    plot_control = NULL
) {
  method <- unique(
    match.arg(
      method,
      choices = c("add", "sparse"),
      several.ok = TRUE
    )
  )

  private_plot <- unique(
    match.arg(
      private_plot,
      choices = c("histogram", "sample"),
      several.ok = TRUE
    )
  )

  pctrl <- .merge_score_plot_control(plot_control)

  group_data <- split_group_input(X, group)
  group_vec <- group_data$group
  X_feat <- group_data$X
  g_levels <- as.character(unique(group_vec))

  if (!is.null(pctrl$group_colors)) {
    col_map <- complete_color_map(
      pctrl$group_colors,
      g_levels
    )
  } else if (.is_color_vec(g_levels)) {
    col_map <- stats::setNames(
      as.character(g_levels),
      g_levels
    )
  } else {
    pal <- grDevices::hcl.colors(
      length(g_levels),
      "Dark3"
    )
    col_map <- stats::setNames(pal, g_levels)
  }

  score_res <- dp_score_group(
    X = X_feat,
    group = group_vec,
    eps = eps,
    delta = delta,
    bins = bins,
    center = center,
    standardize = standardize,
    g_dppca = g_dppca,
    cpp.option = cpp.option,
    axes = axes,
    method = method
  )

  X_score <- as.data.frame(score_res$score)
  colnames(X_score) <- c("pc_x", "pc_y")
  X_score$group <- factor(
    as.character(group_vec),
    levels = g_levels
  )

  xlim <- score_res$frame$xlim
  ylim <- score_res$frame$ylim
  pc_names <- colnames(score_res$score)

  xlab <- if (is.null(pctrl$xlab)) pc_names[1] else pctrl$xlab
  ylab <- if (is.null(pctrl$ylab)) pc_names[2] else pctrl$ylab
  legend_position <- if (pctrl$show_group_legend) "right" else "none"

  p_scatter <- ggplot2::ggplot(
    X_score,
    ggplot2::aes(
      x = .data$pc_x,
      y = .data$pc_y,
      colour = .data$group
    )
  ) +
    ggplot2::geom_point(
      alpha = pctrl$scatter_alpha,
      size = pctrl$scatter_size
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
      base_size = pctrl$base_size,
      legend_position = legend_position
    ) +
    ggplot2::labs(
      x = xlab,
      y = ylab,
      colour = "Group"
    )

  p_scatter <- add_title_dp(
    p_scatter,
    pctrl$scatter_title,
    title_size = pctrl$title_size
  )

  bind_group_histograms <- function(which) {
    dplyr::bind_rows(
      lapply(g_levels, function(g) {
        df <- score_res$groups[[g]][[which]]
        df$group <- g
        df
      })
    )
  }

  hist_plot_out <- list()

  coord_nonprivate_all <- bind_group_histograms("nonprivate")
  hist_plot_out$nonprivate <- make_hist_all_dp(
    coord_nonprivate_all,
    xlim = xlim,
    ylim = ylim,
    col_map = col_map,
    title = pctrl$nonprivate_title,
    xlab = xlab,
    ylab = ylab,
    alpha_range = pctrl$hist_alpha_range,
    base_size = pctrl$base_size,
    title_size = pctrl$title_size,
    legend_position = legend_position
  )

  if ("add" %in% method) {
    coord_add_all <- bind_group_histograms("add")
    hist_plot_out$add <- make_hist_all_dp(
      coord_add_all,
      xlim = xlim,
      ylim = ylim,
      col_map = col_map,
      title = pctrl$add_title,
      xlab = xlab,
      ylab = ylab,
      alpha_range = pctrl$hist_alpha_range,
      base_size = pctrl$base_size,
      title_size = pctrl$title_size,
      legend_position = legend_position
    )
  }

  if ("sparse" %in% method) {
    coord_sparse_all <- bind_group_histograms("sparse")
    hist_plot_out$sparse <- make_hist_all_dp(
      coord_sparse_all,
      xlim = xlim,
      ylim = ylim,
      col_map = col_map,
      title = pctrl$sparse_title,
      xlab = xlab,
      ylab = ylab,
      alpha_range = pctrl$hist_alpha_range,
      base_size = pctrl$base_size,
      title_size = pctrl$title_size,
      legend_position = legend_position
    )
  }

  plot_out <- list(
    scatter = p_scatter,
    nonprivate = hist_plot_out$nonprivate
  )

  if ("add" %in% method) {
    plot_out$add <- hist_plot_out$add
  }
  if ("sparse" %in% method) {
    plot_out$sparse <- hist_plot_out$sparse
  }

  sample_out <- NULL
  sample_plot_out <- NULL
  sctrl <- NULL

  if ("sample" %in% private_plot) {
    sctrl <- .resolve_sampling_control(
      control = sampling_control,
      n = nrow(score_res$score)
    )

    group_sizes <- vapply(
      g_levels,
      function(g) score_res$groups[[g]]$n,
      numeric(1)
    )
    names(group_sizes) <- g_levels

    group_sample_sizes <- allocate_group_sample_sizes(
      group_sizes = group_sizes,
      total_size = sctrl$sample_size
    )

    targets <- c("nonprivate", method)
    sample_out <- list()
    sample_plot_out <- list()

    for (target in targets) {
      pieces_by_method <- stats::setNames(
        vector("list", length(sctrl$method)),
        sctrl$method
      )

      for (sm in sctrl$method) {
        pieces_by_method[[sm]] <- list()
      }

      for (g in g_levels) {
        n_g_sample <- group_sample_sizes[[g]]

        if (n_g_sample == 0L) {
          for (sm in sctrl$method) {
            pieces_by_method[[sm]][[g]] <- data.frame(
              pc_x = numeric(0),
              pc_y = numeric(0),
              group = character(0)
            )
          }
          next
        }

        sampled_g <- sample_private_score_histogram(
          hist_df = score_res$groups[[g]][[target]],
          sample_size = n_g_sample,
          sample_method = sctrl$method,
          bandwidth_scale = sctrl$bandwidth_scale
        )

        for (sm in sctrl$method) {
          df_g <- sampled_g[[sm]]

          # A sparse private histogram can legitimately contain no positive
          # probability mass after thresholding. In that case the sampling
          # helper returns a zero-row data frame. Use a zero-length group
          # vector so that empty sampled results remain valid data frames
          # instead of triggering a replacement-length error.
          df_g$group <- rep(g, nrow(df_g))

          pieces_by_method[[sm]][[g]] <- df_g
        }
      }

      sample_out[[target]] <- list()
      sample_plot_out[[target]] <- list()

      target_label <- switch(
        target,
        nonprivate = "Non-private",
        add = "Add",
        sparse = "Sparse"
      )

      for (sm in sctrl$method) {
        combined <- dplyr::bind_rows(
          pieces_by_method[[sm]]
        )

        combined$group <- factor(
          as.character(combined$group),
          levels = g_levels
        )

        sample_out[[target]][[sm]] <- combined

        sm_label <- if (sm == "center") "Center" else "Uniform"

        sample_plot_out[[target]][[sm]] <-
          make_group_sample_plot_dp(
            sample_df = combined,
            xlim = xlim,
            ylim = ylim,
            col_map = col_map,
            title = paste0(
              target_label,
              " Sample (",
              sm_label,
              ")"
            ),
            xlab = xlab,
            ylab = ylab,
            point_alpha = pctrl$scatter_alpha,
            point_size = pctrl$scatter_size,
            base_size = pctrl$base_size,
            title_size = pctrl$title_size,
            legend_position = legend_position
          )
      }
    }

    plot_out$sample <- sample_plot_out
  }

  # Use one flat row-major patchwork grid so columns align across rows.
  n_cols <- 2 + length(method)
  layout_panels <- list()

  if ("histogram" %in% private_plot) {
    if ("sample" %in% private_plot) {
      layout_panels[[length(layout_panels) + 1L]] <-
        patchwork::plot_spacer()
    } else {
      layout_panels[[length(layout_panels) + 1L]] <- p_scatter
    }

    layout_panels[[length(layout_panels) + 1L]] <-
      hist_plot_out$nonprivate

    for (m in method) {
      layout_panels[[length(layout_panels) + 1L]] <-
        hist_plot_out[[m]]
    }
  }

  if ("sample" %in% private_plot) {
    for (i in seq_along(sctrl$method)) {
      sm <- sctrl$method[i]

      if (i == 1L) {
        layout_panels[[length(layout_panels) + 1L]] <- p_scatter
      } else {
        layout_panels[[length(layout_panels) + 1L]] <-
          patchwork::plot_spacer()
      }

      layout_panels[[length(layout_panels) + 1L]] <-
        sample_plot_out$nonprivate[[sm]]

      for (m in method) {
        layout_panels[[length(layout_panels) + 1L]] <-
          sample_plot_out[[m]][[sm]]
      }
    }
  }

  plot_out$all <- patchwork::wrap_plots(
    layout_panels,
    ncol = n_cols,
    byrow = TRUE
  )

  out <- list(
    score = score_res
  )

  if (!is.null(sample_out)) {
    out$sample <- sample_out
  }

  out$plot <- plot_out
  out$group_colors <- col_map
  out
}
