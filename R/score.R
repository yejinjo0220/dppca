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
#' @param eps Positive number defining the total `epsilon` privacy parameter.
#' @param delta Number in `(0, 1)` defining the total `delta` privacy parameter.
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
#' @param fixed_frame Optional fixed plotting frame. If supplied, it can be a
#'   numeric vector `c(lower, upper)` used for both axes, or a list with numeric
#'   components `xlim` and `ylim`. If `NULL`, a private square frame is estimated
#'   from the score coordinates.
#'
#' @details
#' Let \eqn{v_a} and \eqn{v_b} be the principal component directions selected
#' by `axes = c(a, b)`. After preprocessing, the score point for observation
#' \eqn{i} is \eqn{s_i = (x_i^\top v_a, x_i^\top v_b)}. A non-private score
#' plot would display the points \eqn{s_1, \ldots, s_n} directly. This function
#' instead summarizes their empirical distribution by a two-dimensional histogram
#' and releases private versions of the histogram.
#'
#' If `fixed_frame = NULL`, the plotting frame is constructed privately. The two
#' score coordinates are stacked into one vector, private lower and upper
#' quantiles are estimated using a smooth-sensitivity based quantile mechanism
#' \insertCite{nissim2007smooth}{dppca}, and the resulting interval is used as a
#' common square frame for both axes. If `fixed_frame` is supplied, it is treated
#' as public and no privacy budget is spent on frame construction.
#'
#' The private histogram is computed on the rectangular grid defined by
#' `fixed_frame` or by the private frame and the bin counts in `bins`. Under
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
#' `g_dppca = FALSE` and `fixed_frame = NULL`, half of `eps` and `delta` is used
#' for private frame construction and half for the private histogram. If
#' `g_dppca = TRUE` and `fixed_frame = NULL`, the parameters are split equally
#' among private direction estimation, private frame construction, and private
#' histogram release. If `fixed_frame` is supplied, the frame step is skipped and
#' the remaining steps split the privacy parameters equally.
#'
#' @return A list with components:
#' \describe{
#'   \item{score}{An \eqn{n \times 2} matrix of score coordinates.}
#'   \item{frame}{A list with components `xlim` and `ylim`.}
#'   \item{none}{Data frame for the non-private empirical histogram.}
#'   \item{add}{Data frame for the additive Gaussian private histogram, or
#'   `NULL` if not requested.}
#'   \item{sparse}{Data frame for the sparse private histogram, or `NULL` if
#'   not requested.}
#'   \item{method}{Character vector of private histogram methods used.}
#' }
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
#' \insertRef{nissim2007smooth}{dppca}
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
#' X <- head(gau, 50)
#'
#' set.seed(123)
#' score_out <- dp_score(
#'   X,
#'   eps = 1,
#'   delta = 1e-2,
#'   bins = c(3, 3),
#'   method = "add"
#' )
#' head(score_out$none)
#' head(score_out$add)
#'
#' # Simulated low-rank data.
#' set.seed(123)
#' n <- 50
#' z1 <- rnorm(n)
#' z2 <- rnorm(n)
#' X_sim <- cbind(
#'   x1 = z1 + 0.2 * rnorm(n),
#'   x2 = 0.8 * z1 + 0.2 * rnorm(n),
#'   x3 = z2 + 0.2 * rnorm(n),
#'   x4 = 0.5 * z1 - 0.4 * z2 + 0.2 * rnorm(n)
#' )
#'
#' set.seed(123)
#' dp_score(X_sim, eps = 1, delta = 1e-2, bins = c(3, 3))
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
    axes = c(1, 2),
    fixed_frame = NULL
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
  m_x <- bins[1]
  m_y <- bins[2]
  method <- match.arg(method, choices = c("add", "sparse"), several.ok = TRUE)

  budget <- split_score_privacy_budget(
    eps = eps,
    delta = delta,
    g_dppca = g_dppca,
    fixed_frame = fixed_frame
  )

  score_res <- compute_score_coordinates(
    X = X,
    axes = axes,
    center = center,
    standardize = standardize,
    g_dppca = g_dppca,
    cpp.option = cpp.option,
    eps_dir = budget$eps_dir,
    delta_dir = budget$delta_dir
  )

  frame_out <- dp_frame(
    X = score_res$score,
    eps_frame = budget$eps_frame,
    delta_frame = budget$delta_frame,
    inflate = 0.20,
    q = NULL,
    fixed_frame = fixed_frame
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

  list(
    score = score_res$score,
    frame = frame_out,
    none = hist$none,
    add = hist$add,
    sparse = hist$sparse,
    method = method
  )
}

#' Plot differentially private score histograms
#'
#' This function computes and visualizes two-dimensional principal component
#' score histograms, including the original scatter plot, the non-private
#' empirical histogram, and one or more differentially private histogram
#' estimates. It is a plotting wrapper around [dp_score()] and returns both the
#' computed score output and `ggplot` objects.
#'
#' @inheritParams dp_score
#'
#' @return A list with components:
#' \describe{
#'   \item{score}{The output of [dp_score()].}
#'   \item{plot}{A list containing the scatter plot, histogram panels, and the
#'   combined patchwork plot.}
#' }
#'
#' @seealso
#' [dp_score()] for computing score histograms without plotting.
#' [dp_score_plot_group()] for group-wise score histogram plots.
#'
#' @examples
#' data(gau, package = "dppca")
#' X <- head(gau, 50)
#'
#' set.seed(123)
#' p <- dp_score_plot(
#'   X,
#'   eps = 1,
#'   delta = 1e-2,
#'   bins = c(3, 3),
#'   method = "add"
#' )
#' p$plot$all
#'
#' @export
dp_score_plot <- function(
    X,
    eps,
    delta,
    bins,
    method = c("add", "sparse"),
    center = TRUE,
    standardize = FALSE,
    g_dppca = FALSE,
    cpp.option = FALSE,
    axes = c(1, 2),
    fixed_frame = NULL
) {
  method <- match.arg(method, choices = c("add", "sparse"), several.ok = TRUE)
  color <- "#6A5ACD"

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
    method = method,
    fixed_frame = fixed_frame
  )

  X_score <- as.data.frame(score_res$score)
  colnames(X_score) <- c("pc_x", "pc_y")

  xlim <- score_res$frame$xlim
  ylim <- score_res$frame$ylim
  pc_names <- colnames(score_res$score)

  p_scatter <- ggplot2::ggplot(
    X_score,
    ggplot2::aes(x = .data$pc_x, y = .data$pc_y)
  ) +
    ggplot2::geom_point(alpha = 0.6, size = 1.8, color = color) +
    ggplot2::coord_fixed(xlim = xlim, ylim = ylim) +
    ggplot2::scale_x_continuous(expand = c(0, 0), breaks = pretty(xlim, n = 5)) +
    ggplot2::scale_y_continuous(expand = c(0, 0), breaks = pretty(ylim, n = 5)) +
    theme_dp_base() +
    ggplot2::labs(x = pc_names[1], y = pc_names[2])

  p_scatter <- add_title_dp(p_scatter, "Original Scatter")
  p_none <- make_hist_plot_dp(score_res$none, xlim, ylim, color, "Original Hist")
  p_add <- NULL
  p_sparse <- NULL

  plot_panels <- list(p_scatter, p_none)

  if ("add" %in% method) {
    p_add <- make_hist_plot_dp(score_res$add, xlim, ylim, color, "Add DP Hist")
    plot_panels <- c(plot_panels, list(p_add))
  }

  if ("sparse" %in% method) {
    p_sparse <- make_hist_plot_dp(score_res$sparse, xlim, ylim, color, "Sparse DP Hist")
    plot_panels <- c(plot_panels, list(p_sparse))
  }

  p_all <- patchwork::wrap_plots(plot_panels, nrow = 1)

  list(
    score = score_res,
    plot = list(
      scatter = p_scatter,
      none = p_none,
      add = p_add,
      sparse = p_sparse,
      all = p_all
    )
  )
}

#' Group-wise differentially private score histograms
#'
#' This function computes two-dimensional principal component scores and releases
#' group-wise differentially private histograms on a common score frame and grid.
#' It is useful when observations have group labels and the low-dimensional score
#' distribution should be compared across groups.
#'
#' @inheritParams dp_score
#' @param group Group labels. This can be a vector of length `nrow(X)` or a
#'   single column name in `X`. If a column name is supplied, that column is
#'   used as the group label and removed from the feature matrix.
#'
#' @details
#' The score directions, plotting frame, and histogram grid are shared across all
#' groups. For each group \eqn{g}, the group-specific count in bin \eqn{B_k} is
#' \eqn{c_k^{(g)} = \sum_i 1\{s_i \in B_k, g_i = g\}}. Private histograms are
#' then computed separately for each group on the common grid. Because the groups
#' form a partition of the rows, the group-wise histograms use the same histogram
#' privacy parameters for each group by parallel composition.
#'
#' @return A list with components:
#' \describe{
#'   \item{score}{An \eqn{n \times 2} matrix of score coordinates.}
#'   \item{frame}{A list with components `xlim` and `ylim`.}
#'   \item{groups}{A named list of group-specific histogram outputs.}
#'   \item{method}{Character vector of private histogram methods used.}
#' }
#'
#' @seealso
#' [dp_score_plot_group()] for plotting group-wise score histograms.
#' [dp_score()] for pooled score histograms.
#'
#' @examples
#' data(gau_g, package = "dppca")
#'
#' X <- head(gau_g, 60)
#'
#' set.seed(123)
#' group_out <- dp_score_group(
#'   X,
#'   group = "color",
#'   eps = 1,
#'   delta = 1e-2,
#'   bins = c(3, 3),
#'   method = "add"
#' )
#' names(group_out$groups)
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
    axes = c(1, 2),
    fixed_frame = NULL
) {
  method <- match.arg(method, choices = c("add", "sparse"), several.ok = TRUE)

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
    g_dppca = g_dppca,
    fixed_frame = fixed_frame
  )

  score_res <- compute_score_coordinates(
    X = X_mat,
    axes = axes,
    center = center,
    standardize = standardize,
    g_dppca = g_dppca,
    cpp.option = cpp.option,
    eps_dir = budget$eps_dir,
    delta_dir = budget$delta_dir
  )

  frame_out <- dp_frame(
    X = score_res$score,
    eps_frame = budget$eps_frame,
    delta_frame = budget$delta_frame,
    inflate = 0.20,
    q = NULL,
    fixed_frame = fixed_frame
  )

  grid <- score_histogram_grid(
    xlim = frame_out$xlim,
    ylim = frame_out$ylim,
    m_x = m_x,
    m_y = m_y
  )

  eps_hist_method <- budget$eps_hist / length(method)
  delta_hist_method <- budget$delta_hist / length(method)
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

    groups_out[[g]] <- list(
      n = n_g,
      none = hist$none,
      add = hist$add,
      sparse = hist$sparse
    )
  }

  list(
    score = score_res$score,
    frame = frame_out,
    groups = groups_out,
    method = method
  )
}

#' Plot group-wise differentially private score histograms
#'
#' This function computes and visualizes group-wise differentially private score
#' histograms. It is a plotting wrapper around [dp_score_group()] and returns
#' both the computed group-wise score output and `ggplot` objects.
#'
#' @inheritParams dp_score_group
#'
#' @return A list with components:
#' \describe{
#'   \item{score}{The output of [dp_score_group()].}
#'   \item{plot}{A list containing group-wise histogram plots.}
#'   \item{group_colors}{Named vector of colors used for the groups.}
#' }
#'
#' @seealso
#' [dp_score_group()] for computing group-wise score histograms without
#' plotting.
#' [dp_score_plot()] for pooled score histogram plots.
#'
#' @examples
#' data(gau_g, package = "dppca")
#'
#' X <- head(gau_g, 60)
#'
#' set.seed(123)
#' p <- dp_score_plot_group(
#'   X,
#'   group = "color",
#'   eps = 1,
#'   delta = 1e-5,
#'   bins = c(8, 8),
#'   method = "add"
#' )
#' p$plot$all
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
    fixed_frame = NULL
) {
  method <- match.arg(method, choices = c("add", "sparse"), several.ok = TRUE)

  group_data <- split_group_input(X, group)
  group_vec <- group_data$group
  X_feat <- group_data$X
  g_levels <- as.character(unique(group_vec))

  if (.is_color_vec(g_levels)) {
    col_map <- stats::setNames(as.character(g_levels), g_levels)
  } else {
    pal <- grDevices::hcl.colors(length(g_levels), "Dark3")
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
    method = method,
    fixed_frame = fixed_frame
  )

  X_score <- as.data.frame(score_res$score)
  colnames(X_score) <- c("pc_x", "pc_y")
  X_score$group <- as.character(group_vec)

  xlim <- score_res$frame$xlim
  ylim <- score_res$frame$ylim
  pc_names <- colnames(score_res$score)

  p_scatter <- ggplot2::ggplot(
    X_score,
    ggplot2::aes(x = .data$pc_x, y = .data$pc_y, colour = .data$group)
  ) +
    ggplot2::geom_point(alpha = 0.6, size = 1.8) +
    ggplot2::scale_colour_manual(values = col_map) +
    ggplot2::coord_fixed(xlim = xlim, ylim = ylim) +
    ggplot2::scale_x_continuous(expand = c(0, 0), breaks = pretty(xlim, n = 5)) +
    ggplot2::scale_y_continuous(expand = c(0, 0), breaks = pretty(ylim, n = 5)) +
    theme_dp_base() +
    ggplot2::labs(x = pc_names[1], y = pc_names[2])

  p_scatter <- add_title_dp(p_scatter, "Original Scatter")

  coord_none_all <- dplyr::bind_rows(lapply(names(score_res$groups), function(g) {
    df <- score_res$groups[[g]]$none
    df$group <- g
    df
  }))

  p_none_all <- make_hist_all_dp(coord_none_all, xlim, ylim, col_map, "Original Hist")
  p_add_all <- NULL
  p_sparse_all <- NULL

  plot_panels <- list(p_scatter, p_none_all)

  if ("add" %in% method) {
    coord_add_all <- dplyr::bind_rows(lapply(names(score_res$groups), function(g) {
      df <- score_res$groups[[g]]$add
      df$group <- g
      df
    }))
    p_add_all <- make_hist_all_dp(coord_add_all, xlim, ylim, col_map, "Add DP Hist")
    plot_panels <- c(plot_panels, list(p_add_all))
  }

  if ("sparse" %in% method) {
    coord_sparse_all <- dplyr::bind_rows(lapply(names(score_res$groups), function(g) {
      df <- score_res$groups[[g]]$sparse
      df$group <- g
      df
    }))
    p_sparse_all <- make_hist_all_dp(coord_sparse_all, xlim, ylim, col_map, "Sparse DP Hist")
    plot_panels <- c(plot_panels, list(p_sparse_all))
  }

  p_all <- patchwork::wrap_plots(plot_panels, nrow = 1)

  make_group_layout <- function(which = c("none", "add", "sparse"), title_prefix) {
    which <- match.arg(which)

    group_plots <- lapply(g_levels, function(g) {
      df_g <- score_res$groups[[g]][[which]]
      title_g <- if (.is_color_vec(g_levels)) NULL else paste0(title_prefix, g)

      make_hist_single_dp(
        df_g,
        xlim = xlim,
        ylim = ylim,
        col = col_map[[g]],
        title = title_g
      )
    })
    names(group_plots) <- g_levels

    list(
      all = patchwork::wrap_plots(group_plots, ncol = 3),
      groups = group_plots
    )
  }

  none_layout <- make_group_layout("none", "Original Hist: ")
  add_layout <- if ("add" %in% method) make_group_layout("add", "Add DP Hist: ") else NULL
  sparse_layout <- if ("sparse" %in% method) make_group_layout("sparse", "Sparse DP Hist: ") else NULL

  list(
    score = score_res,
    plot = list(
      scatter = p_scatter,
      none = none_layout,
      add = add_layout,
      sparse = sparse_layout,
      all = p_all
    ),
    group_colors = col_map
  )
}
