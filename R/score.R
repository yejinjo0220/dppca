# ============================================================
# score.R
# Functions for differentially private score histograms
# ============================================================

#' Differentially private score histograms
#'
#' This function projects observations onto supplied private directions or
#' directions estimated by [dp_loading()] and computes two-dimensional score
#' histograms. Observed coordinates and empirical reference histograms are
#' returned only when `non_private = TRUE`.
#'
#' @param X A numeric matrix or data frame with observations in rows and variables
#'   in columns. Non-numeric data-frame columns are excluded with a warning;
#'   numeric values must be finite. At least two observations and two numeric
#'   features are required. If directions must be estimated, the
#'   pairwise-separation condition of [dp_loading()] also applies.
#' @param privacy A [privacy_control()] object allocating budget to `loading`
#'   and `score` (equal shares by default). Each requested histogram method
#'   receives the full histogram allocation.
#' @param bins Positive integer vector of length two giving the number of
#'   histogram bins along the selected score axes.
#' @param method One or both histogram methods: `"add"` (Gaussian additive noise)
#'   and `"sparse"` (thresholded sparse histogram). Default:
#'   `c("add", "sparse")`.
#' @param cpp.option Whether to use the compiled spherical-Kendall implementation
#'   (`TRUE` by default) when estimating private directions; ignored if
#'   `V_dp` is supplied.
#' @param axes Two distinct PC indices for the horizontal and vertical score
#'   axes, respectively (`c(1, 2)` by default). Each must be between `1`
#'   and the number of numeric PCA features.
#' @param V_dp Optional `p` by `p` private PC direction matrix (`NULL` by default),
#'   where `p` is the number of numeric PCA features. Columns must be
#'   orthonormal and rows must follow the feature order in `X`; matching
#'   row names are required when both sets of names exist.
#' @param non_private Whether to compute and store ordinary PCA reference results
#'   (`FALSE` by default). For score analysis, these include mean-centered ordinary
#'   PCA scores and empirical reference histograms.
#' @param group Optional group column name in `X` or vector of group labels
#'   (`NULL` by default). A selected group column is excluded from PCA
#'   features, even when numeric; group labels and membership are assumed
#'   public. All groups share the same private directions,
#'   frame, and histogram grid.
#'
#' @details
#' Scores are projected onto the selected private PC directions. Rather than
#' releasing individual projected coordinates, the function estimates a private
#' two-dimensional plotting frame followed by private histogram probabilities.
#' The frame uses two private coordinate medians and an inflated private
#' 0.99 quantile of distances from their center (inflation factor `1.2`),
#' using unbounded quantile mechanisms.
#'
#' For a fixed frame and directions under replacement adjacency, the histogram
#' count vector has `l1` sensitivity at most `2` and `l2` sensitivity at most
#' `sqrt(2)`. Out-of-frame scores are assigned to boundary bins. Adaptive
#' private directions and frame selection require privacy composition. Input is
#' not automatically centered or standardized; data-dependent mean centering
#' or standardization is not covered by the count sensitivity argument.
#'
#' Two histogram mechanisms are available:
#' \itemize{
#'   \item `"add"`: Add Gaussian noise to all counts, truncate negatives to
#'   zero, and normalize. The classical Gaussian calibration has a standard
#'   sufficient guarantee for `0 < eps_hist < 1`; larger values require a
#'   separate calibration check. All-zero truncated counts cause an error.
#'   \item `"sparse"`: Add Laplace noise to nonzero bin proportions, retain
#'   bins above a privacy-dependent threshold, and normalize positive mass.
#'   If none survive, return an all-zero histogram.
#' }
#'
#' The `score` epsilon allocation is split 35% for the private frame and 65%
#' for the histogram; all score delta goes to the histogram. The frame epsilon
#' is divided equally among two private medians and one radial quantile.
#' Directions and frame are shared across methods. Each requested histogram
#' method receives the full histogram allocation, so jointly releasing `"add"`
#' and `"sparse"` requires composition across methods.
#'
#' For grouped data, the partition and labels are assumed public and fixed;
#' groups share directions, frame and grid. Parallel composition applies to
#' disjoint group histograms under this assumption, while different methods
#' still require composition. `non_private = FALSE` does not privatize group
#' membership. Input dimensions, feature schema, bins, and tuning choices are
#' assumed public. Repeated releases also require composition.
#'
#' The optional ordinary PCA reference uses sample-mean-centered scores and a
#' zero-centered comparison frame. It is not differentially private. See
#' [dp_loading()] for direction-estimation and Gaussian calibration conditions.
#'
#' @return A `dp_score` object containing the private `frame` and requested
#'   histogram methods (`add` and/or `sparse`), or `groups` for grouped fits.
#'   Each histogram has `xmin`, `xmax`, `ymin`, `ymax` and `prob` columns.
#'   `nonprivate` is `NULL` by default; when requested it contains `scores`,
#'   `frame` and the empirical `histogram` (or group histograms under `groups`).
#'   `metadata` contains data dimensions, axes, bins, methods and privacy settings.
#'   Use `summary()` and `plot()` to inspect or visualize stored results.
#'
#' @seealso
#' [plot.dp_score()] for plotting stored results without new private estimates.
#' [dp_pca()] and [plot.dp_pca()] for estimation followed by plots of a stored fit.
#' [plot.dp_score_group()] for displaying a grouped result.
#' [dp_loading()] for private principal component direction estimation.
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
#' # Release a private histogram with fixed preprocessing and no raw scores.
#' set.seed(123)
#' score_gau <- dp_score(
#'   X,
#'   privacy = privacy_control(eps = 1.8, delta = 1e-3),
#'   method = "add",
#'   bins = c(10, 10),
#'   non_private = FALSE
#' )
#'
#' score_gau$frame
#' head(score_gau$add)
#' summary(score_gau)
#' # plot(score_gau, type = "histogram")
#' # To use a group column: dp_score(X, privacy, bins, group = "region")
#'
#' @importFrom VGAM rlaplace
#' @export
dp_score <- function(
    X,
    privacy,
    bins,
    method = c("add", "sparse"),
    cpp.option = TRUE,
    axes = c(1, 2),
    V_dp = NULL,
    non_private = FALSE,
    group = NULL
) {
  validate_logical_value(non_private, "non_private")
  privacy <- .validate_privacy_control(privacy, c("loading", "score"))
  input <- .prepare_pca_input(X, group = group, min_cols = 2L)
  X <- input$X
  validate_score_common(X, bins, cpp.option, axes)

  axes <- as.integer(axes)
  bins <- as.integer(bins)
  method <- unique(match.arg(method, c("add", "sparse"), several.ok = TRUE))
  budget <- split_score_privacy_budget(privacy)

  # Estimate directions and the frame once, whether or not groups are used.
  score_res <- compute_score_coordinates(
    X = X, axes = axes, cpp.option = cpp.option,
    eps_pc = budget$eps_pc, delta_pc = budget$delta_pc, V_dp = V_dp
  )
  frame_out <- dp_frame(
    X = score_res$score, eps_frame = budget$eps_frame, inflate = 0.20
  )

  # Private histogram estimates always use the private projected scores.
  if (is.null(input$group)) {
    hist <- score_histograms(
      X_score = score_res$score,
      xlim = frame_out$xlim, ylim = frame_out$ylim, bins = bins,
      eps_hist = budget$eps_hist, delta_hist = budget$delta_hist,
      method = method, non_private = FALSE
    )
    out <- c(list(frame = frame_out), hist)
  } else {
    groups <- .score_group_histograms(
      scores = score_res$score, group = input$group,
      frame = frame_out, bins = bins,
      eps_hist = budget$eps_hist, delta_hist = budget$delta_hist,
      method = method, non_private = FALSE
    )
    out <- list(frame = frame_out, groups = groups)
  }

  nonprivate <- NULL
  if (non_private) {
    # Ordinary PCA reference, sign-aligned for comparison with private axes.
    V_np <- .nonprivate_pca(X)$directions[, axes, drop = FALSE]
    signs <- ifelse(colSums(V_np * score_res$directions) < 0, -1, 1)
    V_np <- sweep(V_np, 2L, signs, "*")
    X_centered <- sweep(X, 2L, colMeans(X), "-")
    ordinary_scores <- X_centered %*% V_np
    colnames(ordinary_scores) <- paste0("PC", axes)

    # Same bin widths as private, but centered on zero for ordinary PCA.
    ordinary_frame <- .dppca_centered_frame(frame_out)
    nonprivate <- list(scores = ordinary_scores, frame = ordinary_frame)
    if (is.null(input$group)) {
      nonprivate$histogram <- .dppca_reference_histogram(
        ordinary_scores, ordinary_frame, bins
      )
    } else {
      nonprivate$groups <- lapply(names(out$groups), function(g) {
        rows <- input$group == g
        list(
          n = sum(rows),
          histogram = .dppca_reference_histogram(
            ordinary_scores[rows, , drop = FALSE], ordinary_frame, bins
          )
        )
      })
      names(nonprivate$groups) <- names(out$groups)
    }
  }

  out["nonprivate"] <- list(nonprivate)
  metadata <- list(
    n = nrow(X), p = ncol(X), axes = axes, bins = bins, methods = method,
    privacy = privacy, group_name = input$group_name,
    group_levels = input$group_levels,
    excluded_columns = input$excluded_columns
  )
  if (non_private && !is.null(input$group)) {
    metadata$group <- input$group
  }
  out$metadata <- metadata
  class(out) <- if (is.null(input$group)) {
    c("dp_score", "list")
  } else {
    c("dp_score_group", "dp_score", "list")
  }
  out
}

# Internal compatibility view for the existing score plotting/summary helpers.
# Public dp_score() fits keep only the new, organized return structure.
.score_legacy_view <- function(x) {
  out <- x
  out$metadata <- NULL
  if (is.null(x$nonprivate)) {
    out$nonprivate <- NULL
  } else {
    out$score <- x$nonprivate$scores
    out$nonprivate_frame <- x$nonprivate$frame
    if (is.null(x$groups)) {
      out$nonprivate <- x$nonprivate$histogram
    } else {
      for (g in names(out$groups)) {
        out$groups[[g]]$n <- x$nonprivate$groups[[g]]$n
        out$groups[[g]]$nonprivate <- x$nonprivate$groups[[g]]$histogram
      }
      out$nonprivate <- NULL
    }
  }
  legacy_metadata <- x$metadata
  legacy_metadata$method <- legacy_metadata$methods
  attr(out, "metadata") <- legacy_metadata
  out
}

#' Plot a stored private score histogram
#'
#' Reuses the stored frame and histograms. Synthetic sampling is post-processing;
#' no direction, frame or histogram mechanism is run again.
#'
#' @param x A fitted [dp_score()] object.
#' @param type `"histogram"`, `"sample"`, or both (`c("histogram", "sample")`
#'   by default). Samples are synthetic for DP histograms; ordinary-reference
#'   sample plots show centered observed scores unless an explicit
#'   `sample_size` is requested.
#' @param sampling_control Optional [sampling_control()] object used for
#'   `"sample"`; `NULL` uses defaults.
#' @param plot_control Optional [score_plot_control()] appearance settings;
#'   `NULL` uses defaults.
#' @param center Whether to shift displayed private coordinates by the midpoint
#'   of the private frame (`TRUE` by default). Ordinary reference coordinates
#'   remain centered at zero. This changes only the display.
#' @param ... Unused arguments are rejected.
#' @details
#' Non-private reference panels appear automatically when retained at
#' fitting. Combined grouped displays show one legend per histogram method;
#' individual plots retain their legends. Groups with no sampled points are
#' omitted from sampled-point legends.
#'
#' @return Invisibly returns a `dp_score_plot` object with `score`, optional
#'   `sample`, and `plot` components, while displaying `plot$all`.
#' @examples
#' # After estimating once:
#' # summary(score_fit)
#' # plots <- plot(score_fit, type = "histogram")
#' # plot(score_fit, type = "sample",
#' #      sampling_control = sampling_control(sample_size = 500))
#' @export
plot.dp_score <- function(
    x,
    type = c("histogram", "sample"),
    sampling_control = NULL,
    plot_control = NULL,
    center = TRUE,
    ...
) {
  validate_logical_value(center, "center")
  if (length(list(...))) {
    stop(
      "Unused arguments in `...`; use `type`, `sampling_control`, ",
      "`plot_control` or `center`.",
      call. = FALSE
    )
  }
  non_private <- !is.null(x$nonprivate)
  out <- .render_dp_score(
    .score_legacy_view(x), type, sampling_control, plot_control, non_private,
    center = center
  )
  out$score <- x
  print(out$plot$all)
  invisible(out)
}

#' Plot stored group-wise private score histograms
#'
#' Uses the same display options as [plot.dp_score()]. Group sample sizes
#' follow observed proportions if non-private counts were retained; otherwise
#' they are allocated equally across public group labels.
#'
#' @param x A grouped result from [dp_score()] with `group` supplied.
#' @inheritParams plot.dp_score
#' @return Invisibly returns a `dp_score_plot_group` object while displaying
#'   `plot$all`.
#' @export
plot.dp_score_group <- function(
    x,
    type = c("histogram", "sample"),
    sampling_control = NULL,
    plot_control = NULL,
    center = TRUE,
    ...
) {
  validate_logical_value(center, "center")
  if (length(list(...))) {
    stop(
      "Unused arguments in `...`; use `type`, `sampling_control`, ",
      "`plot_control` or `center`.",
      call. = FALSE
    )
  }
  non_private <- !is.null(x$nonprivate)
  out <- .render_dp_score(
    .score_legacy_view(x), type, sampling_control, plot_control, non_private,
    center = center
  )
  out$score <- x
  print(out$plot$all)
  invisible(out)
}

#' Summarize stored score histograms
#'
#' @param object A result from [dp_score()], with or without `group`.
#' @param ... Reserved for future arguments.
#' @return A summary containing the public sample size, axes, grid, private
#'   frame (including its center), allocated privacy budget and a table of
#'   positive-bin counts and probability sums. No data-dependent private
#'   mechanisms are run.
#' @rdname dp_score
#' @export
summary.dp_score <- function(object, ...) {
  .score_display_summary(object)
}

#' @rdname dp_score
#' @export
summary.dp_score_group <- function(object, ...) {
  .score_display_summary(object, grouped = TRUE)
}

# Add display metadata without changing the histogram-summary helper.
.score_display_summary <- function(object, grouped = FALSE) {
  out <- .score_summary(.score_legacy_view(object), grouped = grouped)
  metadata <- object$metadata
  out$p <- metadata$p
  out$methods <- metadata$methods
  budget <- split_score_privacy_budget(metadata$privacy)
  out$budget <- rbind(
    Loading = c(eps = budget$eps_pc, delta = budget$delta_pc),
    `Score frame` = c(eps = budget$eps_frame, delta = 0),
    `Score histogram` = c(eps = budget$eps_hist, delta = budget$delta_hist),
    Total = metadata$privacy$total
  )
  out
}

#' Print a score-histogram summary
#'
#' @param x A score-histogram summary.
#' @param digits Number of significant digits shown in tables.
#' @param ... Passed to table printing.
#' @return Invisibly returns `x`.
#' @rdname dp_score
#' @export
print.summary.dp_score <- function(x, digits = 4, ...) {
  cat("DP Score Histograms Summary\n")
  cat("Data:", x$n, "observations x", x$p, "variables\n")
  cat("Score axes:", paste0("PC", x$axes, collapse = ", "), "\n")
  cat("Grid:", paste(x$bins, collapse = " x "), "\n")
  cat("Methods: score =", paste(x$methods, collapse = ", "), "\n")
  if (!is.null(x$group_name)) {
    cat("Group variable:", x$group_name, "\n")
    cat("Groups:", paste(x$group_levels, collapse = ", "), "\n")
  }
  if (length(x$excluded_columns)) {
    cat("Excluded columns:", paste(x$excluded_columns, collapse = ", "), "\n")
  }

  frame <- data.frame(
    Axis = paste0("PC", x$axes),
    Center = c(mean(x$frame$xlim), mean(x$frame$ylim)),
    Lower = c(x$frame$xlim[1L], x$frame$ylim[1L]),
    Upper = c(x$frame$xlim[2L], x$frame$ylim[2L])
  )
  cat("\nPrivate frame:\n")
  print(frame, digits = digits, row.names = FALSE, ...)

  histograms <- x$histograms
  names(histograms) <- c(
    group = "Group", method = "Method", bins = "Bins",
    positive_bins = "Positive bins", probability_sum = "Probability sum"
  )[names(histograms)]
  cat("\nHistogram summary:\n")
  print(histograms, digits = digits, row.names = FALSE, ...)

  cat("\nAllocated privacy budget",
      if (length(x$methods) > 1L) " (per histogram method)" else "",
      ":\n", sep = "")
  print(x$budget, digits = digits, ...)
  if (length(x$methods) > 1L) {
    cat("\nAll histogram methods use the same private PC directions and frame.\n")
  }
  cat("\nNon-private reference:",
      if (x$non_private) "included" else "not stored", "\n")
  if (!is.null(x$samples)) {
    cat("\nStored synthetic points:\n")
    print(x$samples, digits = digits, row.names = FALSE, ...)
  }
  invisible(x)
}

#' Print stored score histogram results
#'
#' @param x A result from [dp_score()], with or without `group`.
#' @param ... Reserved for S3 generic compatibility.
#' @return Invisibly returns `x`.
#' @rdname dp_score
#' @export
print.dp_score <- function(x, ...) {
  s <- summary(x)
  cat("Differentially Private Score Histograms\n")
  cat("Data:", s$n, "observations x", s$p, "variables\n")
  cat("Score axes:", paste0("PC", s$axes, collapse = ", "), "\n")
  cat("Grid:", paste(s$bins, collapse = " x "), "\n")
  cat("Methods: score =", paste(s$methods, collapse = ", "), "\n")
  if (!is.null(s$group_name)) {
    cat("Group variable:", s$group_name, "\n")
    cat("Groups:", paste(s$group_levels, collapse = ", "), "\n")
  }
  budget <- s$privacy$total
  cat("Total privacy budget",
      if (length(s$methods) > 1L) " (per histogram method)" else "",
      ": eps = ", format(budget[["eps"]], digits = 6),
      "; delta = ", format(budget[["delta"]], digits = 6), "\n", sep = "")
  cat("Non-private reference:",
      if (s$non_private) "included" else "not stored", "\n")
  invisible(x)
}

#' @rdname dp_score
#' @export
print.dp_score_group <- function(x, ...) {
  print.dp_score(x, ...)
}

#' Redraw a stored score comparison plot
#'
#' @param x A plotting result returned by [plot.dp_score()] or
#'   [plot.dp_score_group()].
#' @param ... Optional arguments for [plot.dp_score()] or
#'   [plot.dp_score_group()]. With no arguments, the stored layout is redrawn.
#'   Explicit controls rebuild the display from stored estimates; sampled
#'   displays may generate new synthetic points, with no new private estimates.
#' @return Invisibly returns the original or rebuilt comparison result.
#' @rdname plot.dp_score
#' @export
plot.dp_score_plot <- function(x, ...) {
  if (length(list(...))) {
    return(plot(x$score, ...))
  }
  print(x$plot$all)
  invisible(x)
}

#' @rdname plot.dp_score_group
#' @export
plot.dp_score_plot_group <- function(x, ...) {
  if (length(list(...))) {
    return(plot(x$score, ...))
  }
  print(x$plot$all)
  invisible(x)
}

#' Summarize stored score comparison plots
#'
#' @param object A plotting result returned by [plot.dp_score()] or
#'   [plot.dp_score_group()].
#' @param ... Reserved for future arguments.
#' @return A score-histogram summary with stored synthetic point counts.
#' @rdname plot.dp_score
#' @export
summary.dp_score_plot <- function(object, ...) {
  object$score <- .score_legacy_view(object$score)
  .score_plot_summary(object)
}

#' @rdname plot.dp_score_group
#' @export
summary.dp_score_plot_group <- function(object, ...) {
  object$score <- .score_legacy_view(object$score)
  .score_plot_summary(object)
}

#' Print stored score comparison results
#'
#' @param x A plotting result returned by [plot.dp_score()] or
#'   [plot.dp_score_group()].
#' @param ... Passed to the summary print method.
#' @return Invisibly returns `x`.
#' @rdname plot.dp_score
#' @export
print.dp_score_plot <- function(x, ...) {
  print(summary(x), ...)
  invisible(x)
}

#' @rdname plot.dp_score_group
#' @export
print.dp_score_plot_group <- function(x, ...) {
  print(summary(x), ...)
  invisible(x)
}
