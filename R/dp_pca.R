# dp_pca.R

#' Estimate shared PCA quantities once
#'
#' @param X A numeric matrix or data frame with observations in rows and variables
#'   in columns. Non-numeric data-frame columns are excluded with a warning;
#'   numeric values must be finite. At least two observations and two numeric
#'   features are required. If directions must be estimated, the
#'   pairwise-separation condition of [dp_loading()] also applies.
#' @param privacy A [privacy_control()] object allocating budget to `loading`,
#'   `scree`, and `score` (equal thirds by default).
#' @param scree Required named list with `k = NULL` (all components),
#'   `method = "clipped"`, `control` (required), and `mono = TRUE`.
#'   `k`, when supplied, must be between `2` and the number of numeric features.
#'   Use [clipped_control()], [pmwm_control()], or [huber_control()] for tuning.
#' @param score Named list specifying `method = "sparse"` or `"add"` and
#'   `bins = c(20, 20)` by default. Exactly one histogram method is used.
#' @param axes Two distinct PC indices between `1` and the fitted `k`
#'   (`c(1, 2)` by default).
#' @param cpp.option Whether to use the compiled spherical-Kendall implementation
#'   (`TRUE` by default) when estimating private directions; ignored if
#'   `V_dp` is supplied.
#' @param non_private Whether to compute and store ordinary PCA reference results
#'   (`FALSE` by default).
#' @param V_dp Optional `p` by `p` private PC direction matrix (`NULL` by default),
#'   where `p` is the number of numeric PCA features. Columns must be
#'   orthonormal and rows must follow the feature order in `X`; matching
#'   row names are required when both sets of names exist.
#' @param group Optional group column name in `X` or vector of group labels
#'   (`NULL` by default). A selected group column is excluded from PCA
#'   features, even when numeric; group labels and membership are assumed
#'   public.
#'
#' @return A `dp_pca` object containing private `directions`, `scree`,
#'   `pve`, `hist`, `frame`, and optional `score_groups`.
#'   `nonprivate` contains the optional ordinary PCA reference, or `NULL`
#'   when `non_private = FALSE`. `metadata` contains data dimensions,
#'   privacy allocations, methods, controls, grouping information and axes.
#'   The print method returns the fitted object invisibly.
#'
#' @details
#' A single private PC direction matrix is estimated or supplied and reused
#' for scree and score estimation. Scree values estimate variance along those
#' directions using disjoint paired differences; an odd final observation is
#' unused for scree, while the score histogram uses all rows. PVE is normalized
#' over the retained `k` components, not necessarily all original variables.
#'
#' The selected scree method uses the entire scree allocation. Of the score
#' epsilon budget, 35% goes to the private frame and 65% to the histogram; all
#' score delta goes to the histogram. With supplied `V_dp`, the loading budget
#' records its previous release rather than being reassigned.
#'
#' Grouped data use one direction estimate and private frame. Their disjoint
#' histograms use parallel composition when group labels and membership are
#' public and fixed. The pooled histogram is post-processed from the grouped
#' histograms.
#'
#' Input is not automatically centered or standardized. Data-dependent
#' preprocessing needs separate privacy accounting, and sample size, schema,
#' selected components, and tuning settings are assumed public. The loading
#' conditions in [dp_loading()] also apply. The optional non-private reference
#' uses sample-mean-centered ordinary PCA scores and is not private.
#'
#' @examples
#' data(gau, package = "dppca")
#' set.seed(123)
#' fit <- dp_pca(
#'   gau[1:200, ],
#'   privacy = privacy_control(eps = 5, delta = 1e-4),
#'   scree = list(k = 2, control = clipped_control(C_clip = 3)),
#'   score = list(method = "sparse", bins = c(10, 10)),
#'   non_private = FALSE
#' )
#' fit
#' summary(fit)
#' fit$directions
#' fit$nonprivate # NULL by default
#' fit$metadata$privacy
#'
#' set.seed(123)
#' plots <- plot(fit) # DP panels only
#' plots$private$plots$biplot
#'
#' @export
dp_pca <- function(
    X,
    privacy,
    scree,
    score = list(),
    axes = c(1, 2),
    cpp.option = TRUE,
    non_private = FALSE,
    V_dp = NULL,
    group = NULL
) {
  if (missing(scree)) {
    stop(
      "Supply scree with a method-specific control, for example ",
      "scree = list(control = clipped_control(C_clip = ...))."
    )
  }

  scree_options <- .dppca_options(
    scree,
    list(k = NULL, method = "clipped", control = NULL, mono = TRUE),
    "scree"
  )
  score_options <- .dppca_options(
    score,
    list(method = "sparse", bins = c(20, 20)),
    "score"
  )

  k <- scree_options$k
  mono <- scree_options$mono

  for (options in list(scree_options, score_options)) {
    if (
      !is.character(options$method) || length(options$method) != 1L ||
        is.na(options$method)
    ) {
      stop("Choose exactly one scree method and one score method.")
    }
  }

  scree_method <- match.arg(scree_options$method, c("clipped", "pmwm", "huber"))
  score_method <- match.arg(score_options$method, c("sparse", "add"))
  bins <- score_options$bins

  for (argument in c("cpp.option", "mono", "non_private")) {
    validate_logical_value(get(argument), argument)
  }

  privacy <- .validate_privacy_control(
    privacy, c("loading", "scree", "score")
  )
  budgets <- privacy$budgets

  control <- .merge_scree_control(
    scree_method,
    scree_options$control
  )
  control_constructor <- switch(
    scree_method,
    clipped = clipped_control,
    pmwm = pmwm_control,
    huber = huber_control
  )
  control <- do.call(control_constructor, control)

  if (scree_method == "pmwm" && control$eta >= 0.5) {
    stop("eta must be less than 0.5.")
  }

  prepared <- .prepare_pca_input(X, group = group, min_cols = 2L)
  X_proc <- prepared$X
  n <- nrow(X_proc)
  p <- ncol(X_proc)

  if (p < 2) {
    stop("At least two variables are needed for the five PCA plots.")
  }
  if (is.null(k)) {
    k <- p
  }
  .dppca_number(k, "k", 2, p, integer = TRUE)

  if (
    !is.numeric(axes) || length(axes) != 2L || any(!is.finite(axes)) ||
    any(axes != floor(axes)) || any(axes < 1 | axes > k) || anyDuplicated(axes)
  ) {
    stop("axes must be two distinct integers between 1 and k.")
  }

  if (length(bins) == 1L) {
    bins <- rep(bins, 2)
  }
  validate_bins(bins)

  if (is.null(colnames(X_proc))) {
    colnames(X_proc) <- paste0("X", seq_len(p))
  }
  if (anyDuplicated(colnames(X_proc))) {
    stop("Variable names must be unique.")
  }
  estimator <- switch(
    scree_method,
    clipped = dp_scree_clipped,
    pmwm = dp_scree_pmwm,
    huber = dp_scree_huber
  )

  # Estimate or reuse one full loading matrix and share its projected scores.
  V <- .get_dp_directions(
    X_proc,
    budgets$loading,
    cpp.option,
    V_dp
  )
  Y <- X_proc %*% V

  scree_args <- c(
    list(
      Y = Y[, seq_len(k), drop = FALSE],
      eps = budgets$scree[["eps"]],
      delta = budgets$scree[["delta"]],
      mono = FALSE
    ),
    control
  )
  scree_result <- do.call(estimator, scree_args)
  eigenvalues <- scree_result$scree

  if (any(!is.finite(eigenvalues))) {
    stop("Non-finite eigenvalues from scree estimation.")
  }
  if (mono) {
    eigenvalues <- scree_post_processing(eigenvalues)
  }
  names(eigenvalues) <- colnames(V)[seq_len(k)]

  # Loading has its own allocation; the score budget covers frame and histogram.
  score_budget <- split_score_privacy_budget(privacy)
  frame <- dp_frame(
    Y[, axes, drop = FALSE],
    score_budget$eps_frame,
    inflate = 0.2
  )
  score_groups <- NULL
  if (is.null(prepared$group)) {
    histograms <- score_histograms(
      Y[, axes, drop = FALSE],
      frame$xlim,
      frame$ylim,
      bins,
      score_budget$eps_hist,
      score_budget$delta_hist,
      score_method,
      non_private = FALSE
    )
  } else {
    score_groups <- .score_group_histograms(
      Y[, axes, drop = FALSE], prepared$group, frame, bins,
      score_budget$eps_hist, score_budget$delta_hist, score_method,
      non_private = FALSE
    )
    histograms <- .score_pool_histograms(score_groups, score_method)
  }

  nonprivate <- NULL
  if (non_private) {
    ordinary <- .nonprivate_pca(X_proc)
    ordinary_eigenvalues <- ordinary$eigenvalues
    ordinary_directions <- ordinary$directions

    ordinary_eigenvalues <- stats::setNames(
      pmax(ordinary_eigenvalues[seq_len(k)], 0),
      colnames(V)[seq_len(k)]
    )

    # Align the ordinary reference to the private directions for comparison.
    ordinary_directions <- sweep(
      ordinary_directions,
      2,
      ifelse(colSums(ordinary_directions * V) < 0, -1, 1),
      "*"
    )
    # Ordinary PCA scores are always sample-mean centered.
    X_centered <- sweep(X_proc, 2L, colMeans(X_proc), "-")
    ordinary_scores <- X_centered %*% ordinary_directions[, axes, drop = FALSE]

    # Match the private frame's width, but center the reference grid at (0, 0).
    ordinary_frame <- .dppca_centered_frame(frame)
    ordinary_histogram <- .dppca_reference_histogram(
      ordinary_scores, ordinary_frame, bins
    )

    nonprivate <- list(
      directions = ordinary_directions,
      scree = ordinary_eigenvalues,
      pve = scree_to_pve(ordinary_eigenvalues),
      scores = ordinary_scores,
      hist = ordinary_histogram,
      frame = ordinary_frame,
      variable_sd = apply(X_proc, 2, stats::sd)
    )
    if (!is.null(prepared$group)) {
      nonprivate$group <- prepared$group
      nonprivate$score_groups <- .dppca_reference_groups(
        ordinary_scores, prepared$group, ordinary_frame, bins
      )
    }
  }

  structure(
    list(
      directions = V,
      scree = eigenvalues,
      pve = scree_to_pve(eigenvalues),
      hist = histograms[[score_method]],
      score_groups = score_groups,
      frame = frame,
      metadata = list(
        n = n,
        p = p,
        k = k,
        axes = as.integer(axes),
        group_name = prepared$group_name,
        group_levels = prepared$group_levels,
        excluded_columns = prepared$excluded_columns,
        privacy = privacy,
        score_budget = score_budget,
        methods = list(scree = scree_method, score = score_method),
        controls = list(
          scree = list(k = k, method = scree_method, control = control, mono = mono),
          score = list(method = score_method, bins = bins)
        )
      ),
      nonprivate = nonprivate
    ),
    class = "dp_pca"
  )
}

# Empirical reference grid has the same span as the DP frame, centered at zero.
.dppca_centered_frame <- function(frame) {
  list(
    xlim = c(-1, 1) * diff(frame$xlim) / 2,
    ylim = c(-1, 1) * diff(frame$ylim) / 2
  )
}

# Empirical reference on its own centered grid: right-closed bins, with scores
# outside the frame assigned to the corresponding boundary bins.
.dppca_reference_histogram <- function(
    scores,
    frame,
    bins
) {
  grid <- score_histogram_grid(
    frame$xlim,
    frame$ylim,
    bins[1],
    bins[2]
  )

  bin_index <- function(
      values,
      breaks
  ) {
    index <- cut(values, breaks = breaks, include.lowest = TRUE, labels = FALSE)
    outside <- is.na(index)
    index[outside] <- findInterval(values[outside], breaks, all.inside = TRUE)
    index
  }

  x_index <- bin_index(scores[, 1], grid$x_breaks)
  y_index <- bin_index(scores[, 2], grid$y_breaks)
  histogram <- grid$base_coord
  histogram$prob <- tabulate(
    (y_index - 1L) * bins[1] + x_index,
    nbins = prod(bins)
  ) / nrow(scores)

  histogram
}

# Rebuild ordinary group histograms after sign alignment or display centering.
.dppca_reference_groups <- function(scores, group, frame, bins) {
  group_levels <- unique(as.character(group))
  groups <- lapply(group_levels, function(label) {
    rows <- which(group == label)
    list(
      n = length(rows),
      nonprivate = .dppca_reference_histogram(
        scores[rows, , drop = FALSE], frame, bins
      )
    )
  })
  stats::setNames(groups, group_levels)
}

# Stack same-grid histograms for a group-color overlay.
.dppca_group_histogram <- function(groups, method) {
  pieces <- lapply(names(groups), function(label) {
    histogram <- groups[[label]][[method]]
    histogram$group <- rep(label, nrow(histogram))
    histogram
  })
  combined <- do.call(rbind, pieces)
  combined$group <- factor(combined$group, levels = names(groups))
  combined
}

# Use reference group sizes when supplied, otherwise allocate equally.
# Zero-mass groups remain empty.
.dppca_group_sample <- function(groups, method, sampling, group_sizes = NULL) {
  if (is.null(group_sizes)) {
    group_sizes <- stats::setNames(rep(1, length(groups)), names(groups))
  } else {
    group_sizes <- group_sizes[names(groups)]
  }
  sizes <- allocate_group_sample_sizes(group_sizes, sampling$sample_size)
  pieces <- lapply(names(groups), function(label) {
    if (sizes[[label]] == 0L) {
      return(data.frame(x = numeric(), y = numeric(), group = character()))
    }
    sampled <- sample_private_score_histogram(
      groups[[label]][[method]], sizes[[label]],
      sampling$method, sampling$bandwidth_scale
    )[[1L]]
    data.frame(
      x = sampled$pc_x, y = sampled$pc_y,
      group = rep(label, nrow(sampled))
    )
  })
  do.call(rbind, pieces)
}

#' Plot an integrated PCA fit
#'
#' Uses the stored estimates without estimating loading, scree or score again.
#' Use `set.seed()` before plotting to reproduce synthetic samples.
#'
#' @param x A fitted [dp_pca()] object.
#' @param display `"all"` (default) or a vector selecting `"scree"`,
#'   `"score"`, `"loading"`, `"variable"`, and `"biplot"`.
#' @param scree Plot options: `type = "pve"` or `"scree"`, and
#'   `plot_control = NULL` ([scree_plot_control()]).
#' @param score Plot options: `type = "histogram"` or `"sample"`, and
#'   `plot_control = NULL` ([score_plot_control()]). For grouped plots, the
#'   color map and legend settings also apply to the biplot.
#' @param loading Plot options: `type = "vector"`, `"point"`, or `"heatmap"`,
#'   `heatmap_components = NULL`, and `plot_control = NULL`
#'   ([loading_plot_control()]). Multiple types require
#'   `display = "loading"`.
#' @param variable Plot options: `plot_control = NULL`. Shows labeled arrows
#'   and the correlation-style unit circle.
#' @param biplot Plot options: `alpha = 0` (between 0 and 1), `labels = TRUE`,
#'   and `plot_control = NULL`.
#' @param sampling A [sampling_control()] object with one method, shared by the
#'   score sample plot and biplot. By default, DP sample size is the input
#'   sample size, and non-private plots use observed scores. An explicit
#'   `sample_size` samples both versions from their histograms. For groups,
#'   sizes are allocated across groups; zero-mass groups yield no points.
#' @param center Whether to shift displayed private scores by the midpoint of
#'   the private frame (`TRUE` by default). Shared by score and biplot;
#'   the biplot arrows start at the origin when `TRUE`. Ordinary reference
#'   scores remain centered at zero. This is display-only post-processing.
#' @param variables Variables to display: positive top-N count (`10` by default),
#'   `"all"`, `NULL`, names, or indices. Top-N selection uses norms of private
#'   loading coefficients on the two selected axes. Selection is shared by
#'   loading, variable and biplot panels, including non-private comparisons.
#' @param base_size Base font size (`12` by default); local plot controls
#'   override it.
#' @param ... Unused arguments are rejected.
#'
#' @details
#' Plotting reuses the stored estimates; non-private panels appear
#' automatically if the fit contains them. Synthetic sampling is
#' post-processing, and `set.seed()` can reproduce sampled displays.
#'
#' Variable-plot coordinates use the covariance reconstructed from the `k`
#' private scree estimates and directions; they need not equal observed
#' variable–score correlations. Variables with zero implied variance are omitted.
#'
#' Biplot coordinates divide sampled scores by `d^alpha` and multiply private
#' direction coefficients by `d^alpha`, where
#' `d = sqrt((n - 1) * scree)` on the selected axes. A zero scree value with
#' positive alpha suppresses that axis. Arrows use one display multiplier
#' computed from the private sample and shared with the ordinary comparison;
#' empty or zero-radius samples omit arrows. The comparison uses the same square
#' biplot frame.
#'
#' Centering translates only displayed coordinates. Combined grouped pages
#' show one group legend, while individual ggplot objects retain their own.
#' The returned object stores both the combined pages and individual plots.
#'
#' @return Invisibly returns `private` and, when requested, `nonprivate` lists.
#'   Each contains `plots`, `page`, `sampled_points` (unshifted), `biplot`, `frames` and
#'   `direction_sign`. Biplot metadata includes the alpha-transformed and
#'   displayed coordinates, `scale`, `arrow_scale`, `score_center`,
#'   `score_radius`, `estimated_center`, `centered`, `arrow_start`,
#'   `arrow_endpoints` and display frame. The plots are always drawn;
#'   `out <- plot(fit)` also keeps these results for reuse. Multiple loading
#'   displays with both versions additionally return a `comparison` patchwork.
#' @export
plot.dp_pca <- function(
    x,
    display = "all",
    scree = list(),
    score = list(),
    loading = list(),
    variable = list(),
    biplot = list(),
    sampling = NULL,
    center = TRUE,
    variables = 10,
    base_size = 12,
    ...
) {
  # Resolve plot selections and local options.
  valid_types <- c("scree", "score", "loading", "variable", "biplot")

  if (
    !is.character(display) || !length(display) || anyNA(display) ||
      any(!display %in% c("all", valid_types))
  ) {
    stop("Unknown plot display.", call. = FALSE)
  }

  if ("all" %in% display && length(display) > 1L) {
    stop("Use 'all' alone, or a vector of plot names.", call. = FALSE)
  }

  if (length(list(...))) {
    stop(
      "Unknown plotting argument(s): ",
      paste(names(list(...)), collapse = ", "),
      call. = FALSE
    )
  }

  non_private <- !is.null(x$nonprivate)
  types <- if (identical(display, "all")) valid_types else unique(display)

  sp <- .dppca_options(
    scree,
    list(type = "pve", plot_control = NULL),
    "scree"
  )
  hp <- .dppca_options(
    score,
    list(type = "histogram", plot_control = NULL),
    "score"
  )
  lp <- .dppca_options(
    loading,
    list(
      type = "vector",
      heatmap_components = NULL,
      plot_control = NULL
    ),
    "loading"
  )
  vp <- .dppca_options(
    variable,
    list(plot_control = NULL),
    "variable"
  )
  bp <- .dppca_options(
    biplot,
    list(
      alpha = 0,
      labels = TRUE,
      plot_control = NULL
    ),
    "biplot"
  )

  sp$type <- match.arg(sp$type, c("pve", "scree"))
  hp$type <- match.arg(hp$type, c("histogram", "sample"))
  lp$type <- .validate_loading_display(lp$type)

  if (length(lp$type) > 1L && !identical(types, "loading")) {
    stop(
      "Multiple loading types require display = 'loading'.",
      call. = FALSE
    )
  }

  validate_logical_value(center, "center")
  validate_logical_value(bp$labels, "biplot$labels")
  .dppca_number(bp$alpha, "biplot$alpha", 0, 1)
  .dppca_number(base_size, "base_size", 0, Inf, open = TRUE)

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Install ggplot2 to plot.", call. = FALSE)
  }

  # Resolve appearance and one shared sampling specification.
  pc <- .merge_scree_plot_control(
    sp$plot_control,
    if (sp$type == "pve") "pve" else "scree"
  )
  hc <- .merge_score_plot_control(hp$plot_control)
  grouped <- !is.null(x$score_groups) && length(x$score_groups) > 0L
  group_levels <- if (grouped) names(x$score_groups) else character()
  group_colors <- NULL
  group_legend <- "none"
  if (grouped) {
    group_colors <- if (!is.null(hc$group_colors)) {
      complete_color_map(hc$group_colors, group_levels)
    } else if (.is_color_vec(group_levels)) {
      stats::setNames(group_levels, group_levels)
    } else {
      complete_color_map(NULL, group_levels)
    }
    group_legend <- if (hc$show_group_legend) "right" else "none"
  }
  lc <- .merge_loading_plot_control(lp$plot_control, x$metadata$axes)
  vc <- .dppca_vector_control(vp$plot_control, base_size)
  bc <- .dppca_vector_control(
    bp$plot_control, base_size, biplot = TRUE
  )

  if (!"base_size" %in% names(hp$plot_control)) {
    hc$base_size <- base_size
  }
  if (!"title_size" %in% names(hp$plot_control)) {
    hc$title_size <- hc$base_size + 2
  }
  if (!"base_size" %in% names(lp$plot_control)) {
    lc$base_size <- base_size
  }
  if (!"title_size" %in% names(lp$plot_control)) {
    lc$title_size <- lc$base_size + 2
  }

  heatpcs <- .validate_loading_components(
    lp$heatmap_components,
    ncol(x$directions),
    x$metadata$axes,
    "heatmap_components"
  )

  # Select once from the private directions for every panel and version.
  common_idx <- .select_loading_variables(
    variables, rownames(x$directions), x$directions, x$metadata$axes
  )
  common_names <- rownames(x$directions)[common_idx]

  if (is.null(sampling)) {
    sampling <- sampling_control(sample_size = NULL)
  }

  # NULL means original non-private points, not a synthetic sample.
  sample_nonprivate <- is.list(sampling) && !is.null(sampling$sample_size)
  sampling <- .resolve_sampling_control(sampling, x$metadata$n)

  if (length(sampling$method) != 1L) {
    stop("Choose one shared sampling method.", call. = FALSE)
  }

  need_points <- "biplot" %in% types ||
    ("score" %in% types && hp$type == "sample")
  reference_group_sizes <- if (grouped && !is.null(x$nonprivate)) {
    vapply(x$nonprivate$score_groups, function(g) g$n, numeric(1))
  } else {
    NULL
  }
  versions <- c(if (non_private) "nonprivate", "private")
  result <- list()

  # Compute the DP display first so comparison settings depend only on DP output.
  for (version in c("private", if (non_private) "nonprivate")) {
    private <- version == "private"
    obj <- if (private) x else x$nonprivate
    direction_sign <- rep(1, ncol(obj$directions))

    if (!private) {
      direction_sign <- ifelse(
        colSums(obj$directions * x$directions) < 0,
        -1,
        1
      )
      obj$directions <- sweep(obj$directions, 2, direction_sign, "*")
      obj$scores <- sweep(obj$scores, 2, direction_sign[x$metadata$axes], "*")
      obj$hist <- .dppca_reference_histogram(
        obj$scores,
        obj$frame,
        x$metadata$controls$score$bins
      )
      if (grouped) {
        obj$score_groups <- .dppca_reference_groups(
          obj$scores, obj$group, obj$frame, x$metadata$controls$score$bins
        )
      }
      # Preserve the centered ordinary frame; never overwrite with the DP frame.
    }

    ax <- x$metadata$axes
    V <- obj$directions
    ev <- obj$scree
    prefix <- if (private) "DP" else "Non-private"
    page_title <- if (private) "DP-PCA" else "Non-private PCA"

    # One center per version, shared by all groups.
    estimated_center <- if (private) {
      c(mean(x$frame$xlim), mean(x$frame$ylim))
    } else {
      c(0, 0)  # Ordinary scores were already centered before histogramming.
    }
    score_shift <- if (center) estimated_center else c(0, 0)
    biplot_shift <- score_shift

    # Move display copies only; sample from the original histogram below.
    score_frame <- obj$frame
    score_frame$xlim <- score_frame$xlim - score_shift[1]
    score_frame$ylim <- score_frame$ylim - score_shift[2]
    plot_frames <- list(score = score_frame)
    points <- data.frame(x = numeric(), y = numeric())

    if (need_points) {
      if (grouped) {
        if (private || sample_nonprivate) {
          points <- .dppca_group_sample(
            obj$score_groups,
            if (private) x$metadata$methods$score else "nonprivate",
            sampling,
            group_sizes = reference_group_sizes
          )
        } else {
          points <- data.frame(
            x = obj$scores[, 1], y = obj$scores[, 2], group = obj$group
          )
        }
        points$group <- factor(points$group, levels = group_levels)
      } else if (private || sample_nonprivate) {
        sampled <- sample_private_score_histogram(
          obj$hist,
          sampling$sample_size,
          sampling$method,
          sampling$bandwidth_scale
        )[[1]]
        points <- data.frame(x = sampled$pc_x, y = sampled$pc_y)
      } else {
        points <- data.frame(x = obj$scores[, 1], y = obj$scores[, 2])
      }
    }

    plots <- list()
    coords <- NULL

    for (tp in types) {
      note <- NULL
      panel_base <- base_size
      panel_title_size <- base_size + 2
      style <- if (tp == "biplot") bc else vc
      color <- if (private) style$private_color else style$nonprivate_color

      if (tp %in% c("variable", "biplot")) {
        panel_base <- style$base_size
        panel_title_size <- style$title_size
      } else if (tp == "score") {
        panel_base <- hc$base_size
        panel_title_size <- hc$title_size
      }

      if (tp == "scree") {
        values <- if (sp$type == "pve") scree_to_pve(ev) else ev
        idx <- seq_along(ev)
        dat <- data.frame(pc = idx, value = values[idx])
        series <- if (private) x$metadata$methods$scree else "nonprivate"

        pl <- ggplot2::ggplot(dat, ggplot2::aes(x = .data$pc, y = .data$value)) +
          ggplot2::geom_line(
            color = unname(pc$col[series]),
            linewidth = pc$lwd * 0.4,
            linetype = unname(pc$lty[series])
          ) +
          ggplot2::geom_point(
            color = unname(pc$col[series]),
            shape = unname(pc$pch[series]),
            size = 2 * pc$point_cex
          ) +
          ggplot2::scale_x_continuous(breaks = idx) +
          ggplot2::labs(
            x = "Component",
            y = if (sp$type == "pve") {
              "PVE"
            } else {
              "Eigenvalue"
            }
          )

        if (sum(ev) <= 0) {
          note <- "All estimated eigenvalues are zero"
        }
      } else if (tp == "score") {
        if (hp$type == "histogram") {
          histogram <- if (grouped) {
            .dppca_group_histogram(
              obj$score_groups,
              if (private) x$metadata$methods$score else "nonprivate"
            )
          } else {
            obj$hist
          }
          histogram$xmin <- histogram$xmin - score_shift[1]
          histogram$xmax <- histogram$xmax - score_shift[1]
          histogram$ymin <- histogram$ymin - score_shift[2]
          histogram$ymax <- histogram$ymax - score_shift[2]
          if (grouped) {
            pl <- make_hist_all_dp(
              histogram, score_frame$xlim, score_frame$ylim, group_colors,
              alpha_range = hc$hist_alpha_range, base_size = hc$base_size,
              title_size = hc$title_size, legend_position = group_legend
            )
          } else {
            pl <- make_hist_plot_dp(
              histogram,
              score_frame$xlim,
              score_frame$ylim,
              hc$color,
              alpha_range = hc$hist_alpha_range,
              base_size = hc$base_size,
              title_size = hc$title_size
            )
          }

          # Each histogram retains its own probability-to-opacity scale.
          # make_hist_plot_dp()/make_hist_all_dp() train alpha on that panel.

          if (sum(obj$hist$prob) <= 0) {
            note <- "No bins survived private thresholding"
          }
        } else {
          if (grouped) {
            pl <- make_group_sample_plot_dp(
              data.frame(
                pc_x = points$x - score_shift[1],
                pc_y = points$y - score_shift[2],
                group = points$group
              ),
              score_frame$xlim, score_frame$ylim, group_colors,
              point_alpha = hc$scatter_alpha, point_size = hc$scatter_size,
              base_size = hc$base_size, title_size = hc$title_size,
              legend_position = group_legend
            )
          } else {
            pl <- make_sample_plot_dp(
              data.frame(pc_x = points$x - score_shift[1],
                         pc_y = points$y - score_shift[2]),
              score_frame$xlim,
              score_frame$ylim,
              hc$color,
              point_alpha = hc$scatter_alpha,
              point_size = hc$scatter_size,
              base_size = hc$base_size,
              title_size = hc$title_size
            )
          }

          if (!nrow(points)) {
            note <- "No synthetic points: histogram has zero mass"
          }
        }
      } else if (tp == "loading") {
        for (loading_type in lp$type) {
          components <- if (loading_type == "heatmap") heatpcs else ax
          selected <- .select_loading_variables(
            common_names,
            rownames(V),
            V,
            components
          )
          loading_title <- paste(
            prefix,
            if (loading_type == "heatmap") "Loading Heatmap" else "Loading Plot"
          )
          title_key <- if (private) "private_title" else "nonprivate_title"

          if (
            title_key %in% names(lp$plot_control) &&
              (is.null(lc[[title_key]]) ||
                !lc[[title_key]] %in% c("DP", "Non-private", "Non-Private"))
          ) {
            loading_title <- lc[[title_key]]
          }

          if (loading_type == "heatmap") {
            dat <- .make_loading_heatmap_data(
              V, components, selected, rownames(V)
            )
            fill_limit <- max(abs(dat$loading))
            if (!is.finite(fill_limit) ||
                fill_limit <= sqrt(.Machine$double.eps)) {
              fill_limit <- 1
            }
            pl <- .build_loading_heatmap(
              dat, loading_title, fill_limit, lc
            )
          } else {
            dat <- .make_loading_data(
              V, ax, selected, rownames(V)
            )
            pl <- .build_loading_panel(
              dat,
              if (private) lc$private_color else lc$nonprivate_color,
              loading_title,
              lc$xlab,
              lc$ylab,
              loading_type,
              lc
            )
          }

          key <- if (length(lp$type) == 1L) {
            "loading"
          } else {
            paste0("loading_", loading_type)
          }
          plots[[key]] <- pl
        }
        next
      } else {
        # Variable coordinates and biplot coordinates share the arrow renderer.
        if (tp == "variable") {
          weighted <- sweep(
            V[, seq_along(ev), drop = FALSE],
            2,
            sqrt(pmax(ev, 0)),
            "*"
          )
          denominator <- if (private) {
            sqrt(rowSums(weighted^2))
          } else {
            obj$variable_sd
          }
          A <- weighted[, ax, drop = FALSE]
          good <- is.finite(denominator) & denominator > 0
          A[good, ] <- sweep(
            A[good, , drop = FALSE], 1, denominator[good], "/"
          )
          A[!good, ] <- NA_real_

          if (any(!good)) {
            note <- paste(
              c(note, paste(sum(!good), "zero-variance variables omitted")),
              collapse = "; "
            )
          }
          pl <- ggplot2::ggplot()
        } else {
          # Biplot and score plots use the same display-centering option.
          biplot_points <- points[, c("x", "y"), drop = FALSE]
          biplot_points$x <- biplot_points$x - biplot_shift[1]
          biplot_points$y <- biplot_points$y - biplot_shift[2]
          biplot_frame <- obj$frame
          biplot_frame$xlim <- biplot_frame$xlim - biplot_shift[1]
          biplot_frame$ylim <- biplot_frame$ylim - biplot_shift[2]
          coords <- .dppca_biplot_coordinates(
            points = biplot_points,
            V = V,
            eigenvalues = ev,
            n = x$metadata$n,
            axes = ax,
            alpha = bp$alpha,
            frame = biplot_frame,
            arrow_scale = if (private) NULL else result$private$biplot$arrow_scale
          )

          # Uncentered scores keep their positions; only the arrows move.
          arrow_start <- c(0, 0)
          if (!center) {
            active <- coords$scale > 0
            arrow_start[active] <- estimated_center[active] / coords$scale[active]
          }
          coords$estimated_center <- stats::setNames(estimated_center, c("x", "y"))
          coords$centered <- center
          coords$arrow_start <- stats::setNames(arrow_start, c("x", "y"))
          coords$arrow_endpoints <- sweep(
            coords$display_variables, 2, arrow_start, "+"
          )
          A <- coords$display_variables
          note <- coords$note
          if (grouped) {
            coords$observations$group <- points$group
            coords$display_observations$group <- points$group
          }
          pl <- ggplot2::ggplot(
            coords$display_observations,
            ggplot2::aes(x = .data$x, y = .data$y)
          )
          if (grouped) {
            pl <- pl + ggplot2::geom_point(
              ggplot2::aes(colour = .data$group),
              alpha = bc$point_alpha, size = bc$point_size
            ) + ggplot2::scale_colour_manual(
                values = group_colors,
                breaks = group_levels[group_levels %in% as.character(points$group)],
                drop = TRUE
              )
          } else {
            pl <- pl + ggplot2::geom_point(
              color = bc$point_color,
              alpha = bc$point_alpha,
              size = bc$point_size
            )
          }
        }

        keep <- .dppca_variables(A, rownames(V), common_names)
        if (tp == "biplot" && coords$arrow_scale == 0) {
          keep <- integer()
        }
        dat <- data.frame(
          x = A[keep, 1],
          y = A[keep, 2],
          label = rownames(V)[keep]
        )

        # Variable plots keep origin-based arrows; biplots use their center.
        arrow_mapping <- ggplot2::aes(x = 0, y = 0, xend = .data$x, yend = .data$y)
        if (tp == "biplot") {
          dat$start_x <- rep(arrow_start[1], nrow(dat))
          dat$start_y <- rep(arrow_start[2], nrow(dat))
          dat$x <- dat$x + arrow_start[1]
          dat$y <- dat$y + arrow_start[2]
          arrow_mapping <- ggplot2::aes(
            x = .data$start_x, y = .data$start_y,
            xend = .data$x, yend = .data$y
          )
        }

        pl <- pl +
          ggplot2::geom_hline(
            yintercept = 0, color = "grey75", linewidth = 0.4,
            linetype = "dashed"
          ) +
          ggplot2::geom_vline(
            xintercept = 0, color = "grey75", linewidth = 0.4,
            linetype = "dashed"
          )

        pl <- pl + ggplot2::geom_segment(
          data = dat,
          mapping = arrow_mapping,
          inherit.aes = FALSE,
          color = color,
          linewidth = style$arrow_size,
          arrow = grid::arrow(
            length = grid::unit(0.18, "cm"), type = "closed"
          )
        )

        if ((tp == "variable" || bp$labels) && nrow(dat)) {
          if (requireNamespace("ggrepel", quietly = TRUE)) {
            pl <- pl + ggrepel::geom_text_repel(
              data = dat,
              ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
              inherit.aes = FALSE,
              color = color,
              size = style$label_size,
              max.overlaps = Inf,
              box.padding = 0.4,
              point.padding = 0.2,
              segment.color = "grey60"
            )
          } else {
            pl <- pl + ggplot2::geom_text(
              data = dat,
              ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
              inherit.aes = FALSE,
              color = color,
              size = style$label_size,
              vjust = -0.6
            )
          }
        }

        if (tp == "variable") {
          angle <- seq(0, 2 * pi, length.out = 361)
          circle <- data.frame(x = cos(angle), y = sin(angle))
          pl <- pl +
            ggplot2::geom_path(
              data = circle,
              ggplot2::aes(x = .data$x, y = .data$y),
              inherit.aes = FALSE,
              color = "grey65",
              linetype = 2
            ) +
            ggplot2::coord_equal(
              xlim = c(-1.15, 1.15), ylim = c(-1.15, 1.15)
            )
        } else {
          biplot_frame <- if (private) {
            .dppca_biplot_frame(coords, arrow_start)
          } else {
            result$private$frames$biplot
          }
          coords$frame <- biplot_frame
          plot_frames$biplot <- biplot_frame

          outside <- sum(
            dat$x < biplot_frame$xlim[1] | dat$x > biplot_frame$xlim[2] |
              dat$y < biplot_frame$ylim[1] | dat$y > biplot_frame$ylim[2]
          )
          if (outside > 0) {
            note <- paste(
              c(
                note,
                paste(outside, "arrow endpoints outside the DP comparison frame")
              ),
              collapse = "; "
            )
          }

          pl <- pl +
            ggplot2::coord_equal(
              xlim = biplot_frame$xlim,
              ylim = biplot_frame$ylim,
              expand = FALSE
            ) +
            ggplot2::scale_x_continuous(
              breaks = pretty(biplot_frame$xlim, n = 5)
            ) +
            ggplot2::scale_y_continuous(
              breaks = pretty(biplot_frame$ylim, n = 5)
            )
        }
      }

      # Apply automatic labels first, then non-NULL local overrides.
      if (tp != "scree") {
        pl <- pl + ggplot2::labs(
          x = paste0("PC", ax[1]), y = paste0("PC", ax[2])
        )
      }
      score_title <- if (hp$type == "histogram") {
        paste(prefix, "Histogram")
      } else if (private || sample_nonprivate) {
        paste(prefix, "Sample")
      } else {
        "Non-private Scatter"
      }
      panel_title <- switch(
        tp,
        scree = paste(prefix, "Scree Plot"),
        score = score_title,
        variable = paste(prefix, "Variable Plot"),
        biplot = paste(prefix, "Biplot")
      )

      if (tp != "score") {
        pl <- pl +
          ggplot2::theme_classic(base_size = panel_base) +
          ggplot2::theme(
            panel.border = ggplot2::element_rect(
              colour = if (tp == "scree") "grey50" else "grey80",
              fill = NA,
              linewidth = 0.5
            ),
            axis.line = ggplot2::element_blank()
          )
      }
      if (tp == "biplot") {
        pl <- pl + theme_dp_base(panel_base)
      }

      result_note <- !is.null(note) &&
        grepl("zero|No |no synthetic|omitted|outside", note)
      subtitle <- if (result_note) note else NULL
      pl <- pl +
        ggplot2::labs(title = panel_title, subtitle = subtitle) +
        ggplot2::theme(
          plot.title = ggplot2::element_text(
            face = "bold", size = panel_title_size, hjust = 0.5
          ),
          plot.subtitle = ggplot2::element_text(
            size = panel_base * 0.75, hjust = 0.5
          ),
          plot.caption = ggplot2::element_text(
            size = panel_base * 0.7, hjust = 0.5
          ),
          legend.position = "none",
          aspect.ratio = if (tp == "scree") 1 else NULL,
          plot.margin = ggplot2::margin(8, 12, 8, 12)
        )
      if (grouped && tp %in% c("score", "biplot")) {
        pl <- pl + ggplot2::theme(
          legend.position = if (hc$show_group_legend) "right" else "none",
          legend.text = ggplot2::element_text(size = max(8, 0.8 * panel_base)),
          legend.key.height = grid::unit(3, "mm"),
          legend.key.width = grid::unit(3, "mm")
        ) + ggplot2::guides(
          colour = ggplot2::guide_legend(ncol = 1L),
          fill = ggplot2::guide_legend(ncol = 1L)
        )
        # Group colours identify the groups; omit the often-unhelpful
        # '<supplied vector>' title (and the 'Group' title).
        pl <- pl + ggplot2::labs(fill = NULL, colour = NULL)
      }

      override <- switch(
        tp,
        scree = sp$plot_control,
        score = hp$plot_control,
        loading = lp$plot_control,
        variable = vp$plot_control,
        biplot = bp$plot_control
      )
      label_options <- override[
        intersect(names(override), c("title", "xlab", "ylab"))
      ]
      label_options <- Filter(Negate(is.null), label_options)

      if (length(label_options)) {
        names(label_options)[names(label_options) == "xlab"] <- "x"
        names(label_options)[names(label_options) == "ylab"] <- "y"
        pl <- pl + do.call(ggplot2::labs, label_options)
      }
      if (tp == "score") {
        title_key <- if (!private && hp$type == "sample") {
          "scatter_title"
        } else if (!private) {
          "nonprivate_title"
        } else {
          paste0(x$metadata$methods$score, "_title")
        }
        custom_title <- override[[title_key]]
        default_titles <- c(
          "Non-private Scatter", "Non-private Histogram",
          "Sparse DP Histogram", "Add DP Histogram"
        )
        if (!is.null(custom_title) && !custom_title %in% default_titles) {
          pl <- pl + ggplot2::labs(title = custom_title)
        }
      }
      if (tp == "scree") {
        pl <- pl + ggplot2::theme(
          # Normalize the default 1.5 multiplier to match other panel titles.
          plot.title = ggplot2::element_text(
            size = (base_size + 2) * pc$cex_main / 1.5
          ),
          axis.title = ggplot2::element_text(size = base_size * pc$cex_lab),
          axis.text = ggplot2::element_text(size = base_size * pc$cex_axis)
        )
      }
      plots[[tp]] <- pl
    }

    heatmap_name <- if (
      "loading" %in% types && identical(lp$type, "heatmap")
    ) {
      "loading"
    } else {
      NULL
    }
    page <- .dppca_page(
      plots, page_title, heatmap_name = heatmap_name, base_size = base_size,
      shared_group_legend = grouped && hc$show_group_legend,
      group_count = length(group_levels)
    )
    result[[version]] <- list(
      plots = plots,
      page = page,
      sampled_points = points,
      biplot = coords,
      frames = plot_frames,
      direction_sign = stats::setNames(direction_sign, colnames(V))
    )
  }

  # Drawing uses the familiar non-private -> private order.
  result <- result[versions]

  if (identical(types, "loading") && length(lp$type) > 1L && non_private) {
    if (!requireNamespace("patchwork", quietly = TRUE)) {
      stop("Install patchwork for loading comparison.", call. = FALSE)
    }
    result$comparison <- patchwork::wrap_plots(
      c(unname(result$nonprivate$plots), unname(result$private$plots)),
      ncol = length(lp$type),
      nrow = 2,
      byrow = TRUE,
      guides = "collect"
    ) & ggplot2::theme(legend.position = "right")
    print(result$comparison)
  } else {
    for (version in versions) {
      grid::grid.newpage()
      grid::grid.draw(result[[version]]$page)
    }
  }

  invisible(result)
}

# Transform once; the DP multiplier is reused for the non-private comparison.
.dppca_biplot_coordinates <- function(
    points,
    V,
    eigenvalues,
    n,
    axes,
    alpha,
    frame,
    arrow_scale = NULL
) {
  d <- sqrt((n - 1) * pmax(eigenvalues[axes], 0))
  scale <- if (alpha == 0) rep(1, 2) else d^alpha
  active <- scale > 0
  observations <- as.matrix(points)
  H <- V[, axes, drop = FALSE]

  if (any(active)) {
    observations[, active] <- sweep(
      observations[, active, drop = FALSE], 2, scale[active], "/"
    )
    H[, active] <- sweep(H[, active, drop = FALSE], 2, scale[active], "*")
  }
  if (any(!active)) {
    observations[, !active] <- 0
    H[, !active] <- 0
  }

  if (any(!is.finite(observations)) || any(!is.finite(H))) {
    stop("Non-finite alpha-transformed biplot coordinates.", call. = FALSE)
  }

  score_center <- c(x = NA_real_, y = NA_real_)
  score_radius <- 0

  if (nrow(observations)) {
    score_center <- colMeans(observations)
    centered <- sweep(observations, 2, score_center, "-")
    score_radius <- as.numeric(stats::quantile(
      sqrt(rowSums(centered^2)), probs = 0.9, names = FALSE
    ))
  }

  # All variables determine the multiplier, before selecting plotted labels.
  arrow_radius <- max(sqrt(rowSums(H^2)))

  if (any(!is.finite(c(score_radius, arrow_radius)))) {
    stop("Biplot radii exceed the finite numeric range.", call. = FALSE)
  }

  if (is.null(arrow_scale)) {
    arrow_scale <- if (score_radius > 0 && arrow_radius > 0) {
      score_radius / arrow_radius
    } else {
      0
    }
  }

  displayed_variables <- H * arrow_scale
  observations <- data.frame(x = observations[, 1], y = observations[, 2])
  factor <- numeric(2)
  factor[active] <- 1 / scale[active]
  transformed_frame <- list(
    xlim = frame$xlim * factor[1],
    ylim = frame$ylim * factor[2]
  )

  if (
    !is.finite(arrow_scale) || any(!is.finite(displayed_variables)) ||
      any(!is.finite(unlist(transformed_frame)))
  ) {
    stop("Biplot scaling exceeds the finite numeric range.", call. = FALSE)
  }

  note <- NULL

  if (any(!active)) {
    note <- paste(c(note, "zero-eigenvalue axis suppressed"), collapse = "; ")
  }
  if (!nrow(observations)) {
    note <- paste(c(note, "no synthetic observations"), collapse = "; ")
  }
  if (arrow_scale == 0) {
    note <- paste(c(note, "zero reference radius; arrows omitted"), collapse = "; ")
  }

  list(
    observations = observations,
    variables = H,
    scale = scale,
    display_observations = observations,
    display_variables = displayed_variables,
    arrow_scale = arrow_scale,
    score_center = score_center,
    score_radius = score_radius,
    arrow_radius = arrow_radius,
    transformed_score_frame = transformed_frame,
    note = note
  )
}

# Include the origin, arrow start and DP endpoints in a square display.
.dppca_biplot_frame <- function(
    coordinates,
    arrow_start = c(0, 0)
) {
  frame <- coordinates$transformed_score_frame
  arrows <- sweep(coordinates$display_variables, 2, arrow_start, "+")
  frame$xlim <- range(frame$xlim, 0, arrow_start[1], arrows[, 1])
  frame$ylim <- range(frame$ylim, 0, arrow_start[2], arrows[, 2])

  .dppca_square_frame(frame)
}

#' @param x A `dp_pca` object.
#' @param digits Significant digits for printing, from `1` to `22`.
#' @param ... Additional formatting arguments passed to [print()] by the print
#'   method.
#' @rdname dp_pca
#' @export
print.dp_pca <- function(
    x,
    digits = max(3L, getOption("digits") - 3L),
    ...
) {
  .dppca_number(digits, "digits", 1, 22, integer = TRUE)

  cat("Differentially Private PCA\n")
  cat("Data:", x$metadata$n, "observations x", nrow(x$directions), "variables\n")
  cat(
    "Components:", ncol(x$directions), "PC directions;",
    x$metadata$k, "scree estimates\n"
  )
  cat("Score axes:", paste0("PC", x$metadata$axes, collapse = ", "), "\n")
  if (!is.null(x$metadata$group_name)) {
    cat("Group variable:", x$metadata$group_name, "\n")
    cat("Groups:", paste(x$metadata$group_levels, collapse = ", "), "\n")
  }
  cat(
    "Methods: scree =", x$metadata$methods$scree,
    "; score =", x$metadata$methods$score, "\n"
  )
  cat(
    "Total privacy budget: eps =", format(x$metadata$privacy$total[["eps"]], digits = digits),
    "; delta =", format(x$metadata$privacy$total[["delta"]], digits = digits), "\n"
  )
  cat(
    "Non-private reference:",
    if (is.null(x$nonprivate)) "not stored" else "included", "\n"
  )


  invisible(x)
}

#' Summarize stored DP PCA estimates
#'
#' Summarizes an existing fit without new estimation or random sampling.
#'
#' @param object A [dp_pca()] result.
#' @param ... Additional arguments. The summary calculation takes no extra
#'   options; print methods pass formatting arguments to the printed tables.
#'
#' @return A `summary.dp_pca` object. `importance` contains DP eigenvalues,
#'   PVE and cumulative PVE, with a column for each estimated component.
#'   `budget` contains the allocated eps and delta for loading, scree,
#'   score frame and score histogram, plus the total. Other fields describe
#'   data dimensions, selected axes, methods and whether a non-private
#'   reference is stored. Raw scores and loading matrices are
#'   not copied into the summary.
#'   The print method returns its input invisibly.
#'
#' @details
#' PVE is normalized over the `k` retained components; a total of zero
#' produces zero PVE. When `k` is smaller than the number of variables, a
#' cumulative PVE of one does not describe all original variation.
#'
#' The budget table reports allocated budgets, not cumulative privacy spending.
#' With supplied `V_dp`, the loading allocation records its prior cost. The
#' optional non-private comparison is excluded from the importance table.
#'
#' @rdname dp_pca
#' @export
summary.dp_pca <- function(
    object,
    ...
) {
  if (length(list(...))) {
    stop("summary.dp_pca() takes no additional options.", call. = FALSE)
  }

  eigenvalues <- object$scree
  pve <- scree_to_pve(eigenvalues)
  importance <- rbind(
    Eigenvalue = eigenvalues,
    PVE = pve,
    `Cumulative PVE` = cumsum(pve)
  )
  colnames(importance) <- names(eigenvalues)

  budgets <- object$metadata$privacy$budgets
  score_budget <- object$metadata$score_budget
  budget <- rbind(
    Loading = budgets$loading,
    Scree = budgets$scree,
    `Score frame` = c(eps = score_budget$eps_frame, delta = 0),
    `Score histogram` = c(
      eps = score_budget$eps_hist,
      delta = score_budget$delta_hist
    ),
    Total = object$metadata$privacy$total
  )

  structure(
    list(
      n = object$metadata$n,
      p = nrow(object$directions),
      k = object$metadata$k,
      axes = object$metadata$axes,
      group_name = object$metadata$group_name,
      group_levels = object$metadata$group_levels,
      excluded_columns = object$metadata$excluded_columns,
      methods = object$metadata$methods,
      importance = importance,
      budget = budget,
      nonprivate = !is.null(object$nonprivate)
    ),
    class = "summary.dp_pca"
  )
}

#' @param x A `summary.dp_pca` object.
#' @param digits Significant digits used when printing tables.
#' @rdname dp_pca
#' @export
print.summary.dp_pca <- function(
    x,
    digits = max(3L, getOption("digits") - 3L),
    ...
) {
  .dppca_number(digits, "digits", 1, 22, integer = TRUE)

  cat("DP PCA Summary\n")
  cat("Data:", x$n, "observations x", x$p, "variables\n")
  cat("Components:", x$p, "PC directions;", x$k, "scree estimates\n")
  cat("Score axes:", paste0("PC", x$axes, collapse = ", "), "\n")
  cat(
    "Methods: scree =", x$methods$scree,
    "; score =", x$methods$score, "\n"
  )
  if (!is.null(x$group_name)) {
    cat("Group variable:", x$group_name, "\n")
    cat("Groups:", paste(x$group_levels, collapse = ", "), "\n")
  }

  cat("\nDP component importance:\n")
  print(x$importance, digits = digits, ...)
  cat("\nProportions are normalized over", x$k, "estimated components.\n")

  if (sum(x$importance["Eigenvalue", ]) <= 0) {
    cat("All estimated eigenvalues are zero; proportions are shown as zero.\n")
  }

  cat("\nAllocated privacy budget:\n")
  print(x$budget, digits = digits, ...)
  cat(
    "\nNon-private reference:",
    if (x$nonprivate) "included" else "not stored", "\n"
  )

  invisible(x)
}

# Helpers -------------------------------------------------------------------


.dppca_number <- function(
    x,
    name,
    lo,
    hi,
    open = FALSE,
    integer = FALSE
) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x)) {
    stop("Invalid ", name, ".")
  }

  outside <- if (open) {
    x <= lo || x >= hi
  } else {
    x < lo || x > hi
  }

  if (outside || (integer && x != floor(x))) {
    stop("Invalid ", name, ".")
  }
}

.dppca_variables <- function(
    A,
    nm,
    variables
) {
  good <- which(rowSums(!is.finite(A)) == 0)
  B <- A
  B[!is.finite(B)] <- 0

  selected <- .select_loading_variables(
    variables,
    nm,
    B,
    seq_len(ncol(B))
  )

  intersect(selected, good)
}

# Five plots with a loading heatmap: four panels left, tall heatmap right.
# Heatmap legend widths must never be copied to the other four plots.
.dppca_page <- function(
    plots,
    title,
    heatmap_name = NULL,
    base_size = 12,
    shared_group_legend = FALSE,
    group_count = 0L
) {
  n <- length(plots)
  tall_heatmap <- n == 5L && !is.null(heatmap_name) &&
    heatmap_name %in% names(plots)

  # Keep the group guide attached to the biplot on its right, in one column.
  # Suppress the duplicated score guide only in the composed page; the
  # original individual plots (returned by plot.dp_pca) retain their guides.
  page_plots <- plots
  if (shared_group_legend &&
      all(c("score", "biplot") %in% names(page_plots))) {
    page_plots$score <- page_plots$score +
      ggplot2::theme(legend.position = "none")
  }

  if (tall_heatmap) {
    if (!requireNamespace("patchwork", quietly = TRUE)) {
      stop("Install patchwork to draw the five-panel heatmap layout.",
           call. = FALSE)
    }
    page <- patchwork::wrap_plots(
      A = page_plots$scree,
      B = page_plots$score,
      C = page_plots$variable,
      D = page_plots$biplot,
      E = page_plots[[heatmap_name]],
      design = "ABE\nCDE"
    ) + patchwork::plot_annotation(
      title = title,
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(
          size = base_size + 3, face = "bold", hjust = 0.5
        )
      )
    )
    return(grid::grobTree(patchwork::patchworkGrob(page)))
  }

  nr <- if (n <= 3L) 1L else 2L
  nc <- if (n == 4L) 2L else min(n, 3L)
  # Top row: three equal plots. Bottom row: variable and a wider biplot,
  # whose one-column legend is immediately to the right of its plot panel.
  layout_cols <- if (n == 5L) 6L else nc
  layout <- grid::grid.layout(
    nr + 1L, layout_cols,
    widths = grid::unit(rep(1, layout_cols), "null"),
    heights = grid::unit(
      c(0.35, rep(1, nr)),
      c("inches", rep("null", nr))
    )
  )
  children <- list(
    grid::textGrob(
      title,
      gp = grid::gpar(fontsize = base_size + 3, fontface = "bold"),
      vp = grid::viewport(
        layout.pos.row = 1L, layout.pos.col = seq_len(layout_cols)
      )
    )
  )
  grobs <- lapply(page_plots, ggplot2::ggplotGrob)
  # Do not align heatmap colourbar widths with ordinary plot panels.
  align <- if (!is.null(heatmap_name)) {
    which(names(page_plots) != heatmap_name)
  } else {
    seq_along(grobs)
  }
  if (length(align) > 1L) {
    width_counts <- vapply(grobs[align], function(g) length(g$widths), integer(1))
    if (length(unique(width_counts)) == 1L) {
      # A guide belongs only to its own biplot and should not make every
      # panel reserve the same guide width.
      unguided <- setdiff(align, which(names(page_plots) == "biplot"))
      if (length(unguided) > 1L) {
        aligned_widths <- do.call(
          grid::unit.pmax,
          lapply(grobs[unguided], function(g) g$widths)
        )
        for (i in unguided) grobs[[i]]$widths <- aligned_widths
      }
    }
  }
  for (i in seq_along(page_plots)) {
    if (n == 5L) {
      row <- if (i <= 3L) 2L else 3L
      col <- switch(i, 1:2, 3:4, 5:6, 1:3, 4:6)
    } else {
      row <- 2L + (i - 1L) %/% nc
      col <- 1L + (i - 1L) %% nc
    }
    children[[length(children) + 1L]] <- grid::grobTree(
      grobs[[i]],
      vp = grid::viewport(layout.pos.row = row, layout.pos.col = col)
    )
  }
  grid::gTree(
    children = do.call(grid::gList, children),
    vp = grid::viewport(layout = layout)
  )
}

# Equal x/y spans give a square panel under coord_equal / coord_fixed.
# Only pad the shorter dimension; never rescale the data or force aspect.ratio.
.dppca_square_frame <- function(
    frame
) {
  xl <- as.numeric(frame$xlim)
  yl <- as.numeric(frame$ylim)

  if (
    length(xl) != 2L || length(yl) != 2L ||
      any(!is.finite(c(xl, yl))) || diff(xl) < 0 || diff(yl) < 0
  ) {
    stop("Invalid plot frame.")
  }

  span <- max(diff(xl), diff(yl))
  if (span <= 0) {
    span <- 1
  }

  list(
    xlim = mean(xl) + c(-0.5, 0.5) * span,
    ylim = mean(yl) + c(-0.5, 0.5) * span
  )
}

# Merge only recognized fields, preserving nested control objects and NULL.
.dppca_options <- function(
    x,
    defaults,
    name
) {
  if (
    !is.list(x) || (length(x) && (
      is.null(names(x)) || anyNA(names(x)) || any(names(x) == "") ||
        anyDuplicated(names(x)) || any(!names(x) %in% names(defaults))
    ))
  ) {
    stop(
      name, " must be a named list. Allowed fields: ",
      paste(names(defaults), collapse = ", ")
    )
  }

  for (nm in names(x)) {
    defaults[nm] <- x[nm]
  }

  defaults
}

# Variable/biplot appearance controls (plain named lists; no extra constructor).
# private_color/nonprivate_color style variable points, arrows, and their labels.
# point_color styles the synthetic observation points in biplots.
.dppca_vector_control <- function(
    control,
    base_size,
    biplot = FALSE
) {
  if (is.null(control)) {
    control <- list()
  }

  out <- .dppca_options(
    control,
    list(
      nonprivate_color = if (biplot) "#174A7E" else "black",
      private_color = if (biplot) "#174A7E" else "#1F4E79",
      arrow_size = 0.7,
      point_size = 1.8,
      point_color = if (biplot) "#B0B0B0" else "#6A5ACD",
      point_alpha = 0.6,
      label_size = 4,
      base_size = base_size,
      title_size = base_size + 2,
      title = NULL,
      xlab = NULL,
      ylab = NULL
    ),
    "variable/biplot plot_control"
  )

  .dppca_number(out$base_size, "base_size", 0, Inf, open = TRUE)

  if (!"title_size" %in% names(control)) {
    out$title_size <- out$base_size + 2
  }

  for (nm in c("nonprivate_color", "private_color", "point_color")) {
    .validate_plot_color(out[[nm]], nm)
  }

  for (nm in c("arrow_size", "point_size", "label_size", "title_size")) {
    .dppca_number(out[[nm]], nm, 0, Inf, open = TRUE)
  }

  .dppca_number(out$point_alpha, "point_alpha", 0, 1)

  for (nm in c("title", "xlab", "ylab")) {
    .validate_plot_optional_string(out[[nm]], nm)
  }

  out
}
