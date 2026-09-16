# dp_pca.R

#' Estimate shared PCA quantities once
#'
#' @param X Numeric matrix or data frame with observations in rows.
#' @param privacy Output of [privacy_control()] with loading, scree and score
#'   budgets. The default split allocates one third of eps and delta to each.
#' @param scree Required named list containing a method-specific `control`.
#'   Other fields are `k = NULL` (all variables), `method = "clipped"` and
#'   `mono = TRUE`. Use [clipped_control()], [pmwm_control()] or
#'   [huber_control()] to supply the selected method's tuning parameters.
#' @param score Named list with `method = "sparse"` and `bins = c(20, 20)`.
#'   The other supported histogram method is `"add"`.
#' @param axes Two distinct histogram component indices between `1` and `k`.
#' @param center,standardize Whether to center or standardize the columns.
#' @param cpp.option Whether to use the compiled loading implementation.
#' @param non_private Whether to store ordinary PCA comparison results.
#'
#' @return A `dp_pca` object containing all `directions`, the first `k`
#'   `eigenvalues` and their `pve`, `score_histogram`, `frame`, `privacy`,
#'   `score_budget`, resolved `methods`, `controls`, `preprocessing`, `n`,
#'   `k` and `axes`. `nonprivate` contains the optional reference results,
#'   including observed scores; it is `NULL` when `non_private = FALSE`.
#'
#' @details
#' One private loading estimate is shared by scree and score estimation.
#' The selected scree method receives the full scree allocation and divides
#' it across the first `k` components. Within the score allocation, 35 percent
#' of `eps` is used for the frame and 65 percent for the histogram; all score
#' `delta` is assigned to the histogram. PVE is normalized over the `k`
#' estimated values. Non-private references are comparison information and
#' are not private releases. Budget allocation alone does not establish a
#' privacy guarantee for the complete preprocessing and estimation procedure.
#'
#' @examples
#' data(gau, package = "dppca")
#' set.seed(123)
#' fit <- dp_pca(
#'   gau[1:200, ],
#'   privacy = privacy_control(eps = 5, delta = 1e-4),
#'   scree = list(k = 2, control = clipped_control(C_clip = 3)),
#'   score = list(method = "sparse", bins = c(10, 10))
#' )
#' fit
#' summary(fit)
#' fit$directions
#'
#' set.seed(123)
#' plots <- plot(fit)
#' plots$private$plots$biplot
#'
#' @export
dp_pca <- function(
    X,
    privacy,
    scree,
    score = list(),
    axes = c(1, 2),
    center = TRUE,
    standardize = FALSE,
    cpp.option = TRUE,
    non_private = TRUE
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

  for (argument in c("center", "standardize", "cpp.option", "mono", "non_private")) {
    .dppca_flag(get(argument), argument)
  }

  privacy <- .dppca_fun(".validate_privacy_control")(
    privacy, c("loading", "scree", "score")
  )
  budgets <- privacy$components

  control <- .dppca_fun(".merge_scree_control")(
    scree_method,
    scree_options$control
  )
  control <- do.call(.dppca_fun(paste0(scree_method, "_control")), control)

  if (scree_method == "pmwm" && control$eta >= 0.5) {
    stop("eta must be less than 0.5.")
  }

  X_proc <- .dppca_fun("prep_matrix_for_pca")(X, center, standardize)
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
  .dppca_fun("validate_bins")(bins)

  if (is.null(colnames(X_proc))) {
    colnames(X_proc) <- paste0("V", seq_len(p))
  }
  if (anyDuplicated(colnames(X_proc))) {
    stop("Variable names must be unique.")
  }
  estimator <- .dppca_fun(paste0("dp_scree_", scree_method))

  # Estimate one full loading matrix and reuse its projected scores.
  loadings <- .dppca_fun("dp_loading")(
    X_proc,
    eps = budgets$loading[["eps"]],
    delta = budgets$loading[["delta"]],
    center = FALSE,
    standardize = FALSE,
    cpp.option = cpp.option
  )
  V <- loadings$private
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
    eigenvalues <- .dppca_fun("scree_post_processing")(eigenvalues)
  }
  names(eigenvalues) <- colnames(V)[seq_len(k)]

  # Loading has its own allocation; the score budget covers frame and histogram.
  score_budget <- .dppca_fun("split_score_privacy_budget")(privacy)
  frame <- .dppca_fun("dp_frame")(
    Y[, axes, drop = FALSE],
    score_budget$eps_frame,
    inflate = 0.2
  )
  histograms <- .dppca_fun("score_histograms")(
    Y[, axes, drop = FALSE],
    frame$xlim,
    frame$ylim,
    bins,
    score_budget$eps_hist,
    score_budget$delta_hist,
    score_method
  )

  nonprivate <- NULL
  if (non_private) {
    ordinary_eigenvalues <- eigen(
      stats::cov(X_proc),
      symmetric = TRUE,
      only.values = TRUE
    )$values
    ordinary_directions <- loadings$nonprivate
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
    ordinary_scores <- X_proc %*% ordinary_directions[, axes, drop = FALSE]
    ordinary_histogram <- .dppca_reference_histogram(ordinary_scores, frame, bins)

    nonprivate <- list(
      directions = ordinary_directions,
      eigenvalues = ordinary_eigenvalues,
      pve = .dppca_pve(ordinary_eigenvalues),
      scores = ordinary_scores,
      score_histogram = ordinary_histogram,
      frame = frame,
      variable_sd = apply(X_proc, 2, stats::sd)
    )
  }

  structure(
    list(
      directions = V,
      eigenvalues = eigenvalues,
      pve = .dppca_pve(eigenvalues),
      score_histogram = histograms[[score_method]],
      frame = frame,
      privacy = privacy,
      score_budget = score_budget,
      methods = list(scree = scree_method, score = score_method),
      controls = list(
        scree = list(k = k, method = scree_method, control = control, mono = mono),
        score = list(method = score_method, bins = bins)
      ),
      preprocessing = list(center = center, standardize = standardize),
      n = n,
      k = k,
      axes = as.integer(axes),
      nonprivate = nonprivate
    ),
    class = "dp_pca"
  )
}

# Empirical reference on the private grid: right-closed bins, with scores
# outside the frame assigned to the corresponding boundary bins.
.dppca_reference_histogram <- function(
    scores,
    frame,
    bins
) {
  grid <- .dppca_fun("score_histogram_grid")(
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

#' Plot an integrated PCA fit
#'
#' Uses the stored estimates without estimating loading, scree or score again.
#' Use `set.seed()` before plotting to reproduce synthetic samples.
#'
#' @param x A [dp_pca()] result.
#' @param type `"all"` or an ordered vector containing `"scree"`, `"score"`,
#'   `"loading"`, `"variable"` and/or `"biplot"`.
#' @param non_private Show the stored non-private comparison before DP plots.
#'   Non-private signs are aligned to the DP directions in a plotting copy.
#'   Histogram grids, probability scales and biplot display frames are shared.
#' @param scree List with `display = "pve"` (or `"eigenvalue"`),
#'   `display_k = NULL` and `plot_control = NULL`. Use [scree_plot_control()]
#'   for line, point, title and axis styles. Each panel contains one curve,
#'   so its legend options are not used. PVE uses all estimated eigenvalues.
#' @param score List with `display = "histogram"` (or `"sample"`) and
#'   `plot_control = NULL`. Appearance uses [score_plot_control()].
#' @param loading List with `display = "point"`, `heatmap_components = NULL`
#'   and `plot_control = NULL`. Display also accepts `"vector"` or `"heatmap"`.
#'   Multiple displays are available with `type = "loading"` only.
#'   Appearance uses [loading_plot_control()], including `heatmap_colors`.
#' @param variable List with `display = "vector"` (or `"point"`),
#'   `circle = TRUE` and `plot_control = NULL`. DP correlations are
#'   reconstructed using all estimated eigenvalues.
#' @param biplot List with `alpha = 0` and `plot_control = NULL`.
#'   Alpha is between zero and one. The display is always square: the shorter
#'   axis range is padded while preserving equal physical units on both axes.
#'   Tick labels show the actual alpha-transformed score coordinates.
#' @param sampling A [sampling_control()] result with one method, shared by
#'   score and biplot. By default DP panels sample `n` synthetic points and
#'   non-private panels use the original scores. An explicit `sample_size`
#'   samples both versions from their stored histograms. Histogram panels
#'   themselves do not change.
#' @param variables Common variable selection: top count, `"all"`, `NULL`,
#'   names or indices. Common top-N selection uses raw loading on the fit axes
#'   within each version. The loading, variable and biplot lists can override
#'   `variables` and `labels`; local top-N selection uses that panel's geometry.
#' @param labels Show variable labels.
#' @param base_size Common font size. Local plot controls take precedence.
#' @param ... Unrecognized arguments are rejected.
#'
#' @details
#' Biplot scores are divided by `d^alpha` and loading columns multiplied by
#' `d^alpha`, where `d = sqrt((n - 1) * eigenvalue)` on the selected axes.
#' Each version uses its own eigenvalues. A zero eigenvalue with positive
#' alpha suppresses that axis and is reported in the plot subtitle.
#'
#' A single arrow multiplier is computed from the DP panel: the 90th percentile
#' of distances from the transformed synthetic-point mean, divided by the
#' largest norm of all transformed DP variable vectors. There is no additional
#' 0.8 factor. Points keep their transformed coordinates and arrows start at
#' the origin. The same multiplier is used for the non-private comparison.
#' Empty or zero-radius DP samples give a zero multiplier and omit arrows.
#'
#' The biplot display first extends the alpha-transformed DP frame to include
#' the origin and all DP arrow endpoints, then pads the shorter axis to make
#' a square. It does not change the estimated frame or histogram.
#' Non-private plots reuse this display frame;
#' endpoints outside it are reported in the subtitle. Changing `variables`
#' alone does not change the DP multiplier or display frame.
#'
#' Variable and biplot `plot_control` accept plain lists with
#' `nonprivate_color`, `private_color`, `arrow_size`, `point_size`,
#' `point_alpha`, `label_size`, `base_size`, `title_size`, `title`, `xlab`
#' and `ylab`. `point_color` sets the synthetic/observation-point color in
#' biplots. Variable points and arrows use the corresponding version color.
#' A NULL title or axis label keeps the automatic text. Biplot defaults are
#' `point_color = "#B0B0B0"` and both version colors `"#174A7E"`.
#'
#' @return Invisibly returns `private` and, when requested, `nonprivate` lists.
#'   Each contains `plots`, `page`, `sampled_points`, `biplot`, `frames` and
#'   `direction_sign`. Biplot metadata includes the alpha-transformed and
#'   displayed coordinates, `scale`, `arrow_scale`, `score_center`,
#'   `score_radius`, and display frame. The plots are always drawn;
#'   `out <- plot(fit)` also keeps these results for reuse. Multiple loading
#'   displays with both versions additionally return a `comparison` patchwork.
#' @export
plot.dp_pca <- function(
    x,
    type = "all",
    non_private = !is.null(x$nonprivate),
    scree = list(),
    score = list(),
    loading = list(),
    variable = list(),
    biplot = list(),
    sampling = NULL,
    variables = 10,
    labels = TRUE,
    base_size = 12,
    ...
) {
  # Resolve plot selections and local options.
  valid_types <- c("scree", "score", "loading", "variable", "biplot")

  if (
    !is.character(type) || !length(type) || anyNA(type) ||
      any(!type %in% c("all", valid_types))
  ) {
    stop("Unknown plot type.", call. = FALSE)
  }

  if ("all" %in% type && length(type) > 1L) {
    stop("Use 'all' alone, or a vector of plot names.", call. = FALSE)
  }

  if (length(list(...))) {
    stop(
      "Unknown plotting argument(s): ",
      paste(names(list(...)), collapse = ", "),
      call. = FALSE
    )
  }

  types <- if (identical(type, "all")) valid_types else unique(type)

  sp <- .dppca_options(
    scree,
    list(display = "pve", display_k = NULL, plot_control = NULL),
    "scree"
  )
  hp <- .dppca_options(
    score,
    list(display = "histogram", plot_control = NULL),
    "score"
  )
  lp <- .dppca_options(
    loading,
    list(
      display = "point",
      heatmap_components = NULL,
      plot_control = NULL,
      variables = variables,
      labels = labels
    ),
    "loading"
  )
  vp <- .dppca_options(
    variable,
    list(
      display = "vector",
      circle = TRUE,
      plot_control = NULL,
      variables = variables,
      labels = labels
    ),
    "variable"
  )
  bp <- .dppca_options(
    biplot,
    list(
      alpha = 0,
      plot_control = NULL,
      variables = variables,
      labels = labels
    ),
    "biplot"
  )

  sp$display <- match.arg(sp$display, c("pve", "eigenvalue"))
  hp$display <- match.arg(hp$display, c("histogram", "sample"))
  vp$display <- match.arg(vp$display, c("vector", "point"))
  lp$display <- .dppca_fun(".validate_loading_display")(lp$display)

  if (length(lp$display) > 1L && !identical(types, "loading")) {
    stop(
      "Multiple loading displays require type = 'loading'.",
      call. = FALSE
    )
  }

  .dppca_flag(non_private, "non_private")
  .dppca_flag(labels, "labels")
  .dppca_flag(vp$circle, "variable$circle")
  .dppca_number(bp$alpha, "biplot$alpha", 0, 1)
  .dppca_number(base_size, "base_size", 0, Inf, open = TRUE)

  for (options in list(lp, vp, bp)) {
    .dppca_flag(options$labels, "labels")
  }

  scree_n <- if (is.null(sp$display_k)) x$k else sp$display_k
  .dppca_number(scree_n, "scree$display_k", 1, x$k, integer = TRUE)

  if (non_private && is.null(x$nonprivate)) {
    stop(
      "No non-private reference was stored. Refit with non_private = TRUE.",
      call. = FALSE
    )
  }

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Install ggplot2 to plot.", call. = FALSE)
  }

  # Resolve appearance and one shared sampling specification.
  pc <- .dppca_fun(".merge_scree_plot_control")(
    sp$plot_control,
    if (sp$display == "pve") "pve" else "scree"
  )
  hc <- .dppca_fun(".merge_score_plot_control")(hp$plot_control)
  lc <- .dppca_fun(".merge_loading_plot_control")(lp$plot_control, x$axes)
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

  heatpcs <- .dppca_fun(".validate_loading_components")(
    lp$heatmap_components,
    ncol(x$directions),
    x$axes,
    "heatmap_components"
  )

  if (is.null(sampling)) {
    sampling <- .dppca_fun("sampling_control")(sample_size = NULL)
  }

  # NULL means original non-private points, not a synthetic sample.
  sample_nonprivate <- is.list(sampling) && !is.null(sampling$sample_size)
  sampling <- .dppca_fun(".resolve_sampling_control")(sampling, x$n)

  if (length(sampling$method) != 1L) {
    stop("Choose one shared sampling method.", call. = FALSE)
  }

  need_points <- "biplot" %in% types ||
    ("score" %in% types && hp$display == "sample")
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
      obj$scores <- sweep(obj$scores, 2, direction_sign[x$axes], "*")
      obj$score_histogram <- .dppca_reference_histogram(
        obj$scores,
        x$frame,
        x$controls$score$bins
      )
      obj$frame <- x$frame
    }

    ax <- x$axes
    V <- obj$directions
    ev <- obj$eigenvalues
    prefix <- if (private) "DP" else "Non-private"
    page_title <- if (private) "Private PCA" else "Non-private PCA"
    score_frame <- x$frame
    plot_frames <- list(score = score_frame)
    points <- data.frame(x = numeric(), y = numeric())

    if (need_points) {
      if (private || sample_nonprivate) {
        sampled <- .dppca_fun("sample_private_score_histogram")(
          obj$score_histogram,
          sampling$sample_size,
          sampling$method,
          sampling$bandwidth_scale
        )[[1]]
        points <- data.frame(x = sampled$pc_x, y = sampled$pc_y)
      } else {
        points <- data.frame(x = obj$scores[, 1], y = obj$scores[, 2])
      }
    }

    common_idx <- .dppca_fun(".select_loading_variables")(
      variables,
      rownames(V),
      V,
      ax
    )
    common_names <- rownames(V)[common_idx]
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

      local_options <- switch(tp, loading = lp, variable = vp, biplot = bp)
      local_supplied <- switch(
        tp,
        loading = loading,
        variable = variable,
        biplot = biplot,
        list()
      )
      local_variables <- if ("variables" %in% names(local_supplied)) {
        local_options$variables
      } else {
        common_names
      }
      local_labels <- if (is.null(local_options)) labels else local_options$labels

      if (tp == "scree") {
        values <- if (sp$display == "pve") .dppca_pve(ev) else ev
        idx <- seq_len(scree_n)
        dat <- data.frame(pc = idx, value = values[idx])
        series <- if (private) x$methods$scree else "nonprivate"

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
            y = if (sp$display == "pve") {
              "Proportion of Variance Explained"
            } else {
              "Eigenvalue"
            }
          )

        if (sp$display == "pve" && scree_n < x$k) {
          pl <- pl + ggplot2::labs(
            caption = paste("PVE denominator: all", x$k, "estimated PCs")
          )
        }
        if (sum(ev) <= 0) {
          note <- "All estimated eigenvalues are zero"
        }
      } else if (tp == "score") {
        if (hp$display == "histogram") {
          pl <- .dppca_fun("make_hist_plot_dp")(
            obj$score_histogram,
            score_frame$xlim,
            score_frame$ylim,
            hc$color,
            alpha_range = hc$hist_alpha_range,
            base_size = hc$base_size,
            title_size = hc$title_size
          )

          # Both panels use the DP probability-to-opacity scale.
          alpha_max <- max(x$score_histogram$prob, na.rm = TRUE)
          if (!is.finite(alpha_max) || alpha_max <= 0) {
            alpha_max <- 1
          }
          alpha_scale <- pl$scales$get_scales("alpha")
          alpha_scale$limits <- c(0, alpha_max)
          alpha_scale$oob <- function(z, range, ...) {
            pmin(pmax(z, range[1]), range[2])
          }

          if (sum(obj$score_histogram$prob) <= 0) {
            note <- "No bins survived private thresholding"
          }
        } else {
          pl <- .dppca_fun("make_sample_plot_dp")(
            data.frame(pc_x = points$x, pc_y = points$y),
            score_frame$xlim,
            score_frame$ylim,
            hc$color,
            point_alpha = hc$scatter_alpha,
            point_size = hc$scatter_size,
            base_size = hc$base_size,
            title_size = hc$title_size
          )

          if (!nrow(points)) {
            note <- "No synthetic points: histogram has zero mass"
          }
        }
      } else if (tp == "loading") {
        for (display in lp$display) {
          components <- if (display == "heatmap") heatpcs else ax
          selected_variables <- if (
            display == "heatmap" && !"variables" %in% names(loading)
          ) {
            variables
          } else {
            local_variables
          }
          selected <- .dppca_fun(".select_loading_variables")(
            selected_variables,
            rownames(V),
            V,
            components
          )
          loading_title <- paste(
            prefix,
            if (display == "heatmap") "Loading Heatmap" else "Loading Plot"
          )
          title_key <- if (private) "private_title" else "nonprivate_title"

          if (
            title_key %in% names(lp$plot_control) &&
              (is.null(lc[[title_key]]) ||
                !lc[[title_key]] %in% c("DP", "Non-private", "Non-Private"))
          ) {
            loading_title <- lc[[title_key]]
          }

          if (display == "heatmap") {
            dat <- .dppca_fun(".make_loading_heatmap_data")(
              V, components, selected, rownames(V)
            )
            pl <- .dppca_fun(".build_loading_heatmap")(
              dat, loading_title, 1, lc
            )
          } else {
            dat <- .dppca_fun(".make_loading_data")(
              V, ax, selected, rownames(V)
            )
            limit <- max(0.1, abs(c(dat$x, dat$y))) * 1.2
            pl <- .dppca_fun(".build_loading_panel")(
              dat,
              if (private) lc$private_color else lc$nonprivate_color,
              loading_title,
              lc$xlab,
              lc$ylab,
              local_labels,
              limit,
              display,
              lc
            )
          }

          key <- if (length(lp$display) == 1L) {
            "loading"
          } else {
            paste0("loading_", display)
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
          coords <- .dppca_biplot_coordinates(
            points = points,
            V = V,
            eigenvalues = ev,
            n = x$n,
            axes = ax,
            alpha = bp$alpha,
            frame = score_frame,
            arrow_scale = if (private) NULL else result$private$biplot$arrow_scale
          )
          A <- coords$display_variables
          note <- coords$note
          pl <- ggplot2::ggplot(
            coords$display_observations,
            ggplot2::aes(x = .data$x, y = .data$y)
          ) +
            ggplot2::geom_point(
              color = bc$point_color,
              alpha = bc$point_alpha,
              size = bc$point_size
            )
        }

        keep <- .dppca_variables(A, rownames(V), local_variables)
        if (tp == "biplot" && coords$arrow_scale == 0) {
          keep <- integer()
        }
        dat <- data.frame(
          x = A[keep, 1],
          y = A[keep, 2],
          label = rownames(V)[keep]
        )

        pl <- pl +
          ggplot2::geom_hline(
            yintercept = 0, color = "grey75", linewidth = 0.4,
            linetype = "dashed"
          ) +
          ggplot2::geom_vline(
            xintercept = 0, color = "grey75", linewidth = 0.4,
            linetype = "dashed"
          )

        if (tp == "variable" && vp$display == "point") {
          pl <- pl + ggplot2::geom_point(
            data = dat,
            ggplot2::aes(x = .data$x, y = .data$y),
            inherit.aes = FALSE,
            color = color,
            alpha = style$point_alpha,
            size = style$point_size
          )
        } else {
          pl <- pl + ggplot2::geom_segment(
            data = dat,
            ggplot2::aes(x = 0, y = 0, xend = .data$x, yend = .data$y),
            inherit.aes = FALSE,
            color = color,
            linewidth = style$arrow_size,
            arrow = grid::arrow(
              length = grid::unit(0.18, "cm"), type = "closed"
            )
          )
        }

        if (local_labels && nrow(dat)) {
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

        if (tp == "variable" && vp$circle) {
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
        } else if (tp == "variable") {
          limit <- max(0.1, abs(A[is.finite(A)])) * 1.2
          pl <- pl + ggplot2::coord_equal(
            xlim = c(-limit, limit), ylim = c(-limit, limit)
          )
        } else {
          biplot_frame <- if (private) {
            .dppca_biplot_frame(coords)
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
      score_title <- if (hp$display == "histogram") {
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
        pl <- pl + .dppca_fun("theme_dp_base")(panel_base)
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

      override <- switch(
        tp,
        scree = sp$plot_control,
        score = hp$plot_control,
        local_options$plot_control
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
        title_key <- if (!private && hp$display == "sample") {
          "scatter_title"
        } else if (!private) {
          "nonprivate_title"
        } else {
          paste0(x$methods$score, "_title")
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
          plot.title = ggplot2::element_text(size = base_size * pc$cex_main),
          axis.title = ggplot2::element_text(size = base_size * pc$cex_lab),
          axis.text = ggplot2::element_text(size = base_size * pc$cex_axis)
        )
      }
      plots[[tp]] <- pl
    }

    heatmap_name <- if (
      "loading" %in% types && identical(lp$display, "heatmap")
    ) {
      "loading"
    } else {
      NULL
    }
    page <- .dppca_page(
      plots, page_title, heatmap_name = heatmap_name, base_size = base_size
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

  if (identical(types, "loading") && length(lp$display) > 1L && non_private) {
    if (!requireNamespace("patchwork", quietly = TRUE)) {
      stop("Install patchwork for loading comparison.", call. = FALSE)
    }
    result$comparison <- patchwork::wrap_plots(
      c(unname(result$nonprivate$plots), unname(result$private$plots)),
      ncol = length(lp$display),
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

# Include the origin and DP arrows, then pad to an equal-unit square display.
.dppca_biplot_frame <- function(
    coordinates
) {
  frame <- coordinates$transformed_score_frame
  arrows <- coordinates$display_variables
  frame$xlim <- range(frame$xlim, 0, arrows[, 1])
  frame$ylim <- range(frame$ylim, 0, arrows[, 2])

  .dppca_square_frame(frame)
}

#' @export
print.dp_pca <- function(
    x,
    digits = max(3L, getOption("digits") - 3L),
    ...
) {
  .dppca_number(digits, "digits", 1, 22, integer = TRUE)

  cat("Differentially Private PCA\n")
  cat("Data:", x$n, "observations x", nrow(x$directions), "variables\n")
  cat(
    "Components:", ncol(x$directions), "loadings;",
    x$k, "scree estimates\n"
  )
  cat("Score axes:", paste0("PC", x$axes, collapse = ", "), "\n")
  cat(
    "Methods: scree =", x$methods$scree,
    "; score =", x$methods$score, "\n"
  )
  cat(
    "Preprocessing: center =", x$preprocessing$center,
    "; standardize =", x$preprocessing$standardize, "\n"
  )
  cat(
    "Allocated budget: eps =", format(x$privacy$total[["eps"]], digits = digits),
    "; delta =", format(x$privacy$total[["delta"]], digits = digits), "\n"
  )
  cat(
    "Non-private reference:",
    if (is.null(x$nonprivate)) "not stored" else "included", "\n"
  )

  displayed <- seq_len(min(6L, length(x$eigenvalues)))
  cat("\nDP eigenvalues")
  if (length(displayed) < length(x$eigenvalues)) {
    cat(
      " (first ", length(displayed), " of ", length(x$eigenvalues), ")",
      sep = ""
    )
  }
  cat(":\n")
  print(x$eigenvalues[displayed], digits = digits, ...)
  cat("\nUse summary(x) for component and budget tables; plot(x) for plots.\n")

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
#'   data dimensions, selected axes, methods, preprocessing and whether a
#'   non-private reference is stored. Raw scores and loading matrices are
#'   not copied into the summary.
#'
#' @details
#' PVE is normalized over the `k` estimated eigenvalues. If `k` is smaller
#' than the number of variables, cumulative PVE of one does not mean that
#' the first `k` components explain all of the data's variance. If all
#' estimated eigenvalues are zero, PVE and cumulative PVE are shown as zero.
#'
#' The budget table reports allocations, not measured privacy loss. In
#' particular, the pure-DP frame uses no delta; a scree method may leave
#' some of its allocated delta unused. Non-private comparison values are
#' not included in the importance table.
#' @export
summary.dp_pca <- function(
    object,
    ...
) {
  if (length(list(...))) {
    stop("summary.dp_pca() takes no additional options.", call. = FALSE)
  }

  eigenvalues <- object$eigenvalues
  pve <- .dppca_pve(eigenvalues)
  importance <- rbind(
    Eigenvalue = eigenvalues,
    PVE = pve,
    `Cumulative PVE` = cumsum(pve)
  )
  colnames(importance) <- names(eigenvalues)

  budgets <- object$privacy$components
  score_budget <- object$score_budget
  budget <- rbind(
    Loading = budgets$loading,
    Scree = budgets$scree,
    `Score frame` = c(eps = score_budget$eps_frame, delta = 0),
    `Score histogram` = c(
      eps = score_budget$eps_hist,
      delta = score_budget$delta_hist
    ),
    Total = object$privacy$total
  )

  structure(
    list(
      n = object$n,
      p = nrow(object$directions),
      k = object$k,
      axes = object$axes,
      methods = object$methods,
      preprocessing = object$preprocessing,
      importance = importance,
      budget = budget,
      nonprivate = !is.null(object$nonprivate)
    ),
    class = "summary.dp_pca"
  )
}

#' @param x A `summary.dp_pca` object.
#' @param digits Significant digits used when printing tables.
#' @rdname summary.dp_pca
#' @export
print.summary.dp_pca <- function(
    x,
    digits = max(3L, getOption("digits") - 3L),
    ...
) {
  .dppca_number(digits, "digits", 1, 22, integer = TRUE)

  cat("DP PCA summary\n")
  cat("Data:", x$n, "observations x", x$p, "variables\n")
  cat(
    "Methods: scree =", x$methods$scree,
    "; score =", x$methods$score, "\n"
  )
  cat("Score axes:", paste0("PC", x$axes, collapse = ", "), "\n")

  cat("\nDP component importance:\n")
  print(x$importance, digits = digits, ...)
  cat("\nProportions are normalized over", x$k, "estimated components.\n")
  if (x$k < x$p) {
    cat("Unestimated components are not included in the denominator.\n")
  }
  if (sum(x$importance["Eigenvalue", ]) <= 0) {
    cat("All estimated eigenvalues are zero; proportions are shown as zero.\n")
  }

  cat("\nAllocated privacy budget:\n")
  print(x$budget, digits = digits, ...)
  cat(
    "\nNon-private reference:",
    if (x$nonprivate) "included in fit" else "not stored", "\n"
  )

  invisible(x)
}

# Helpers -------------------------------------------------------------------

# Lookup supports both devtools::load_all() and source() with installed dppca.
.dppca_fun <- function(
    name
) {
  fun <- get0(
    name,
    envir = environment(.dppca_fun),
    mode = "function",
    inherits = TRUE
  )

  if (is.null(fun) && requireNamespace("dppca", quietly = TRUE)) {
    fun <- get0(
      name,
      envir = asNamespace("dppca"),
      mode = "function",
      inherits = FALSE
    )
  }

  if (is.null(fun)) {
    stop(
      "Missing existing dppca function: ", name,
      ". Load the supplied package R files with devtools::load_all() first."
    )
  }

  fun
}

.dppca_flag <- function(
    x,
    name
) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop(name, " must be TRUE or FALSE.")
  }
}

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

.dppca_pve <- function(
    x
) {
  if (sum(x) > 0) {
    x / sum(x)
  } else {
    stats::setNames(rep(0, length(x)), names(x))
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

  selected <- .dppca_fun(".select_loading_variables")(
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
    base_size = 12
) {
  n <- length(plots)
  tall_heatmap <- n == 5L && !is.null(heatmap_name) &&
    heatmap_name %in% names(plots)

  if (tall_heatmap) {
    plots <- c(
      plots[setdiff(names(plots), heatmap_name)],
      plots[heatmap_name]
    )
  }

  nr <- if (n <= 3L) {
    1L
  } else {
    2L
  }
  nc <- if (tall_heatmap) {
    3L
  } else if (n == 4L) {
    2L
  } else {
    min(n, 3L)
  }
  layout_cols <- if (n == 5L && !tall_heatmap) {
    6L
  } else {
    nc
  }
  column_widths <- if (tall_heatmap) {
    c(1, 1, 1.15)
  } else {
    rep(1, layout_cols)
  }

  layout <- grid::grid.layout(
    nr + 1L,
    layout_cols,
    widths = grid::unit(column_widths, "null"),
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
        layout.pos.row = 1,
        layout.pos.col = seq_len(layout_cols)
      )
    )
  )
  grobs <- lapply(plots, ggplot2::ggplotGrob)

  # Align only ordinary panels. A heatmap has its own row labels and legend.
  align <- if (tall_heatmap) {
    seq_len(4L)
  } else if (!is.null(heatmap_name)) {
    which(names(plots) != heatmap_name)
  } else {
    seq_along(grobs)
  }

  if (length(align) > 1L) {
    width_counts <- vapply(
      grobs[align],
      function(g) {
        length(g$widths)
      },
      integer(1)
    )

    if (length(unique(width_counts)) == 1L) {
      widths <- do.call(
        grid::unit.pmax,
        lapply(grobs[align], function(g) {
          g$widths
        })
      )

      for (i in align) {
        grobs[[i]]$widths <- widths
      }
    }
  }

  for (i in seq_along(plots)) {
    if (tall_heatmap) {
      row <- if (i == 5L) {
        2:3
      } else {
        2L + (i - 1L) %/% 2L
      }
      col <- if (i == 5L) {
        3L
      } else {
        1L + (i - 1L) %% 2L
      }
    } else if (n == 5L) {
      row <- if (i <= 3L) {
        2L
      } else {
        3L
      }
      col <- switch(i, 1:2, 3:4, 5:6, 2:3, 4:5)
    } else {
      row <- 2L + (i - 1L) %/% nc
      col <- 1L + (i - 1L) %% nc
    }

    children[[i + 1L]] <- grid::grobTree(
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
    .dppca_fun(".validate_plot_color")(out[[nm]], nm)
  }

  for (nm in c("arrow_size", "point_size", "label_size", "title_size")) {
    .dppca_number(out[[nm]], nm, 0, Inf, open = TRUE)
  }

  .dppca_number(out$point_alpha, "point_alpha", 0, 1)

  for (nm in c("title", "xlab", "ylab")) {
    .dppca_fun(".validate_plot_optional_string")(out[[nm]], nm)
  }

  out
}
