# ============================================================
# loading.R
# Functions for differentially private loading plot
# ============================================================

#' Estimate non-private and differentially private loadings
#'
#' @param X Numeric matrix or data frame; rows are observations.
#' @param eps Positive privacy parameter epsilon.
#' @param delta Privacy parameter delta in `(0, 1)`.
#' @param center Whether to center columns. Default: `TRUE`.
#' @param standardize Whether to divide columns by their sample standard
#'   deviations. Default: `FALSE`.
#' @param cpp.option Use the compiled spherical Kendall implementation from
#'   dppca. Default: `TRUE`. Use `FALSE` for the R implementation in this file.
#'
#' @return A list with two numeric `p` by `p` matrices:
#'   \itemize{
#'     \item `nonprivate`: ordinary PCA directions from the sample covariance.
#'     \item `private`: DP directions from the noisy spherical Kendall matrix.
#'   }
#'   Rows are variables; columns are `PC1`, ..., `PCp`. Both matrices contain
#'   orthonormal directions, without eigenvalue or correlation scaling.
#' @details Both estimates use the same preprocessed data. All private
#'   directions come from one noisy spherical Kendall matrix; there is no `k`
#'   or direction-estimation toggle. The full result includes a non-private
#'   reference. Sample standardization uses data-dependent scales and requires
#'   separate privacy analysis.
#' @examples
#' set.seed(123)
#' X <- matrix(rnorm(150), 50, 3)
#' result <- dp_loading(X, eps = 2, delta = 1e-3, cpp.option = FALSE)
#' result$nonprivate
#' result$private
#' @export
dp_loading <- function(
    X,
    eps,
    delta,
    center = TRUE,
    standardize = FALSE,
    cpp.option = TRUE
) {
  # Prepare the data and validate privacy parameters.
  X <- .loading_matrix(X, center, standardize)

  .loading_positive(eps, "eps")
  .loading_positive(delta, "delta")

  if (delta >= 1) {
    stop("`delta` must be less than 1.", call. = FALSE)
  }

  .loading_flag(cpp.option, "cpp.option")

  n <- nrow(X)
  p <- ncol(X)

  # Add symmetric Gaussian noise to the spherical Kendall matrix.
  kendall <- .loading_kendall(X, cpp.option)
  sigma <- 4 * sqrt(2 * log(1.25 / delta)) / (n * eps)
  z <- stats::rnorm(p * (p + 1) / 2, sd = sigma)

  noise <- matrix(0, p, p)
  diag(noise) <- z[seq_len(p)]
  noise[lower.tri(noise)] <- z[-seq_len(p)] / sqrt(2)
  noise[upper.tri(noise)] <- t(noise)[upper.tri(noise)]

  # Estimate all directions for both versions.
  V_private <- eigen(kendall + noise, symmetric = TRUE)$vectors
  V_nonprivate <- eigen(stats::cov(X), symmetric = TRUE)$vectors

  loading_names <- list(colnames(X), paste0("PC", seq_len(p)))
  dimnames(V_nonprivate) <- loading_names
  dimnames(V_private) <- loading_names

  list(
    nonprivate = V_nonprivate,
    private = V_private
  )
}

#' Build loading plots directly from data
#'
#' @inheritParams dp_loading
#' @param axes Two distinct PCs for vector and point displays.
#' @param display One or more of `"vector"`, `"point"`, `"heatmap"`.
#' @param heatmap_components PCs in the heatmap; `NULL` uses `axes`.
#' @param variables A positive number selects top variables separately in each
#'   panel. Use `"all"` or `NULL`, variable names, or multiple column indices
#'   for other selections. Default: top 10.
#' @param labels Show variable labels in vector and point displays.
#' @param align_sign Align only non-private plotting directions to the private
#'   directions. Default: `TRUE`. Returned loading matrices are unchanged.
#' @param plot_control Optional output of [loading_plot_control()].
#' @return A list with `loading` (the result of [dp_loading()]) and `plot`
#'   (display-named `nonprivate` and `private` lists, plus `all`). Print
#'   `result$plot$all` to draw the combined patchwork.
#' @details One display has non-private left and private right. Multiple
#'   displays have non-private above private, ordered vector, point, heatmap.
#'   Non-private panels and comparison-based plot ranges are reference output;
#'   the combined result is not a private-only release.
#' @examples
#' \dontrun{
#' set.seed(123)
#' X <- matrix(rnorm(150), 50, 3)
#' out <- dp_loading_plot(X, eps = 2, delta = 1e-3, cpp.option = FALSE)
#' out$loading$nonprivate
#' out$loading$private
#' out$plot$all
#' }
#' @export
dp_loading_plot <- function(
    X,
    eps,
    delta,
    axes = c(1, 2),
    display = "vector",
    heatmap_components = NULL,
    center = TRUE,
    standardize = FALSE,
    cpp.option = TRUE,
    variables = 10,
    labels = TRUE,
    align_sign = TRUE,
    plot_control = NULL
) {
  # Prepare the data and check the plot settings.
  X <- .loading_matrix(X, center, standardize)

  if (ncol(X) < 2L) {
    stop("Loading plots need at least two variables.", call. = FALSE)
  }

  if (anyDuplicated(colnames(X))) {
    stop("Variable names must be unique.", call. = FALSE)
  }

  axes <- .validate_loading_axes(axes, ncol(X))
  display <- .validate_loading_display(display)
  components <- .validate_loading_components(
    components = heatmap_components,
    p = ncol(X),
    default = axes,
    argument = "heatmap_components"
  )

  .loading_flag(labels, "labels")
  .loading_flag(align_sign, "align_sign")

  control <- .merge_loading_plot_control(plot_control, axes)
  .loading_plot_packages(display, labels)

  # X is already preprocessed, so do not center or standardize it again.
  loading <- dp_loading(
    X = X,
    eps = eps,
    delta = delta,
    center = FALSE,
    standardize = FALSE,
    cpp.option = cpp.option
  )

  .plot_loading_result(
    loading = loading,
    axes = axes,
    display = display,
    heatmap_components = components,
    variables = variables,
    labels = labels,
    align_sign = align_sign,
    plot_control = control
  )
}

#' Control loading plot appearance
#'
#' @param nonprivate_color,private_color Colors for vectors and points.
#' @param arrow_size,point_size,label_size Vector width, point size, label size.
#' @param heatmap_values Print coefficients inside heatmap cells.
#' @param heatmap_digits Number of decimal places in heatmap cells.
#' @param heatmap_text_size Heatmap cell text size.
#' @param base_size,title_size Base font and title sizes.
#' @param nonprivate_title,private_title Panel titles; `NULL` hides them.
#' @param xlab,ylab Axis labels; `NULL` uses the selected PC names.
#' @param heatmap_colors Named character vector with `low` and `high` colors.
#'   Both heatmaps use these colors, with white at zero.
#' @return A named list of plot settings.
#' @export
loading_plot_control <- function(
    nonprivate_color = "black",
    private_color = "#1F4E79",
    arrow_size = 0.7,
    point_size = 2.8,
    label_size = 4,
    heatmap_values = TRUE,
    heatmap_digits = 2L,
    heatmap_text_size = 3.5,
    base_size = 12,
    title_size = 12,
    nonprivate_title = "Non-private",
    private_title = "DP",
    xlab = NULL,
    ylab = NULL,
    heatmap_colors = c(low = "#2C5D8A", high = "#B44A4A")
) {
  control <- list(
    nonprivate_color = nonprivate_color,
    private_color = private_color,
    arrow_size = arrow_size,
    point_size = point_size,
    label_size = label_size,
    heatmap_values = heatmap_values,
    heatmap_digits = heatmap_digits,
    heatmap_text_size = heatmap_text_size,
    heatmap_colors = heatmap_colors,
    base_size = base_size,
    title_size = title_size,
    nonprivate_title = nonprivate_title,
    private_title = private_title,
    xlab = xlab,
    ylab = ylab
  )

  size_options <- c(
    "arrow_size",
    "point_size",
    "label_size",
    "heatmap_text_size",
    "base_size",
    "title_size"
  )

  for (name in size_options) {
    .loading_positive(control[[name]], name)
  }

  .loading_flag(heatmap_values, "heatmap_values")

  if (
    !is.numeric(heatmap_digits) ||
      length(heatmap_digits) != 1L ||
      !is.finite(heatmap_digits) ||
      heatmap_digits < 0 ||
      heatmap_digits != floor(heatmap_digits)
  ) {
    stop("`heatmap_digits` must be a non-negative integer.", call. = FALSE)
  }

  for (name in c("nonprivate_color", "private_color")) {
    value <- control[[name]]

    if (!is.character(value) || length(value) != 1L || is.na(value)) {
      stop("`", name, "` must be a color string.", call. = FALSE)
    }

    tryCatch(
      grDevices::col2rgb(value),
      error = function(e) {
        stop("Invalid color for `", name, "`.", call. = FALSE)
      }
    )
  }

  if (
    !is.character(heatmap_colors) ||
      length(heatmap_colors) != 2L ||
      anyNA(heatmap_colors) ||
      !setequal(names(heatmap_colors), c("low", "high"))
  ) {
    stop(
      "`heatmap_colors` must be a character vector named 'low' and 'high'.",
      call. = FALSE
    )
  }

  tryCatch(
    grDevices::col2rgb(heatmap_colors),
    error = function(e) {
      stop("Invalid color in `heatmap_colors`.", call. = FALSE)
    }
  )

  for (name in c("nonprivate_title", "private_title", "xlab", "ylab")) {
    value <- control[[name]]

    if (
      !is.null(value) &&
        (!is.character(value) || length(value) != 1L || is.na(value))
    ) {
      stop("`", name, "` must be NULL or one string.", call. = FALSE)
    }
  }

  control
}

# Input preparation and estimation -------------------------------------------

#' Preprocess a matrix for principal component analysis
#'
#' @param X A numeric matrix or data frame.
#' @param center A logical value indicating whether to center columns.
#' @param standardize A logical value indicating whether to scale columns.
#'
#' @return A numeric matrix after preprocessing.
#' @noRd
prep_matrix_for_pca <- function(
    X,
    center = TRUE,
    standardize = FALSE
) {
  .loading_flag(center, "center")
  .loading_flag(standardize, "standardize")

  if (is.data.frame(X)) {
    numeric_cols <- vapply(X, is.numeric, logical(1))
    if (!all(numeric_cols)) {
      stop("All columns of `X` must be numeric.", call. = FALSE)
    }
  }

  X <- as.matrix(X)

  if (!is.numeric(X)) {
    stop("`X` must be a numeric matrix or data frame.", call. = FALSE)
  }

  storage.mode(X) <- "double"

  if (nrow(X) < 2) {
    stop("`X` must have at least two rows.", call. = FALSE)
  }
  if (ncol(X) < 1) {
    stop("`X` must have at least one column.", call. = FALSE)
  }
  if (anyNA(X) || any(!is.finite(X))) {
    stop("`X` must contain only finite, non-missing values.", call. = FALSE)
  }

  X_proc <- X

  if (center) {
    X_proc <- scale(X_proc, center = TRUE, scale = FALSE)
  }

  if (standardize) {
    sds <- apply(X_proc, 2, stats::sd)
    if (any(!is.finite(sds)) || any(sds <= 0)) {
      stop(
        "All columns of `X` must have positive finite standard deviations ",
        "when `standardize = TRUE`.",
        call. = FALSE
      )
    }
    X_proc <- sweep(X_proc, 2, sds, "/")
  }

  X_proc
}

.loading_flag <- function(
    x,
    name
) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop("`", name, "` must be TRUE or FALSE.", call. = FALSE)
  }
}

.loading_positive <- function(
    x,
    name
) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0) {
    stop("`", name, "` must be a positive finite number.", call. = FALSE)
  }
}

.loading_matrix <- function(
    X,
    center,
    standardize
) {
  .loading_flag(center, "center")
  .loading_flag(standardize, "standardize")

  if (is.data.frame(X) && !all(vapply(X, is.numeric, logical(1)))) {
    stop("All columns of `X` must be numeric.", call. = FALSE)
  }

  X <- as.matrix(X)

  if (
    !is.numeric(X) ||
      nrow(X) < 2L ||
      ncol(X) < 1L ||
      any(!is.finite(X))
  ) {
    stop(
      "`X` must be finite numeric data with at least two rows and one column.",
      call. = FALSE
    )
  }

  storage.mode(X) <- "double"

  if (is.null(colnames(X))) {
    colnames(X) <- paste0("X", seq_len(ncol(X)))
  }

  if (center) {
    X <- scale(X, center = TRUE, scale = FALSE)
  }

  if (standardize) {
    sds <- apply(X, 2L, stats::sd)

    if (any(!is.finite(sds)) || any(sds <= 0)) {
      stop(
        "Standardization requires positive finite column standard deviations.",
        call. = FALSE
      )
    }

    X <- sweep(X, 2L, sds, "/")
  }

  X
}

.loading_kendall <- function(
    X,
    cpp.option
) {
  if (cpp.option) {
    # Find the C++ backend in this package or an installed dppca package.
    cpp <- get0("tau_sph_cpp", mode = "function", inherits = TRUE)

    if (is.null(cpp) && requireNamespace("dppca", quietly = TRUE)) {
      cpp <- get0(
        "tau_sph_cpp",
        envir = asNamespace("dppca"),
        mode = "function",
        inherits = FALSE
      )
    }

    if (is.null(cpp)) {
      stop(
        "`cpp.option = TRUE` needs dppca's compiled backend. Use FALSE for R only.",
        call. = FALSE
      )
    }

    return(cpp(X))
  }

  # R implementation: accumulate normalized pairwise differences.
  n <- nrow(X)
  kendall <- matrix(0, ncol(X), ncol(X))

  for (i in seq_len(n - 1L)) {
    differences <- sweep(
      X[(i + 1L):n, , drop = FALSE],
      MARGIN = 2L,
      STATS = X[i, ],
      FUN = "-"
    )
    lengths <- sqrt(rowSums(differences^2))

    if (
      any(!is.finite(lengths)) ||
        any(lengths <= sqrt(.Machine$double.eps))
    ) {
      stop(
        "Spherical Kendall directions require distinct, ",
        "numerically separated observations.",
        call. = FALSE
      )
    }

    unit_differences <- sweep(differences, 1L, lengths, "/")
    kendall <- kendall + crossprod(unit_differences)
  }

  kendall * (2 / (n * (n - 1)))
}

# Display selection ----------------------------------------------------------

.validate_loading_components <- function(
    components,
    p,
    default,
    argument = "components"
) {
  if (is.null(components)) {
    components <- default
  }

  if (
    !is.numeric(components) ||
      !length(components) ||
      any(!is.finite(components)) ||
      any(components != floor(components)) ||
      any(components < 1 | components > p)
  ) {
    stop(
      "`", argument, "` must contain integers between 1 and ", p, ".",
      call. = FALSE
    )
  }

  unique(as.integer(components))
}

.validate_loading_axes <- function(
    axes,
    p
) {
  if (!is.numeric(axes) || length(axes) != 2L || anyDuplicated(axes)) {
    stop("`axes` must contain two distinct PC indices.", call. = FALSE)
  }

  .validate_loading_components(axes, p, NULL, "axes")
}

.validate_loading_display <- function(
    display
) {
  valid <- c("vector", "point", "heatmap")

  if (
    !is.character(display) ||
      !length(display) ||
      anyNA(display) ||
      any(!display %in% valid)
  ) {
    stop("`display` must contain vector, point, or heatmap.", call. = FALSE)
  }

  valid[valid %in% display]
}

.select_loading_variables <- function(
    variables,
    variable_names,
    V_loading,
    components
) {
  p <- length(variable_names)

  if (is.null(variables)) {
    return(seq_len(p))
  }

  if (
    is.character(variables) &&
      length(variables) == 1L &&
      identical(tolower(variables), "all")
  ) {
    return(seq_len(p))
  }

  # A single number selects the variables with the largest loading magnitudes.
  if (is.numeric(variables) && length(variables) == 1L) {
    if (
      !is.finite(variables) ||
        variables < 1 ||
        variables != floor(variables)
    ) {
      stop("A single `variables` value must be a positive integer.", call. = FALSE)
    }

    magnitude <- sqrt(rowSums(V_loading[, components, drop = FALSE]^2))
    ranked <- order(magnitude, decreasing = TRUE)

    return(ranked[seq_len(min(variables, p))])
  }

  if (is.character(variables) && length(variables) && !anyNA(variables)) {
    unknown <- setdiff(variables, variable_names)

    if (length(unknown)) {
      stop(
        "Unknown variables: ", paste(unknown, collapse = ", "),
        call. = FALSE
      )
    }

    return(match(unique(variables), variable_names))
  }

  if (is.numeric(variables) && length(variables) > 1L) {
    return(.validate_loading_components(variables, p, NULL, "variables"))
  }

  stop(
    "Use a top-variable count, names, column indices, or 'all' for `variables`.",
    call. = FALSE
  )
}

.merge_loading_plot_control <- function(
    plot_control,
    axes
) {
  defaults <- loading_plot_control()

  if (!is.null(plot_control)) {
    option_names <- names(plot_control)

    if (
      !is.list(plot_control) ||
        (length(plot_control) && (
          is.null(option_names) ||
            anyNA(option_names) ||
            anyDuplicated(option_names) ||
            any(!option_names %in% names(defaults))
        ))
    ) {
      stop(
        "`plot_control` must be a named list of loading_plot_control options.",
        call. = FALSE
      )
    }

    defaults[option_names] <- plot_control
  }

  control <- do.call(loading_plot_control, defaults)

  if (is.null(control$xlab)) {
    control$xlab <- paste0("PC", axes[1])
  }

  if (is.null(control$ylab)) {
    control$ylab <- paste0("PC", axes[2])
  }

  control
}

.loading_plot_packages <- function(
    display,
    labels
) {
  packages <- c("ggplot2", "patchwork")

  if (labels && any(display %in% c("vector", "point"))) {
    packages <- c(packages, "ggrepel")
  }

  for (package in packages) {
    if (!requireNamespace(package, quietly = TRUE)) {
      stop("Install '", package, "' to create loading plots.", call. = FALSE)
    }
  }
}

# Build plots from estimated matrices ----------------------------------------

.plot_loading_result <- function(
    loading,
    axes,
    display,
    heatmap_components,
    variables,
    labels,
    align_sign,
    plot_control
) {
  V_private <- loading$private
  V_nonprivate <- loading$nonprivate
  variable_names <- rownames(V_private)
  control <- .merge_loading_plot_control(plot_control, axes)

  # Align a plotting copy only; keep the returned loading matrices unchanged.
  if (align_sign) {
    signs <- ifelse(colSums(V_nonprivate * V_private) < 0, -1, 1)
    V_nonprivate <- sweep(V_nonprivate, 2L, signs, "*")
  }

  matrices <- list(
    nonprivate = V_nonprivate,
    private = V_private
  )
  plots <- list(
    nonprivate = list(),
    private = list()
  )
  selected <- list()

  for (kind in display) {
    components <- if (kind == "heatmap") {
      heatmap_components
    } else {
      axes
    }

    # Select variables separately for each version.
    indices <- lapply(matrices, function(V) {
      .select_loading_variables(
        variables = variables,
        variable_names = variable_names,
        V_loading = V,
        components = components
      )
    })

    selected[[kind]] <- lapply(indices, function(i) variable_names[i])

    data <- lapply(names(matrices), function(version) {
      V <- matrices[[version]]

      if (kind == "heatmap") {
        .make_loading_heatmap_data(
          V = V,
          components = components,
          selected = indices[[version]],
          variable_names = variable_names
        )
      } else {
        .make_loading_data(
          V = V,
          axes = axes,
          selected = indices[[version]],
          variable_names = variable_names
        )
      }
    })
    names(data) <- names(matrices)

    # Use the same scale for the non-private and private panels.
    values <- unlist(lapply(data, function(d) {
      if (kind == "heatmap") {
        d$loading
      } else {
        c(d$x, d$y)
      }
    }))
    limit <- max(abs(values))

    if (!is.finite(limit) || limit <= sqrt(.Machine$double.eps)) {
      limit <- 1
    }

    for (version in names(matrices)) {
      base_title <- control[[paste0(version, "_title")]]
      suffix <- switch(
        kind,
        vector = "Vector",
        point = "Point",
        heatmap = "Heatmap"
      )
      title <- if (is.null(base_title)) {
        NULL
      } else {
        paste0(base_title, " (", suffix, ")")
      }

      if (kind == "heatmap") {
        plots[[version]][[kind]] <- .build_loading_heatmap(
          dat = data[[version]],
          title = title,
          fill_limit = limit,
          control = control
        )
      } else {
        plots[[version]][[kind]] <- .build_loading_panel(
          dat = data[[version]],
          color = control[[paste0(version, "_color")]],
          title = title,
          xlab = control$xlab,
          ylab = control$ylab,
          labels = labels,
          axis_limit = 1.2 * limit,
          display = kind,
          control = control
        )
      }
    }
  }

  # One display: non-private left, private right.
  # Multiple displays: non-private row above private row.
  single <- length(display) == 1L
  combined <- patchwork::wrap_plots(
    c(unname(plots$nonprivate), unname(plots$private)),
    ncol = if (single) 2L else length(display),
    nrow = if (single) 1L else 2L,
    byrow = TRUE,
    guides = "collect"
  ) & ggplot2::theme(legend.position = "right")

  attr(combined, "axes") <- axes
  attr(combined, "display") <- display
  attr(combined, "heatmap_components") <- heatmap_components

  point_kind <- intersect(c("vector", "point"), display)
  attr(combined, "selected_variables") <- if (length(point_kind)) {
    selected[[point_kind[1]]]
  } else {
    selected$heatmap
  }

  if ("heatmap" %in% display) {
    attr(combined, "selected_heatmap_variables") <- selected$heatmap
  }

  list(
    loading = loading,
    plot = c(plots, list(all = combined))
  )
}

.make_loading_data <- function(
    V,
    axes,
    selected,
    variable_names
) {
  data.frame(
    variable = variable_names[selected],
    x = V[selected, axes[1]],
    y = V[selected, axes[2]],
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}

.make_loading_heatmap_data <- function(
    V,
    components,
    selected,
    variable_names
) {
  data.frame(
    variable = factor(
      rep(variable_names[selected], times = length(components)),
      levels = rev(variable_names[selected])
    ),
    component = factor(
      rep(paste0("PC", components), each = length(selected)),
      levels = paste0("PC", components)
    ),
    loading = as.numeric(V[selected, components, drop = FALSE])
  )
}

.build_loading_panel <- function(
    dat,
    color,
    title,
    xlab,
    ylab,
    labels,
    axis_limit,
    display,
    control
) {
  plot <- ggplot2::ggplot(dat) +
    ggplot2::geom_hline(
      yintercept = 0,
      color = "grey80",
      linewidth = 0.4,
      linetype = "dashed"
    ) +
    ggplot2::geom_vline(
      xintercept = 0,
      color = "grey80",
      linewidth = 0.4,
      linetype = "dashed"
    )

  if (display == "vector") {
    plot <- plot + ggplot2::geom_segment(
      ggplot2::aes(x = 0, y = 0, xend = .data$x, yend = .data$y),
      color = color,
      linewidth = control$arrow_size,
      lineend = "round",
      arrow = grid::arrow(
        length = grid::unit(0.18, "cm"),
        type = "closed"
      )
    )
  } else {
    plot <- plot + ggplot2::geom_point(
      ggplot2::aes(x = .data$x, y = .data$y),
      color = color,
      size = control$point_size,
      alpha = 0.9
    )
  }

  plot <- plot +
    ggplot2::coord_equal(
      xlim = c(-axis_limit, axis_limit),
      ylim = c(-axis_limit, axis_limit),
      expand = FALSE
    ) +
    ggplot2::labs(title = title, x = xlab, y = ylab) +
    ggplot2::theme_classic(base_size = control$base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        size = control$title_size,
        face = "bold",
        hjust = 0.5
      ),
      panel.border = ggplot2::element_rect(
        color = "grey82",
        fill = NA,
        linewidth = 0.5
      )
    )

  if (labels) {
    plot <- plot + ggrepel::geom_text_repel(
      ggplot2::aes(x = .data$x, y = .data$y, label = .data$variable),
      color = color,
      size = control$label_size,
      max.overlaps = Inf,
      box.padding = 0.35,
      point.padding = 0.2,
      min.segment.length = 0,
      segment.color = "grey65",
      segment.size = 0.35,
      seed = 123
    )
  }

  plot
}

.build_loading_heatmap <- function(
    dat,
    title,
    fill_limit,
    control
) {
  dat$loading_label <- formatC(
    dat$loading,
    format = "f",
    digits = control$heatmap_digits
  )
  dat$label_color <- ifelse(
    abs(dat$loading) >= 0.55 * fill_limit,
    "white",
    "grey20"
  )

  plot <- ggplot2::ggplot(
    dat,
    ggplot2::aes(
      x = .data$component,
      y = .data$variable,
      fill = .data$loading
    )
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.6) +
    ggplot2::scale_fill_gradient2(
      low = control$heatmap_colors[["low"]],
      mid = "white",
      high = control$heatmap_colors[["high"]],
      midpoint = 0,
      limits = c(-fill_limit, fill_limit),
      name = "Loading"
    )

  if (control$heatmap_values) {
    plot <- plot +
      ggplot2::geom_text(
        ggplot2::aes(label = .data$loading_label, color = .data$label_color),
        size = control$heatmap_text_size,
        show.legend = FALSE
      ) +
      ggplot2::scale_color_identity()
  }

  plot +
    ggplot2::labs(title = title, x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = control$base_size) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(
        size = control$title_size,
        face = "bold",
        hjust = 0.5,
        margin = ggplot2::margin(b = 6)
      ),
      axis.text.x = ggplot2::element_text(face = "bold"),
      axis.text.y = ggplot2::element_text(margin = ggplot2::margin(r = 4)),
      legend.title = ggplot2::element_text(face = "bold")
    )
}
