# ============================================================
# loading.R
# Private PC direction estimation and loading visualization
# ============================================================

#' Estimate differentially private PC directions
#'
#' Estimate principal component (PC) directions using a noisy spherical Kendall
#' matrix, with optional ordinary PCA directions for comparison.
#'
#' @param X A numeric matrix or data frame with observations in rows and variables
#'   in columns. Non-numeric data-frame columns are excluded with a warning;
#'   numeric values must be finite. At least two observations and one numeric
#'   feature are required. Observations must be sufficiently separated to compute
#'   normalized pairwise differences.
#' @param privacy A [privacy_control()] object allocating privacy budget to
#'   `loading` alone. Either specify total `eps` and `delta` or an absolute
#'   `loading` budget.
#' @param cpp.option Whether to use the compiled spherical-Kendall implementation
#'   (`TRUE` by default); `FALSE` uses the R implementation.
#' @param non_private Whether to compute and store ordinary PCA reference results
#'   (`FALSE` by default).
#'
#' @return A `dp_loading` object with the following components:
#'   \describe{
#'     \item{`directions`}{A `p` by `p` matrix of private PC directions.}
#'     \item{`nonprivate`}{A `p` by `p` matrix of ordinary PC directions,
#'       or `NULL` when `non_private = FALSE`.}
#'     \item{`metadata`}{Sample size, variable count and privacy parameters.}
#'   }
#'   Rows represent variables, and columns represent `PC1`, ..., `PCp`.
#'   The columns are orthonormal eigenvectors, without eigenvalue or
#'   correlation scaling.
#'
#' @details
#' The private PC directions are eigenvectors of a spherical Kendall matrix
#' perturbed with symmetric Gaussian noise. Directions are ordered by the noisy
#' eigenvalues; they need not be eigenvectors of the sample covariance matrix.
#' The spherical Kendall construction uses normalized pairwise differences and
#' is invariant to common translations.
#'
#' With `non_private = TRUE`, ordinary covariance-based directions are stored
#' after sign alignment with the private directions. This comparison is not a
#' private release.
#'
#' Under fixed-size replacement adjacency, the matrix has Frobenius sensitivity
#' at most `4 / n`. The diagonal Gaussian noise standard deviation is
#' `4 * sqrt(2 * log(1.25 / delta)) / (n * eps)`; off-diagonal entries use this
#' scale divided by `sqrt(2)`. This classical calibration has a standard
#' sufficient guarantee for `0 < eps < 1`; larger `eps` requires an additional
#' calibration check, which is not performed by the function.
#'
#' Duplicate or numerically indistinguishable observations are rejected, but
#' these checks do not privatize input validation. Data schema and sample size
#' are assumed public. Data-dependent preprocessing requires separate privacy
#' analysis and accounting.
#'
#' @examples
#' set.seed(123)
#' X <- matrix(rnorm(150), nrow = 50, ncol = 3)
#' fit <- dp_loading(
#'   X, privacy = privacy_control(eps = 0.5, delta = 1e-3),
#'   cpp.option = FALSE
#' )
#' print(fit)
#' summary(fit)
#' fit$directions
#' fit_with_reference <- dp_loading(
#'   X, privacy = privacy_control(loading = c(eps = 0.5, delta = 1e-3)),
#'   cpp.option = FALSE, non_private = TRUE
#' )
#' fit_with_reference$nonprivate
#' \dontrun{
#' loading_plot <- plot(fit, type = c("vector", "point", "heatmap"),
#'                      heatmap_components = 1:3, variables = "all")
#' loading_plot                 # Display the saved plot
#' loading_plot$loading$directions # Access the stored directions
#' loading_plot$plot$all       # Access the combined plot
#' }
#' @export
dp_loading <- function(
    X,
    privacy,
    cpp.option = TRUE,
    non_private = FALSE
) {
  # Validate the supplied data without centering or standardizing it.
  X <- .loading_matrix(X)

  if (missing(privacy)) {
    stop("Supply `privacy = privacy_control(...)`.", call. = FALSE)
  }

  # The implicit default split is for integrated PCA, not for dp_loading().
  implicit <- isTRUE(attr(privacy, "implicit_split"))
  required <- if (implicit) c("loading", "scree", "score") else "loading"
  checked <- .validate_privacy_control(privacy, required)
  budget <- if (implicit) checked$total else checked$budgets$loading
  eps <- budget[["eps"]]
  delta <- budget[["delta"]]

  .loading_positive(eps, "eps")
  .loading_positive(delta, "delta")

  if (delta >= 1) {
    stop("`delta` must be less than 1.", call. = FALSE)
  }

  validate_logical_value(cpp.option, "cpp.option")
  validate_logical_value(non_private, "non_private")

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

  # Estimate private directions and the optional ordinary reference.
  V_private <- eigen(kendall + noise, symmetric = TRUE)$vectors
  V_nonprivate <- if (non_private) {
    V <- eigen(stats::cov(X), symmetric = TRUE)$vectors

    # Align the non-private PC signs with the private directions before storage.
    signs <- ifelse(colSums(V * V_private) < 0, -1, 1)
    sweep(V, 2L, signs, "*")
  } else {
    NULL
  }

  loading_names <- list(colnames(X), paste0("PC", seq_len(p)))

  dimnames(V_private) <- loading_names
  if (non_private) {
    dimnames(V_nonprivate) <- loading_names
  }

  structure(
    list(
      directions = V_private,
      nonprivate = V_nonprivate,
      metadata = list(n = n, p = p, privacy = c(eps = eps, delta = delta))
    ),
    class = "dp_loading"
  )
}

#' Print private PC direction results
#'
#' Display a concise overview of a [dp_loading()] fit, including its data
#' dimensions, privacy parameters, and whether a non-private reference is stored.
#'
#' @param x A `dp_loading` or `summary.dp_loading` object, depending
#'   on the S3 method.
#' @param ... Additional arguments, reserved for S3 generic compatibility.
#' @return Invisibly returns `x`.
#' @rdname dp_loading
#' @export
print.dp_loading <- function(x, ...) {
  cat("Differentially Private PC Directions\n")
  if (!is.null(x$metadata$n)) {
    cat("Data:", x$metadata$n, "observations x",
        nrow(x$directions), "variables\n")
  } else {
    cat("Variables:", nrow(x$directions), "\n")
  }
  cat("Components:", ncol(x$directions), "PC directions\n")
  privacy <- x$metadata$privacy
  if (!is.null(privacy)) {
    cat("Total privacy budget: eps =", format(privacy[["eps"]], digits = 6),
        "; delta =", format(privacy[["delta"]], digits = 6), "\n")
  }
  cat("Non-private reference:",
      if (is.null(x$nonprivate)) "not stored" else "included", "\n")
  invisible(x)
}

#' Summarize private PC directions
#'
#' Summarize a fitted [dp_loading()] object without recomputing its estimates
#' or consuming additional privacy budget.
#'
#' @param object A `dp_loading` object.
#' @param ... Additional arguments, reserved for S3 generic compatibility.
#' @return A `summary.dp_loading` object containing the stored `private` and
#'   `nonprivate` direction matrices, with sample size and privacy budget
#'   metadata. The print method displays up to six variables and four PCs;
#'   the full matrices remain stored. The non-private matrix remains `NULL`
#'   if it was not stored.
#' @rdname dp_loading
#' @export
summary.dp_loading <- function(object, ...) {
  structure(
    list(
      directions = object$directions,
      nonprivate = object$nonprivate,
      metadata = object$metadata
    ),
    class = "summary.dp_loading"
  )
}

#' Print summarized PC directions
#'
#' @param x A `summary.dp_loading` object.
#' @param digits Number of significant digits used for printing.
#' @param ... Additional arguments.
#' @return Invisibly returns `x`.
#' @rdname dp_loading
#' @export
print.summary.dp_loading <- function(x, digits = 3, ...) {
  cat("DP PC Directions Summary\n")
  if (!is.null(x$metadata$n)) {
    cat("Data:", x$metadata$n, "observations x",
        nrow(x$directions), "variables\n")
  } else {
    cat("Variables:", nrow(x$directions), "\n")
  }
  cat("Components:", ncol(x$directions), "PC directions\n")

  rows <- seq_len(min(6L, nrow(x$directions)))
  cols <- seq_len(min(4L, ncol(x$directions)))
  caption <- paste0(" (first ", length(rows), " variables x first ",
                    length(cols), " PCs):\n")
  cat("\nPrivate PC directions", caption)
  print(x$directions[rows, cols, drop = FALSE], digits = digits, ...)
  if (!is.null(x$nonprivate)) {
    cat("\nNon-private PC directions", caption)
    print(x$nonprivate[rows, cols, drop = FALSE], digits = digits, ...)
  }

  privacy <- x$metadata$privacy
  if (!is.null(privacy)) {
    cat("\nTotal privacy budget: eps =", format(privacy[["eps"]], digits = 6),
        "; delta =", format(privacy[["delta"]], digits = 6), "\n")
  }
  cat("Non-private reference:",
      if (is.null(x$nonprivate)) "not stored" else "included", "\n")
  cat(if (!is.null(x$nonprivate)) {
    "Full direction matrices are available in $directions and $nonprivate.\n"
  } else {
    "The full direction matrix is available in $directions.\n"
  })
  invisible(x)
}

#' Plot stored private PC directions
#'
#' Create loading plots from an existing [dp_loading()] fit. Plotting uses the
#' stored direction matrices and does not re-estimate PCs or spend additional
#' privacy budget.
#'
#' @param x A `dp_loading` object.
#' @param axes Two distinct PC indices for vector and point plots, displayed
#'   on fixed axes from -1.2 to 1.2. Default: `c(1, 2)`.
#' @param type One or more of `"vector"`, `"point"`, and `"heatmap"`.
#'   The default shows all three types. With only one stored component, the
#'   default shows a heatmap because vector and point plots need two axes.
#' @param heatmap_components PC indices to include in a heatmap; `NULL`
#'   (default) uses `axes`.
#' @param variables Variables to display. A single positive integer selects
#'   the variables with the largest Euclidean norms across the displayed PC
#'   coefficients (top 10 by default). Use `"all"` or `NULL` for every
#'   variable, or provide variable names or a vector of column indices.
#'   Selection is performed separately for each version and display.
#' @param plot_control Optional appearance settings from
#'   [loading_plot_control()] or a named list of its options.
#' @param ... Additional arguments, reserved for S3 generic compatibility.
#' @return Draws the combined plot and invisibly returns a `dp_loading_plot`
#'   object with `loading` (the original fit) and `plot` (individual plots and
#'   the combined `plot$all`). Printing the saved object displays the plot.
#' @details
#' Vector and point plots use fixed axes from -1.2 to 1.2 and show variable
#' labels. Each heatmap scales colors over its own displayed coefficients,
#' so comparable coefficients can have different colors between versions.
#' An ordinary PCA reference is shown automatically only when it was stored
#' at fitting; it is not a private release.
#'
#' @export
plot.dp_loading <- function(
    x,
    axes = c(1, 2),
    type = c("vector", "point", "heatmap"),
    heatmap_components = NULL,
    variables = 10,
    plot_control = NULL,
    ...
) {
  if (length(list(...))) {
    stop("Use `type`, the named loading plot options or `plot_control`.", call. = FALSE)
  }
  if (missing(type) && ncol(x$directions) == 1L) {
    type <- "heatmap"
  }
  result <- .plot_loading_result(
    loading = x,
    axes = axes,
    display = type,
    heatmap_components = heatmap_components,
    variables = variables,
    plot_control = plot_control
  )
  print(result$plot$all)
  invisible(result)
}

#' Print a saved loading plot
#'
#' Draw the combined plot stored in a `dp_loading_plot` object. The original
#' directions remain available in `x$loading`.
#'
#' @param x A `dp_loading_plot` object.
#' @param ... Additional arguments, reserved for S3 generic compatibility.
#' @return Invisibly returns `x`.
#' @rdname plot.dp_loading
#' @export
print.dp_loading_plot <- function(x, ...) {
  print(x$plot$all)
  invisible(x)
}

#' Summarize a loading plot result
#'
#' Summarize the fitted PC directions stored in a `dp_loading_plot` object.
#'
#' @param object A `dp_loading_plot` object.
#' @param ... Arguments passed to [summary.dp_loading()].
#' @return A `summary.dp_loading` object.
#' @rdname plot.dp_loading
#' @export
summary.dp_loading_plot <- function(object, ...) {
  summary(object$loading, ...)
}

#' Redraw or customize a loading plot result
#'
#' Display a previously constructed plot or redraw it using new plotting
#' options without re-estimating the PC directions.
#'
#' @param x A `dp_loading_plot` object.
#' @param ... Options passed to [plot.dp_loading()]. With no additional
#'   options, redraw the stored plot using its original settings.
#' @return Invisibly returns a `dp_loading_plot` object.
#' @rdname plot.dp_loading
#' @export
plot.dp_loading_plot <- function(x, ...) {
  if (!length(list(...))) {
    print(x$plot$all)
    return(invisible(x))
  }
  plot(x$loading, ...)
}

#' Control loading plot appearance
#'
#' Specify colors, text sizes, labels, and heatmap settings for loading plots.
#' Pass the resulting list through the `plot_control` argument of
#' [plot.dp_loading()].
#'
#' @param nonprivate_color,private_color Colors of the non-private and private
#'   vectors or points, respectively.
#' @param arrow_size Width of arrows in vector plots.
#' @param point_size Size of points in point plots.
#' @param label_size Size of variable labels in vector and point plots.
#' @param heatmap_values Whether to display loading coefficients in heatmap
#'   cells. Default: `TRUE`.
#' @param heatmap_digits Number of decimal places for heatmap labels.
#' @param heatmap_text_size Size of coefficient labels in heatmaps.
#' @param base_size Base font size for all loading plots.
#' @param title_size Font size for panel titles.
#' @param nonprivate_title,private_title Text used for the corresponding panel
#'   titles. The display name is appended in parentheses; `NULL` hides a title.
#' @param xlab,ylab Axis labels for vector and point plots. `NULL` uses the
#'   selected PC names.
#' @param heatmap_colors Named character vector with `low` and `high` colors
#'   for negative and positive loading coefficients. White represents zero.
#'   Private and non-private heatmaps use separate, symmetric color ranges
#'   based on their own displayed coefficients.
#' @return A named list of validated loading plot settings.
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

  validate_logical_value(heatmap_values, "heatmap_values")

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

#' Select numeric PCA variables without centering or scaling
#'
#' @param X A real numeric matrix or a data frame.
#' @param group NULL, a column name, or one group label per observation.
#' @param warn Whether to report excluded columns and the selected group.
#' @param min_cols Minimum number of numeric PCA variables required.
#' @return A list containing `X`, `group`, `group_name`, `group_levels`,
#'   and `excluded_columns`. Column types, names, and group labels are
#'   assumed public; the warning is not a private release mechanism.
#' @noRd
.prepare_pca_input <- function(X, group = NULL, warn = TRUE, min_cols = 1L) {
  validate_logical_value(warn, "warn")
  if (!is.numeric(min_cols) || is.complex(min_cols) || length(min_cols) != 1L ||
      !is.finite(min_cols) || min_cols < 1L || min_cols != floor(min_cols)) {
    stop("`min_cols` must be a positive integer.", call. = FALSE)
  }
  if (!is.data.frame(X) && !is.matrix(X)) {
    stop("`X` must be a numeric matrix or data frame.", call. = FALSE)
  }
  if (nrow(X) < 2L) {
    stop("`X` must have at least two rows.", call. = FALSE)
  }
  if (is.matrix(X) && (!is.numeric(X) || is.complex(X))) {
    stop(
      "`X` must be a real numeric matrix. For mixed columns, supply a ",
      "data frame so numeric columns can be selected.",
      call. = FALSE
    )
  }

  column_names <- colnames(X)
  if (is.null(column_names)) {
    column_names <- if (ncol(X)) paste0("X", seq_len(ncol(X))) else character()
    colnames(X) <- column_names
  }
  if (anyNA(column_names) || any(!nzchar(column_names)) ||
      anyDuplicated(column_names)) {
    stop("`X` must have unique, non-missing column names.", call. = FALSE)
  }
  numeric_cols <- if (is.data.frame(X)) {
    vapply(X, function(column) {
      is.numeric(column) && !is.complex(column) && is.null(dim(column))
    }, logical(1))
  } else {
    rep(TRUE, ncol(X))
  }
  group_values <- group_name <- group_levels <- NULL
  group_column <- NULL

  if (!is.null(group)) {
    if (is.character(group) && length(group) == 1L && !is.na(group)) {
      if (!group %in% column_names) {
        stop("Group column `", group, "` was not found in `X`.", call. = FALSE)
      }
      group_column <- group
      group_name <- group
      group_values <- if (is.data.frame(X)) X[[group]] else X[, group]
    } else {
      group_name <- "<supplied vector>"
      group_values <- group
      # Recover a useful label when the supplied vector uniquely matches a
      # column. This changes metadata only, not feature selection.
      matching_columns <- column_names[vapply(seq_along(column_names), function(i) {
        column <- if (is.data.frame(X)) X[[i]] else X[, i]
        identical(unname(column), unname(group_values))
      }, logical(1))]
      if (length(matching_columns) == 1L) {
        group_name <- matching_columns
      }
    }
    if (!is.atomic(group_values) || !is.null(dim(group_values)) ||
        is.complex(group_values) || length(group_values) != nrow(X) ||
        anyNA(group_values)) {
      stop(
        "`group` must select a column or supply one non-missing label per row.",
        call. = FALSE
      )
    }
    group_values <- as.character(group_values)
    if (anyNA(group_values) || any(!nzchar(trimws(group_values)))) {
      stop("Group labels must be non-missing and non-empty.", call. = FALSE)
    }
    group_levels <- unique(group_values)
  }

  # A selected group is a label, even if the column stores numeric codes.
  keep <- numeric_cols & !(column_names %in% group_column)
  excluded_columns <- column_names[!keep]
  if (!any(keep)) {
    stop("No numeric PCA variables remain after excluding non-numeric and group columns.",
         call. = FALSE)
  }
  if (sum(keep) < min_cols) {
    stop("At least ", min_cols, " numeric PCA variables are required after exclusions.",
         call. = FALSE)
  }
  X <- as.matrix(X[, keep, drop = FALSE])
  storage.mode(X) <- "double"
  if (anyNA(X) || any(!is.finite(X))) {
    stop("Numeric PCA variables must contain only finite, non-missing values.",
         call. = FALSE)
  }

  if (warn && (length(excluded_columns) || !is.null(group_values))) {
    messages <- c(
      paste0("PCA: ", ncol(X), " numeric variables used."),
      paste0("Excluded from PCA (", length(excluded_columns), "): ",
             if (length(excluded_columns)) paste(excluded_columns, collapse = ", ") else "none",
             ".")
    )
    if (!is.null(group_values)) {
      messages <- c(messages,
        paste0("Grouping variable: ", group_name),
        paste0("Groups (", length(group_levels), "): ",
               paste(group_levels, collapse = ", "))
      )
    }
    warning(paste(messages, collapse = "\n"), call. = FALSE)
  }

  list(
    X = X,
    group = group_values,
    group_name = group_name,
    group_levels = group_levels,
    excluded_columns = excluded_columns
  )
}

#' Prepare the shared numeric matrix used by estimation routines
#' @inheritParams .prepare_pca_input
#' @return A double matrix of the selected numeric variables.
#' @noRd
prep_matrix_for_pca <- function(X) {
  .prepare_pca_input(X)$X
}


.loading_positive <- function(
    x,
    name
) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0) {
    stop("`", name, "` must be a positive finite number.", call. = FALSE)
  }
}

.loading_matrix <- function(X) {
  X <- prep_matrix_for_pca(X)

  if (is.null(colnames(X))) {
    colnames(X) <- paste0("X", seq_len(ncol(X)))
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

#' Validate a full external loading matrix
#'
#' @param V_dp Finite real `p` by `p` matrix with orthonormal columns.
#' @param X Input data with `p` variables in the same order as the rows
#'   of `V_dp`. Names are checked when both inputs have variable names.
#' @return Validated directions named by the variables and `PC1`, ..., `PCp`.
#' @details Shape and orthogonality checks do not verify privacy provenance.
#' @noRd
.validate_dp_directions <- function(V_dp, X) {
  V_dp <- as.matrix(V_dp)
  p <- ncol(X)

  if (!is.numeric(V_dp) || is.complex(V_dp) ||
      !identical(dim(V_dp), c(p, p)) ||
      any(!is.finite(V_dp))) {
    stop(
      "V_dp must be a finite real p by p matrix.",
      call. = FALSE
    )
  }

  if (max(abs(crossprod(V_dp) - diag(p))) > 1e-6) {
    stop(
      "Columns of V_dp must be orthonormal.",
      call. = FALSE
    )
  }

  if (!is.null(rownames(V_dp)) &&
      !is.null(colnames(X)) &&
      !identical(rownames(V_dp), colnames(X))) {
    stop(
      "Rows of V_dp must match the columns of X in order.",
      call. = FALSE
    )
  }

  dimnames(V_dp) <- list(
    colnames(X),
    paste0("PC", seq_len(p))
  )

  V_dp
}


#' Resolve shared private directions
#'
#' @param X Validated input data matrix.
#' @param budget Loading allocation as a named `c(eps, delta)` vector.
#' @param cpp.option Use the compiled implementation when estimating loadings.
#' @param V_dp Optional full direction matrix on the same variable scales as
#'   `X`. When supplied, estimation is skipped and `budget` records its prior
#'   cost; the caller is responsible for its privacy provenance.
#' @return A validated full direction matrix. No ordinary PCA is computed.
#' @noRd
.get_dp_directions <- function(
    X, budget, cpp.option, V_dp = NULL
) {
  if (is.null(V_dp)) {
    V_dp <- dp_loading(
      X,
      privacy = privacy_control(loading = budget),
      cpp.option = cpp.option,
      non_private = FALSE
    )$directions
  }

  .validate_dp_directions(V_dp, X)
}

.nonprivate_pca <- function(X) {
  eig <- eigen(stats::cov(X), symmetric = TRUE)

  dimnames(eig$vectors) <- list(
    colnames(X),
    paste0("PC", seq_len(ncol(X)))
  )

  list(
    directions = eig$vectors,
    eigenvalues = pmax(eig$values, 0)
  )
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
    stop("`type` must contain vector, point, or heatmap.", call. = FALSE)
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

.loading_plot_packages <- function() {
  packages <- c("ggplot2", "patchwork")

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
    plot_control,
    non_private = !is.null(loading$nonprivate)
) {
  V_private <- loading$directions
  p <- ncol(V_private)
  if (anyDuplicated(rownames(V_private))) {
    stop("Variable names must be unique.", call. = FALSE)
  }
  display <- .validate_loading_display(display)
  # A heatmap can display one component, even for a one-variable fit.
  if (identical(display, "heatmap") && p == 1L && identical(axes, c(1, 2))) {
    axes <- c(1L, 1L)
  } else {
    axes <- .validate_loading_axes(axes, p)
  }
  heatmap_components <- .validate_loading_components(
    heatmap_components, p, axes, "heatmap_components"
  )
  validate_logical_value(non_private, "non_private")
  if (non_private && is.null(loading$nonprivate)) {
    stop("No non-private reference was stored. Refit with non_private = TRUE.",
         call. = FALSE)
  }
  .loading_plot_packages()

  V_nonprivate <- if (non_private) loading$nonprivate else NULL
  variable_names <- rownames(V_private)
  if (is.null(variable_names)) {
    variable_names <- paste0("X", seq_len(nrow(V_private)))
  }
  control <- .merge_loading_plot_control(plot_control, axes)

  matrices <- list(private = V_private)
  if (non_private) {
    matrices <- list(nonprivate = V_nonprivate, private = V_private)
  }
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

    # Heatmap limits are computed separately for each version below.
    # Vector and point plots use a fixed [-1.2, 1.2] square.

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
        fill_limit <- max(abs(data[[version]]$loading))
        if (!is.finite(fill_limit) ||
            fill_limit <= sqrt(.Machine$double.eps)) {
          fill_limit <- 1
        }
        plots[[version]][[kind]] <- .build_loading_heatmap(
          dat = data[[version]],
          title = title,
          fill_limit = fill_limit,
          control = control
        )
      } else {
        plots[[version]][[kind]] <- .build_loading_panel(
          dat = data[[version]],
          color = control[[paste0(version, "_color")]],
          title = title,
          xlab = control$xlab,
          ylab = control$ylab,
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
    ncol = if (single) length(matrices) else length(display),
    nrow = if (single) 1L else length(matrices),
    byrow = TRUE,
    guides = "keep"
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

  structure(
    list(loading = loading, plot = c(plots, list(all = combined))),
    class = "dp_loading_plot"
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
      xlim = c(-1.2, 1.2),
      ylim = c(-1.2, 1.2),
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

  if (requireNamespace("ggrepel", quietly = TRUE)) {
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
  } else {
    plot <- plot + ggplot2::geom_text(
      ggplot2::aes(x = .data$x, y = .data$y, label = .data$variable),
      color = color,
      size = control$label_size,
      vjust = -0.6
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
      name = NULL
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
