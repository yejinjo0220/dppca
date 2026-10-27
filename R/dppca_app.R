#' Launch the dppca Shiny app
#'
#' @description
#' Launches an interactive 'shiny' application for exploring differentially
#' private PCA visualizations. The PCA tab uses [dp_pca()] for integrated
#' estimation, summaries, and scree, score, loading, variable, and biplot
#' displays, including grouped score plots when grouping is enabled.
#' The Loading, Scree, and Score tabs provide standalone analyses, including
#' method comparisons and grouped score plots.
#'
#' @param X A numeric matrix or data frame with observations in rows and variables
#'   in columns. Non-numeric data-frame columns are excluded with a warning;
#'   numeric values must be finite. Input is not automatically centered or
#'   standardized. If provided, the app uses this dataset initially. At least
#'   two observations and two numeric features are required.
#' @param group Optional group column name in `X` or vector of group labels
#'   (`NULL` by default). A selected group column is excluded from PCA
#'   features, even when numeric; group labels and membership are assumed
#'   public. Only used with supplied `X`; activate grouping
#'   in the PCA or Score tab to use these labels.
#'
#' @return Invisibly returns the value produced by \code{shiny::runApp()} after
#'   the application exits.
#'
#' @details
#' The app accepts built-in data, supplied `X`, or uploaded CSV files.
#' Grouping is available for PCA and Score analyses; a selected group column is
#' excluded from numeric PCA features. Loading and Scree analyses use all
#' numeric feature columns.
#'
#' Running an analysis stores its estimates for subsequent plotting. Display
#' settings can redraw an existing fit without consuming another privacy budget.
#' Non-private reference output is optional in the app and is not a private
#' release. Budget tables display allocations but do not track cumulative
#' privacy cost across runs.
#'
#' Input is not automatically centered or standardized. Any data-dependent
#' preprocessing requires separate privacy accounting; see [dp_pca()] and
#' [dp_score()] for applicable conditions.
#'
#' @examples
#' if (interactive()) {
#' # Launch the app with built-in example datasets.
#' dppca_app()
#'
#' # Launch the app with a user-supplied numeric dataset.
#' data(gau, package = "dppca")
#' dppca_app(gau)
#'
#' # Launch the app with group labels stored in a column.
#' data(gau_g, package = "dppca")
#' dppca_app(gau_g, group = "group")
#' }
#'
#' @export
dppca_app <- function(X = NULL, group = NULL) {
  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("Package 'shiny' is required to launch the dppca app.", call. = FALSE)
  }

  if (!is.null(X)) {
    if (!(is.data.frame(X) || is.matrix(X))) {
      stop("'X' must be a numeric matrix or data frame.", call. = FALSE)
    }
    if (nrow(X) < 2L) {
      stop("'X' must have at least two rows.", call. = FALSE)
    }
    if (!is.null(group)) {
      if (is.character(group) && length(group) == 1L && is.data.frame(X)) {
        if (!group %in% colnames(X)) {
          stop("Column '", group, "' was not found in X.", call. = FALSE)
        }
      } else if (length(group) != nrow(X)) {
        stop("'group' must be a vector of length nrow(X), or a column name in X.", call. = FALSE)
      }
    }
  }

  old_X <- getOption("dppca.app.X", NULL)
  old_group <- getOption("dppca.app.group", NULL)

  options(dppca.app.X = X, dppca.app.group = group)

  on.exit({
    options(dppca.app.X = old_X, dppca.app.group = old_group)
  }, add = TRUE)

  app_dir <- system.file("shiny/dppca-app", package = "dppca")
  if (identical(app_dir, "")) {
    stop("Could not find the Shiny app directory.", call. = FALSE)
  }

  invisible(shiny::runApp(app_dir, display.mode = "normal"))
}
