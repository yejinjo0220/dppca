#' Launch the dppca Shiny app
#'
#' @description
#' Launches an interactive Shiny application for exploring differentially
#' private PCA visualizations. The app provides interfaces for
#' [dp_scree_plot()] and [dp_score_plot()], including private scree plot
#' methods and histogram-based private score plot methods.
#'
#' If `X` is supplied, the app starts with `X` as the initial dataset. If
#' `group` is supplied, the score plot can use the group labels for coloring.
#' The argument `group` can be either a vector of length `nrow(X)` or the name
#' of a column in `X`. When `group` is a column name, that column is used as
#' group labels and removed from the PCA feature matrix.
#'
#' @param X Optional numeric matrix or data frame. If supplied, the app opens
#'   with this data as the initial data source.
#' @param group Optional group labels. This can be either a vector of length
#'   `nrow(X)` or a single column name in `X`.
#'
#' @return No return value. This function is called for its side effect of
#'   launching a Shiny application.
#'
#' @examples
#' if (interactive()) {
#'   dppca_app()
#'
#'   data(gau, package = "dppca")
#'   dppca_app(gau)
#'
#'   data(gau_g, package = "dppca")
#'   dppca_app(gau_g, group = "group")
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
