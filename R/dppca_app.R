#' Launch the dppca Shiny app
#'
#' @description
#' Opens an interactive Shiny app for exploring \code{dp_scree_plot()} and
#' \code{dp_score_plot()} with the updated scree and score method options.
#' If \code{X} is supplied, the app opens using that
#' data. If \code{group} is supplied, the score plot uses the grouped version
#' when the group labels are available. If \code{group = "group"}, the
#' column named \code{group} is used as group labels and removed from the
#' PCA feature matrix. In the app, the scree menu uses checkboxes and can overlay
#' \code{method = "clipped"}, \code{method = "pmwm"}, \code{method = "huber"},
#' or a vector such as \code{method = c("clipped", "pmwm", "huber")}. The
#' score menu can pass \code{method = "add"},
#' \code{method = "sparse"}, or \code{method = c("add", "sparse")} and
#' bin counts through \code{bins = c(m_x, m_y)} to the score plotting
#' functions. The app uses a taller, centered scree plot panel
#' so that the scree curves are not visually flattened on wide screens.
#'
#' @param X Optional numeric matrix or data frame. If supplied, the app opens
#'   with this data selected as the default data source.
#' @param group Optional group labels. This can be either a vector of length
#'   \code{nrow(X)} or a single column name in \code{X}.
#'
#' @return Invisibly launches a Shiny application.
#'
#' @examples
#' \dontrun{
#' dppca_app()
#' dppca_app(gau)
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
  if (app_dir == "") {
    stop("Could not find the Shiny app directory.", call. = FALSE)
  }

  shiny::runApp(app_dir, display.mode = "normal")
}
