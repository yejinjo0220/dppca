#' Run the DP PCA Shiny App
#'
#' Launches the interactive Shiny application for DP PCA visualization.
#'
#' @export
run_dp_score_plot <- function() {
  app_dir <- system.file("shiny-apps", "dp_score_plot", package = "dppca")
  if (app_dir == "") {
    stop("Could not find Shiny app directory. Try reinstalling the package.", call. = FALSE)
  }
  shiny::runApp(app_dir, display.mode = "normal")
}
