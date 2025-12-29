#' Run DP PCA score-plot Shiny app for user-supplied data
#'
#' @param X Numeric matrix/data.frame (n x p).
#' @param G Optional group vector (length n) or a column name in X.
#' @param title App title.
#' @export
run_dp_app_data <- function(X, G = NULL, title = "DP PCA Score Plot (User data)") {

  if (missing(X)) stop("X is required.", call. = FALSE)

  X <- as.data.frame(X)

  # Allow G to be a column name in X
  if (is.character(G) && length(G) == 1L && G %in% names(X)) {
    group_vec <- X[[G]]
    X <- X[, setdiff(names(X), G), drop = FALSE]
  } else {
    group_vec <- G
  }

  X_mat <- as.matrix(X)
  if (!is.numeric(X_mat)) stop("X must be numeric.", call. = FALSE)
  if (ncol(X_mat) < 2) stop("X must have at least 2 numeric columns.", call. = FALSE)

  n <- nrow(X_mat)
  if (!is.null(group_vec) && length(group_vec) != n) {
    stop("If provided, length(G) must equal nrow(X).", call. = FALSE)
  }

  # ---- helpers (keep inside so it ships cleanly) ----
  parse_ratio <- function(txt) {
    txt <- gsub("\\s+", "", txt)
    if (txt == "") return(NULL)
    parts <- strsplit(txt, ",")[[1]]
    parts <- parts[parts != ""]
    if (!length(parts)) return(NULL)
    vals <- suppressWarnings(as.numeric(parts))
    if (any(is.na(vals))) return(NULL)
    vals
  }

  # ---- UI ----
  ui <- shiny::fluidPage(
    shiny::titlePanel(title),

    shiny::sidebarLayout(
      shiny::sidebarPanel(

        shiny::checkboxInput("center", "Center variables (center = TRUE)", value = TRUE),
        shiny::checkboxInput("scale",  "Scale to unit variance (scale. = TRUE)", value = FALSE),
        shiny::checkboxInput("dp_pca_flag", "Use DP PCA", value = FALSE),

        shiny::tags$hr(),

        shiny::helpText("default is eps_total = 3, delta_total ≈ 1 / 10^{ceil(log10 n)}."),
        shiny::numericInput("eps", "Total epsilon (eps_total):", value = 4, min = 0.0001, step = 0.1),
        shiny::numericInput("delta", "Total delta (delta_total):", value = 1e-4, min = 0, step = 1e-5),

        shiny::textInput("eps_ratio", "eps_ratio (PCA, Quantile, Hist): ex. 0.3, 0.2, 0.5", value = ""),
        shiny::textInput("delta_ratio", "delta_ratio (PCA, Quantile, Hist):", value = ""),

        shiny::tags$hr(),

        shiny::selectInput(
          "bin_method", "Bin method:",
          choices  = c("J (Jing rule)"="J", "W (Wasserman-Lei rule)"="W", "Manual (use m_x, m_y)"="none"),
          selected = "J"
        ),

        shiny::conditionalPanel(
          condition = "input.bin_method == 'none'",
          shiny::numericInput("m_x", "Number of bins in x (m_x):", value = 10, min = 1, step = 1),
          shiny::numericInput("m_y", "Number of bins in y (m_y):", value = 10, min = 1, step = 1)
        ),

        shiny::tags$hr(),

        shiny::selectInput(
          "mechanism", "Mechanism:",
          choices  = c("all", "none", "add", "sparse"),
          selected = "all"
        ),

        shiny::checkboxInput("sampling", "Draw synthetic samples (sampling)", value = TRUE)
      ),

      shiny::mainPanel(
        shiny::plotOutput("dpPlot", height = "700px")
      )
    )
  )

  # ---- server ----
  server <- function(input, output, session) {

    # set default eps/delta when app starts (based on n)
    shiny::observeEvent(TRUE, {
      pow <- ceiling(log10(n))
      if (!is.finite(pow)) pow <- 2
      delta_default <- 10^(-pow)
      shiny::updateNumericInput(session, "delta", value = delta_default)
      # eps는 네가 4로 쓰고 있으니 유지, 3으로 바꾸고 싶으면 여기서 update 하면 됨
    }, once = TRUE)

    hist_res_reactive <- shiny::reactive({

      shiny::req(!is.na(input$eps), input$eps > 0)
      shiny::req(!is.na(input$delta), input$delta >= 0)

      eps_ratio   <- parse_ratio(input$eps_ratio)
      delta_ratio <- parse_ratio(input$delta_ratio)

      m_x <- if (input$bin_method == "none") input$m_x else NULL
      m_y <- if (input$bin_method == "none") input$m_y else NULL

      common_args <- list(
        center       = input$center,
        scale.       = input$scale,
        dp_pca_flag  = input$dp_pca_flag,
        cpp.option   = FALSE,
        eps_total    = input$eps,
        delta_total  = input$delta,
        eps_ratio    = eps_ratio,
        delta_ratio  = delta_ratio,
        inflate      = 0.10,
        q_frame      = NULL,
        m_x          = m_x,
        m_y          = m_y,
        bin_method   = input$bin_method,
        mechanism    = input$mechanism,
        sampling     = input$sampling
      )

      if (is.null(group_vec)) {
        do.call(
          dppca::dp_score_plot,
          c(list(X = X_mat, color = "#6A5ACD"), common_args)
        )
      } else {
        do.call(
          dppca::dp_score_plot_group,
          c(list(X = X_mat, G = group_vec), common_args)
        )
      }
    })

    # update bin label using current result
    shiny::observe({
      res <- hist_res_reactive()
      if (is.null(res$breaks) || is.null(res$breaks$x) || is.null(res$breaks$y)) return()

      mx <- length(res$breaks$x) - 1
      my <- length(res$breaks$y) - 1

      if (input$bin_method == "J") {
        label_J <- sprintf("J (Jing rule / %d × %d bins)", mx, my)
        label_W <- "W (Wasserman–Lei rule)"
      } else if (input$bin_method == "W") {
        label_J <- "J (Jing rule)"
        label_W <- sprintf("W (Wasserman–Lei rule / %d × %d bins)", mx, my)
      } else {
        label_J <- "J (Jing rule)"
        label_W <- "W (Wasserman–Lei rule)"
      }

      shiny::updateSelectInput(
        session, "bin_method",
        choices = setNames(
          c("J","W","none"),
          c(label_J, label_W, "Manual (use m_x, m_y)")
        ),
        selected = input$bin_method
      )
    })

    output$dpPlot <- shiny::renderPlot({
      hist_res_reactive()$plot$all
    })
  }

  shiny::runApp(shiny::shinyApp(ui, server), display.mode = "normal")
}
