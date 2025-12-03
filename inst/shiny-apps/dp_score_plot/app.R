
library(shiny)
library(dppca)

parse_ratio <- function(txt) {
  txt <- gsub("\\s+", "", txt)
  if (txt == "") return(NULL)

  parts <- strsplit(txt, ",")[[1]]
  parts <- parts[parts != ""]

  if (!length(parts)) return(NULL)

  vals <- suppressWarnings(as.numeric(parts))
  if (any(is.na(vals))) {
    return(NULL)
  }
  vals
}


# UI ------------------------------------------------------------------
ui <- fluidPage(
  titlePanel("DP PCA Score Plot"),

  sidebarLayout(
    sidebarPanel(

      # data selection
      selectInput(
        inputId = "dataset",
        label   = "Dataset:",
        choices = c(
          "eur_map (no group)"                  = "eur_map",
          "eur_map (group = color)"             = "eur_map_g",
          "gaussian mixture (no group)"         = "gau",
          "gaussian mixture (group = color)"    = "gau_g",
          "adult (numeric only)"                = "adult"
        ),
        selected = "eur_map"
      ),

      # center / scale option
      checkboxInput(
        inputId = "center",
        label   = "Center variables (center = TRUE)",
        value   = TRUE
      ),
      checkboxInput(
        inputId = "scale",
        label   = "Scale to unit variance (scale. = TRUE)",
        value   = FALSE
      ),

      numericInput(
        inputId = "pc_x",
        label   = "PC index on x-axis:",
        value   = 1,
        min     = 1,
        step    = 1
      ),
      numericInput(
        inputId = "pc_y",
        label   = "PC index on y-axis:",
        value   = 2,
        min     = 1,
        step    = 1
      ),

      checkboxInput(
        inputId = "dp_pca_flag",
        label   = "Use DP PCA",
        value   = FALSE
      ),

      tags$hr(),

      # total eps, delta ----
      helpText("default is eps_total = 3, delta_total ≈ 1 / 10^{ceil(log10 n)}."),
      numericInput(
        inputId = "eps",
        label   = "Total epsilon (eps_total):",
        value   = 4,
        min     = 0.0001,
        step    = 0.1
      ),
      numericInput(
        inputId = "delta",
        label   = "Total delta (delta_total):",
        value   = 1e-4,
        min     = 0,
        step    = 1e-5
      ),

      # ---- (3) eps_ratio / delta_ratio ----
      textInput(
        inputId = "eps_ratio",
        label   = "eps_ratio (PCA, Quantile, Hist): ex. 0.3, 0.2, 0.5",
        value   = ""
      ),
      textInput(
        inputId = "delta_ratio",
        label   = "delta_ratio (PCA, Quantile, Hist):",
        value   = ""
      ),

      tags$hr(),

      # bin method + m_x, m_y
      selectInput(
        inputId = "bin_method",
        label   = "Bin method:",
        choices = c(
          "J (Jing rule)"               = "J",
          "W (Wasserman rule)"      = "W",
          "Manual (use m_x, m_y)"       = "none"
        ),
        selected = "J"
      ),

      conditionalPanel(
        condition = "input.bin_method == 'none'",
        numericInput(
          inputId = "m_x",
          label   = "Number of bins in x (m_x):",
          value   = 10,
          min     = 1,
          step    = 1
        ),
        numericInput(
          inputId = "m_y",
          label   = "Number of bins in y (m_y):",
          value   = 10,
          min     = 1,
          step    = 1
        )
      ),

      tags$hr(),

      # Mechanism and sampling
      selectInput(
        inputId = "mechanism",
        label   = "Mechanism:",
        choices = c("all", "none", "add", "sparse"),
        selected = "all"
      ),

      checkboxInput(
        inputId = "sampling",
        label   = "Draw synthetic samples (sampling)",
        value   = TRUE
      )
    ),

    mainPanel(
      plotOutput("dpPlot", height = "700px")
    )
  )
)

# SERVER
server <- function(input, output, session) {

  # dataset에 따라 X, G 결정 ----
  dataset_reactive <- reactive({
    switch(input$dataset,
           "eur_map" = {
             list(X = as.matrix(eur_map), G = NULL)
           },
           "eur_map_g" = {
             list(
               X = as.matrix(eur_map_g[ , !(names(eur_map_g) %in% "color") ]),
               G = eur_map_g$color
             )
           },
           "gau" = {
             list(X = as.matrix(gau), G = NULL)
           },
           "gau_g" = {
             list(
               X = as.matrix(gau_g[ , !(names(gau_g) %in% "color") ]),
               G = gau_g$color
             )
           },
           "adult" = {
             list(X = as.matrix(adult), G = NULL)
           }
    )
  })

  # dataset에 따라 pc_x, pc_y의 가능한 범위(최대) 설정
  observeEvent(dataset_reactive(), {
    dat <- dataset_reactive()
    X <- dat$X
    p <- ncol(X)
    if (is.null(p) || is.na(p) || p < 2) return()

    # x축: 기본 1
    updateNumericInput(
      session, "pc_x",
      value = min(1, p),
      min   = 1,
      max   = p
    )
    # y축: 기본 2 (단 p=1이면 1로)
    updateNumericInput(
      session, "pc_y",
      value = if (p >= 2) 2 else 1,
      min   = 1,
      max   = p
    )
  })

  # dataset 선택 시 eps_total, delta_total 추천값 설정
  observeEvent(dataset_reactive(), {
    dat <- dataset_reactive()
    X <- dat$X
    n <- nrow(X)
    if (is.null(n) || is.na(n) || n <= 0) return()

    eps_default <- 4

    # delta_default = 1 / 10^{ceil(log10(n))}
    pow <- ceiling(log10(n))
    if (!is.finite(pow)) pow <- 2
    delta_default <- 10^(-pow)

    updateNumericInput(session, "eps",   value = eps_default)
    updateNumericInput(session, "delta", value = delta_default)
  }, ignoreInit = FALSE)

  # DP result
  hist_res_reactive <- reactive({

    # eps, delta가 “입력 중”일 때(NA 등)에는 계산하지 않음
    req(!is.na(input$eps), input$eps > 0)
    req(!is.na(input$delta), input$delta >= 0)

    dat <- dataset_reactive()
    X <- dat$X
    G <- dat$G

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
      inflate      = 0.20,
      q_frame      = NULL,
      m_x          = m_x,
      m_y          = m_y,
      bin_method   = input$bin_method,
      mechanism    = input$mechanism,
      sampling     = input$sampling,
      axes         = c(input$pc_x, input$pc_y)
    )

    if (is.null(G)) {
      do.call(
        dp_score_plot,
        c(list(X = X, color = "#6A5ACD"), common_args)
      )
    } else {
      do.call(
        dp_score_plot_group,
        c(list(X = X, G = G), common_args)
      )
    }
  })

  # bin_method label update
  observe({
    res <- hist_res_reactive()
    if (is.null(res$breaks) || is.null(res$breaks$x) || is.null(res$breaks$y)) {
      return()
    }

    m_x <- length(res$breaks$x) - 1
    m_y <- length(res$breaks$y) - 1

    if (input$bin_method == "J") {
      label_J <- sprintf("J (Jing rule / %d × %d bins)", m_x, m_y)
      label_W <- "W (Wasserman rule)"
    } else if (input$bin_method == "W") {
      label_J <- "J (Jing rule)"
      label_W <- sprintf("W (Wasserman rule / %d × %d bins)", m_x, m_y)
    } else {
      label_J <- "J (Jing rule)"
      label_W <- "W (Wasserman rule)"
    }

    choices <- setNames(
      c("J", "W", "none"),
      c(label_J, label_W, "Manual (use m_x, m_y)")
    )

    updateSelectInput(
      session, "bin_method",
      choices  = choices,
      selected = input$bin_method
    )
  })

  # Draw plot
  output$dpPlot <- renderPlot({
    res <- hist_res_reactive()
    res$plot$all
  })
}

shinyApp(ui = ui, server = server)
