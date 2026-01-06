library(shiny)
library(shinycssloaders)
library(shinyjs)
library(dppca)
#----------------------------
# helpers
#----------------------------

`%||%` <- function(a, b) if (!is.null(a)) a else b

parse_ratio <- function(txt) {
  txt <- gsub("\\s+", "", txt)
  if (txt == "") return(NULL)

  parts <- strsplit(txt, ",", fixed = TRUE)[[1]]
  parts <- parts[parts != ""]
  if (!length(parts)) return(NULL)

  parse_one <- function(s) {
    if (grepl("^[+-]?[0-9]*\\.?[0-9]+/[+-]?[0-9]*\\.?[0-9]+$", s)) {
      ab <- strsplit(s, "/", fixed = TRUE)[[1]]
      a <- suppressWarnings(as.numeric(ab[1]))
      b <- suppressWarnings(as.numeric(ab[2]))
      if (is.na(a) || is.na(b) || b == 0) return(NA_real_)
      return(a / b)
    }
    suppressWarnings(as.numeric(s))
  }

  vals <- vapply(parts, parse_one, numeric(1))
  if (any(is.na(vals))) return(NULL)
  vals
}

ratio_needed_len <- function(dp_pca_flag) if (isTRUE(dp_pca_flag)) 3L else 2L

ratio_ok_or_null <- function(x, need_len) {
  if (is.null(x)) return(TRUE)
  is.numeric(x) && length(x) == need_len && all(is.finite(x))
}

#----------------------------
# UI
#----------------------------
ui <- fluidPage(
  useShinyjs(),

  titlePanel("DP PCA Score Plot"),

  sidebarLayout(
    sidebarPanel(

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

      checkboxInput("center", "Center variables (center = TRUE)", TRUE),
      checkboxInput("scale",  "Scale to unit variance (scale. = TRUE)", FALSE),

      checkboxInput("dp_pca_flag", "Use DP PCA", FALSE),

      tags$hr(),

      helpText("default is eps_total = 4, delta_total ≈ 1 / 10^{ceil(log10 n)}."),
      numericInput("eps",   "Total epsilon (eps_total):", value = 4,    min = 0.0001, step = 0.1),
      numericInput("delta", "Total delta (delta_total):", value = 1e-4, min = 0,      step = 1e-5),

      tags$hr(),

      uiOutput("eps_ratio_ui"),
      uiOutput("delta_ratio_ui"),

      tags$hr(),

      selectInput(
        inputId = "bin_method",
        label   = "Bin method:",
        choices = c(
          "J (Jing rule)"               = "J",
          "W (Wasserman-Lei rule)"      = "W",
          "Manual (use m_x, m_y)"       = "none"
        ),
        selected = "J"
      ),

      conditionalPanel(
        condition = "input.bin_method == 'none'",
        numericInput("m_x", "Number of bins in x (m_x):", value = 10, min = 1, step = 1),
        numericInput("m_y", "Number of bins in y (m_y):", value = 10, min = 1, step = 1)
      ),

      tags$hr(),

      selectInput(
        inputId = "mechanism",
        label   = "Mechanism:",
        choices = c("all", "none", "add", "sparse"),
        selected = "all"
      ),
      checkboxInput("sampling", "Draw synthetic samples (sampling)", TRUE),

      tags$hr(),

      actionButton("run", "Run", class = "btn-primary", width = "100%"),

      tags$hr(),
      helpText("Tip: change settings, then click Run. Errors will be shown above the plot.")
    ),

    mainPanel(
      uiOutput("error_ui"),
      withSpinner(plotOutput("dpPlot", height = "700px"), type = 6)
    )
  )
)

#----------------------------
# SERVER
#----------------------------
server <- function(input, output, session) {

  # (1) dp_pca_flag에 따라 ratio 입력 라벨 자동 변경
  output$eps_ratio_ui <- renderUI({
    need_len <- ratio_needed_len(input$dp_pca_flag)
    label <- if (need_len == 2L) {
      "eps_ratio (Quantile, Hist) [2 values]: ex. 0.6, 0.4 or 3, 2"
    } else {
      "eps_ratio (PCA, Quantile, Hist) [3 values]: ex. 0.3, 0.2, 0.5 or 3, 2, 5"
    }
    textInput("eps_ratio", label = label, value = input$eps_ratio %||% "")
  })

  output$delta_ratio_ui <- renderUI({
    need_len <- ratio_needed_len(input$dp_pca_flag)
    label <- if (need_len == 2L) {
      "delta_ratio (Quantile, Hist) [2 values]: ex. 0.6, 0.4"
    } else {
      "delta_ratio (PCA, Quantile, Hist) [3 values]: ex. 0.3, 0.2, 0.5"
    }
    textInput("delta_ratio", label = label, value = input$delta_ratio %||% "")
  })

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

  # dataset 선택 시 eps_total, delta_total 추천값 설정
  observeEvent(dataset_reactive(), {
    dat <- dataset_reactive()
    X <- dat$X
    n <- nrow(X)
    if (is.null(n) || is.na(n) || n <= 0) return()

    eps_default <- 4
    pow <- ceiling(log10(n))
    if (!is.finite(pow)) pow <- 2
    delta_default <- 10^(-pow)

    updateNumericInput(session, "eps",   value = eps_default)
    updateNumericInput(session, "delta", value = delta_default)
  }, ignoreInit = FALSE)

  # (2) 결과 저장용 reactiveVal: ok/res/msg
  result_state <- reactiveVal(NULL)

  # (3) Run 버튼: 누르는 동안 비활성화 + 계산 + 다시 활성화
  observeEvent(input$run, {
    shinyjs::disable("run")
    on.exit(shinyjs::enable("run"), add = TRUE)

    out <- isolate({
      tryCatch({
        req(!is.na(input$eps), input$eps > 0)
        req(!is.na(input$delta), input$delta >= 0)

        dat <- dataset_reactive()
        X <- dat$X
        G <- dat$G

        need_len <- ratio_needed_len(input$dp_pca_flag)

        eps_ratio_raw   <- parse_ratio(input$eps_ratio)
        delta_ratio_raw <- parse_ratio(input$delta_ratio)

        # 길이 mismatch(예: "0"만 입력 중)이면 NULL로 보내서 디폴트 분배 사용
        eps_ratio   <- if (ratio_ok_or_null(eps_ratio_raw, need_len)) eps_ratio_raw else NULL
        delta_ratio <- if (ratio_ok_or_null(delta_ratio_raw, need_len)) delta_ratio_raw else NULL

        m_x <- if (input$bin_method == "none") input$m_x else NULL
        m_y <- if (input$bin_method == "none") input$m_y else NULL

        common_args <- list(
          center       = input$center,
          scale.       = input$scale,
          dp_pca_flag  = input$dp_pca_flag,
          cpp.option   = TRUE,
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

        res <- if (is.null(G)) {
          do.call(dp_score_plot, c(list(X = X, color = "#6A5ACD"), common_args))
        } else {
          do.call(dp_score_plot_group, c(list(X = X, G = G), common_args))
        }

        list(ok = TRUE, res = res, msg = NULL)

      }, error = function(e) {
        list(ok = FALSE, res = NULL, msg = conditionMessage(e))
      })
    })

    result_state(out)
  }, ignoreInit = TRUE)

  # (4) 에러 박스: 플롯 위에 빨간 경고 박스로 표시
  output$error_ui <- renderUI({
    x <- result_state()
    if (is.null(x) || isTRUE(x$ok)) return(NULL)

    div(
      class = "alert alert-danger",
      tags$b("Error: "),
      tags$span(x$msg)
    )
  })

  # (5) Plot 영역: Run 전 안내 / 성공 시 플롯 / 에러 시 (플롯은 비워두거나 안내)
  output$dpPlot <- renderPlot({
    x <- result_state()

    if (is.null(x)) {
      plot.new()
      text(0.5, 0.5, "Set options, then click 'Run'.")
      return(invisible(NULL))
    }

    if (!isTRUE(x$ok)) {
      plot.new()
      text(0.5, 0.5, "Cannot draw plot due to error.\nSee the red message above.")
      return(invisible(NULL))
    }

    x$res$plot$all
  })

  # (선택) bin_method label update: 성공 결과 있을 때만
  observe({
    x <- result_state()
    if (is.null(x) || !isTRUE(x$ok)) return()

    res <- x$res
    if (is.null(res$breaks) || is.null(res$breaks$x) || is.null(res$breaks$y)) return()

    m_x <- length(res$breaks$x) - 1
    m_y <- length(res$breaks$y) - 1

    if (input$bin_method == "J") {
      label_J <- sprintf("J (Jing rule / %d × %d bins)", m_x, m_y)
      label_W <- "W (Wasserman–Lei rule)"
    } else if (input$bin_method == "W") {
      label_J <- "J (Jing rule)"
      label_W <- sprintf("W (Wasserman–Lei rule / %d × %d bins)", m_x, m_y)
    } else {
      label_J <- "J (Jing rule)"
      label_W <- "W (Wasserman–Lei rule)"
    }

    choices <- setNames(
      c("J", "W", "none"),
      c(label_J, label_W, "Manual (use m_x, m_y)")
    )

    updateSelectInput(session, "bin_method",
                      choices  = choices,
                      selected = input$bin_method)
  })
}

shinyApp(ui = ui, server = server)

