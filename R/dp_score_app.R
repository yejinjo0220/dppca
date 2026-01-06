#' Run DP PCA score-plot Shiny app for user-supplied data
#'
#' @param X Numeric matrix/data.frame (n x p).
#' @param G Optional group vector (length n) or a column name in X.
#' @param title App title.
#' @export
dp_score_app <- function(X, G = NULL, title = "DP PCA Score Plot (User data)") {

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

  # ---- helpers (inside for clean shipping) ----

  `%||%` <- function(a, b) if (!is.null(a)) a else b

  # fraction(예: 1/8) 지원 + comma-separated
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
    # NULL이면 사용자 미입력 -> OK (내부 디폴트 분배 사용)
    if (is.null(x)) return(TRUE)
    is.numeric(x) && length(x) == need_len && all(is.finite(x))
  }

  # ---- UI ----
  ui <- shiny::fluidPage(
    shinyjs::useShinyjs(),

    shiny::titlePanel(title),

    shiny::sidebarLayout(
      shiny::sidebarPanel(
        shiny::checkboxInput("center", "Center variables (center = TRUE)", value = TRUE),
        shiny::checkboxInput("scale",  "Scale to unit variance (scale. = TRUE)", value = FALSE),
        shiny::checkboxInput("dp_pca_flag", "Use DP PCA", value = FALSE),

        shiny::tags$hr(),

        shiny::helpText("default is eps_total = 4, delta_total ≈ 1 / 10^{ceil(log10 n)}."),
        shiny::numericInput("eps",   "Total epsilon (eps_total):", value = 4,    min = 0.0001, step = 0.1),
        shiny::numericInput("delta", "Total delta (delta_total):", value = 1e-4, min = 0,      step = 1e-5),

        shiny::tags$hr(),

        shiny::uiOutput("eps_ratio_ui"),
        shiny::uiOutput("delta_ratio_ui"),

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

        shiny::checkboxInput("sampling", "Draw synthetic samples (sampling)", value = TRUE),

        shiny::tags$hr(),

        shiny::actionButton("run", "Run", class = "btn-primary", width = "100%"),

        shiny::tags$hr(),
        shiny::helpText("Tip: change settings, then click Run. Errors will be shown above the plot.")
      ),

      shiny::mainPanel(
        shiny::uiOutput("error_ui"),
        shinycssloaders::withSpinner(shiny::plotOutput("dpPlot", height = "700px"), type = 6)
      )
    )
  )

  # ---- server ----
  server <- function(input, output, session) {

    # (1) dp_pca_flag에 따라 ratio 입력 라벨 자동 변경
    output$eps_ratio_ui <- shiny::renderUI({
      need_len <- ratio_needed_len(input$dp_pca_flag)
      label <- if (need_len == 2L) {
        "eps_ratio (Quantile, Hist) [2 values]: ex. 0.6, 0.4 or 3, 2"
      } else {
        "eps_ratio (PCA, Quantile, Hist) [3 values]: ex. 0.3, 0.2, 0.5 or 3, 2, 5"
      }
      shiny::textInput("eps_ratio", label = label, value = input$eps_ratio %||% "")
    })

    output$delta_ratio_ui <- shiny::renderUI({
      need_len <- ratio_needed_len(input$dp_pca_flag)
      label <- if (need_len == 2L) {
        "delta_ratio (Quantile, Hist) [2 values]: ex. 0.6, 0.4"
      } else {
        "delta_ratio (PCA, Quantile, Hist) [3 values]: ex. 0.3, 0.2, 0.5"
      }
      shiny::textInput("delta_ratio", label = label, value = input$delta_ratio %||% "")
    })

    # (2) 앱 시작 시 delta 추천값 자동 세팅
    shiny::observeEvent(TRUE, {
      pow <- ceiling(log10(n))
      if (!is.finite(pow)) pow <- 2
      delta_default <- 10^(-pow)
      shiny::updateNumericInput(session, "delta", value = delta_default)
    }, once = TRUE)

    # (3) 결과 저장 state: ok/res/msg
    result_state <- shiny::reactiveVal(NULL)

    # (4) Run 버튼: 누르는 동안 비활성화 + 실행 + 다시 활성화
    shiny::observeEvent(input$run, {
      shinyjs::disable("run")
      on.exit(shinyjs::enable("run"), add = TRUE)

      out <- shiny::isolate({
        tryCatch({
          shiny::req(!is.na(input$eps), input$eps > 0)
          shiny::req(!is.na(input$delta), input$delta >= 0)

          need_len <- ratio_needed_len(input$dp_pca_flag)

          eps_ratio_raw   <- parse_ratio(input$eps_ratio)
          delta_ratio_raw <- parse_ratio(input$delta_ratio)

          # 길이 mismatch(예: "0"만 입력 중)이면 NULL로 보내서 내부 디폴트 분배 사용
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

          res <- if (is.null(group_vec)) {
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

          list(ok = TRUE, res = res, msg = NULL)

        }, error = function(e) {
          list(ok = FALSE, res = NULL, msg = conditionMessage(e))
        })
      })

      result_state(out)
    }, ignoreInit = TRUE)

    # (5) 에러 박스: 플롯 위에 빨간 경고 박스로 표시
    output$error_ui <- shiny::renderUI({
      x <- result_state()
      if (is.null(x) || isTRUE(x$ok)) return(NULL)

      shiny::div(
        class = "alert alert-danger",
        shiny::tags$b("Error: "),
        shiny::tags$span(x$msg)
      )
    })

    # (6) Plot: Run 전 안내 / 에러 안내 / 성공 플롯
    output$dpPlot <- shiny::renderPlot({
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

    # (7) bin_method 라벨 업데이트: 성공 결과 있을 때만
    shiny::observe({
      x <- result_state()
      if (is.null(x) || !isTRUE(x$ok)) return()

      res <- x$res
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
          c("J", "W", "none"),
          c(label_J, label_W, "Manual (use m_x, m_y)")
        ),
        selected = input$bin_method
      )
    })
  }

  shiny::runApp(shiny::shinyApp(ui, server), display.mode = "normal")
}
