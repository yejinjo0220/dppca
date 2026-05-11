# Shiny demo app for the dppca package
# Put this file at: inst/shiny/dppca-app/app.R

library(shiny)
library(dppca)

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------

.load_dppca_data <- function(name) {
  env <- new.env(parent = emptyenv())
  utils::data(list = name, package = "dppca", envir = env)
  if (!exists(name, envir = env, inherits = FALSE)) {
    stop("Could not load dataset '", name, "' from package dppca.")
  }
  get(name, envir = env, inherits = FALSE)
}

.default_delta <- function(n) {
  # Use a valid DP delta in (0, 1). This corresponds to 10^{-ceiling(log10(n))}.
  10^(-ceiling(log10(n)))
}

.group_vector <- function(G) {
  if (is.null(G)) return(NULL)
  if (is.data.frame(G)) {
    if (ncol(G) < 1L) return(NULL)
    return(G[[1L]])
  }
  if (is.matrix(G)) {
    if (ncol(G) < 1L) return(NULL)
    return(G[, 1L])
  }
  as.vector(G)
}

.infer_group_column <- function(dat) {
  if (!is.data.frame(dat)) dat <- as.data.frame(dat)
  if (ncol(dat) < 2L) return(NULL)

  nm <- names(dat)
  lower_nm <- tolower(nm)
  preferred <- match(c("color", "colour", "group", "groups", "label", "labels", "class", "cluster", "g"), lower_nm)
  preferred <- preferred[!is.na(preferred)]
  if (length(preferred) > 0L) return(nm[preferred[1L]])

  non_num <- which(!vapply(dat, is.numeric, logical(1)))
  if (length(non_num) > 0L) return(nm[non_num[1L]])

  low_card <- which(vapply(dat, function(x) {
    ux <- unique(stats::na.omit(x))
    is.numeric(x) && length(ux) >= 2L && length(ux) <= min(20L, floor(length(x) / 3)) && all(abs(ux - round(ux)) < 1e-8)
  }, logical(1)))
  if (length(low_card) > 0L && ncol(dat) >= 3L) return(nm[low_card[length(low_card)]])

  NULL
}

.split_grouped_data <- function(dat, group = NULL) {
  dat <- as.data.frame(dat)

  if (!is.null(group)) {
    if (is.character(group) && length(group) == 1L && group %in% names(dat)) {
      G <- dat[[group]]
      X <- dat[, setdiff(names(dat), group), drop = FALSE]
      return(list(X = X, G = G))
    }
    G <- .group_vector(group)
    return(list(X = dat, G = G))
  }

  group_col <- .infer_group_column(dat)
  if (!is.null(group_col)) {
    G <- dat[[group_col]]
    X <- dat[, setdiff(names(dat), group_col), drop = FALSE]
    return(list(X = X, G = G))
  }

  list(X = dat, G = NULL)
}

.load_grouped_example <- function(grouped_name, fallback_name, label, group_col = "color") {
  grouped <- .load_dppca_data(grouped_name)

  # For dppca example datasets, gau_g and eur_map_g are the corresponding
  # feature data with an additional group/color column. Use that column as
  # group labels, and remove it from the PCA feature matrix.
  if (is.data.frame(grouped) || is.matrix(grouped)) {
    grouped_df <- as.data.frame(grouped)
    if (!is.null(group_col) && group_col %in% names(grouped_df)) {
      split <- .split_grouped_data(grouped_df, group_col)
    } else {
      split <- .split_grouped_data(grouped_df)
    }
    return(list(X = split$X, G = split$G, label = label))
  }

  # Fallback behavior: if the grouped object is only a label vector,
  # pair it with the corresponding feature matrix.
  fallback_X <- .load_dppca_data(fallback_name)
  list(X = fallback_X, G = .group_vector(grouped), label = label)
}

.load_example_dataset <- function(name) {
  switch(
    name,
    gau = list(X = .load_dppca_data("gau"), G = NULL, label = "Gaussian"),
    gau_grouped = .load_grouped_example("gau_g", "gau", "Gaussian with groups", group_col = "color"),
    eur_map = list(X = .load_dppca_data("eur_map"), G = NULL, label = "Europe map"),
    eur_map_grouped = .load_grouped_example("eur_map_g", "eur_map", "Europe map with groups", group_col = "color"),
    adult = list(X = .load_dppca_data("adult"), G = NULL, label = "Adult"),
    stop("Unknown example dataset: ", name)
  )
}

.to_numeric_matrix <- function(dat) {
  dat <- as.data.frame(dat)

  # Keep only numeric columns. This makes the app robust for example data
  # that may include labels or other non-numeric columns.
  is_num <- vapply(dat, is.numeric, logical(1))
  dat_num <- dat[, is_num, drop = FALSE]

  if (ncol(dat_num) < 2L) {
    stop("The selected data must contain at least two numeric columns.")
  }
  if (nrow(dat_num) < 2L) {
    stop("The selected data must contain at least two rows.")
  }

  X <- as.matrix(dat_num)
  storage.mode(X) <- "double"
  X
}

.selected_scree_methods <- function(input) {
  method <- input$scree_method
  if (is.null(method) || !length(method)) {
    return(character(0))
  }
  unique(method)
}

.make_scree_control <- function(input) {
  method <- .selected_scree_methods(input)

  clipped <- clipped_control(
    C_clip = input$C_clip
  )

  pmwm <- pmwm_control(
    beta = input$pmwm_beta,
    a = input$pmwm_a,
    b = input$pmwm_b,
    trim_const = input$pmwm_trim_const,
    eta = input$pmwm_eta,
    split_mode = input$pmwm_split_mode
  )

  huber <- huber_control(
    mu0 = input$huber_mu0,
    eta0 = input$huber_eta0,
    T = input$huber_T,
    M = input$huber_M,
    k_min_m2 = input$huber_k_min_m2,
    k_max_m2 = input$huber_k_max_m2,
    m2_frac = input$huber_m2_frac
  )

  controls <- list(clipped = clipped, pmwm = pmwm, huber = huber)

  if (length(method) == 1L) {
    controls[[method]]
  } else {
    controls[method]
  }
}

# ------------------------------------------------------------
# Data supplied by dppca_app(X, group = ...)
# ------------------------------------------------------------

.supplied_X <- getOption("dppca.app.X", NULL)
.supplied_group <- getOption("dppca.app.group", NULL)
.has_supplied_data <- !is.null(.supplied_X)

# ------------------------------------------------------------
# UI
# ------------------------------------------------------------

ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      body {
        background: #f7f8fb;
      }
      .app-title {
        margin-bottom: 14px;
      }
      .well {
        background: #ffffff;
        border: 1px solid #e5e7eb;
        border-radius: 10px;
        box-shadow: 0 1px 3px rgba(15, 23, 42, 0.06);
      }
      .main-plot-area {
        max-width: 1160px;
        margin: 0 auto;
      }
      .plot-card {
        background: #ffffff;
        border: 1px solid #e5e7eb;
        border-radius: 12px;
        box-shadow: 0 2px 8px rgba(15, 23, 42, 0.06);
        padding: 16px 18px 18px 18px;
        margin-bottom: 24px;
      }
      .plot-card-scree {
        max-width: 980px;
        margin-left: auto;
        margin-right: auto;
      }
      .plot-card-score {
        max-width: 1120px;
        margin-left: auto;
        margin-right: auto;
      }
      .plot-card h4 {
        margin-top: 0;
        margin-bottom: 12px;
        font-weight: 600;
        color: #1f2937;
      }
      .plot-card-note {
        margin-top: -4px;
        margin-bottom: 10px;
        color: #6b7280;
        font-size: 12px;
      }
      details.control-details {
        margin-top: 12px;
        margin-bottom: 10px;
        border: 1px solid #ddd;
        border-radius: 4px;
        background: #ffffff;
      }
      details.control-details summary {
        cursor: pointer;
        font-weight: 600;
        padding: 8px 10px;
        list-style: none;
        user-select: none;
      }
      details.control-details summary::-webkit-details-marker {
        display: none;
      }
      details.control-details summary::before {
        content: '▶ ';
        display: inline-block;
        width: 1.2em;
      }
      details.control-details[open] summary::before {
        content: '▼ ';
      }
      details.control-details .control-body {
        padding: 0 10px 10px 10px;
      }
    "))
  ),

  div(class = "app-title", h2("dppca: Differentially Private PCA Visualization")),

  sidebarLayout(
    sidebarPanel(
      width = 4,

      h4("Data"),
      radioButtons(
        "data_source", "Data source",
        choices = if (.has_supplied_data) {
          c("Provided data" = "provided", "Example data" = "example", "Upload CSV" = "upload")
        } else {
          c("Example data" = "example", "Upload CSV" = "upload")
        },
        selected = if (.has_supplied_data) "provided" else "example"
      ),
      conditionalPanel(
        condition = "input.data_source == 'provided'",
        helpText("Using the data supplied through dppca_app(X, group = ...).")
      ),
      conditionalPanel(
        condition = "input.data_source == 'example'",
        selectInput(
          "example_dataset", "Example dataset",
          choices = c(
            "Gaussian" = "gau",
            "Gaussian with groups" = "gau_grouped",
            "Europe map" = "eur_map",
            "Europe map with groups" = "eur_map_grouped",
            "Adult" = "adult"
          ),
          selected = "gau_grouped"
        )
      ),
      conditionalPanel(
        condition = "input.data_source == 'upload'",
        fileInput("csv_file", "Upload CSV file", accept = c(".csv")),
        checkboxInput("csv_header", "Header", value = TRUE)
      ),

      hr(),
      h4("Common PCA / privacy options"),
      numericInput("k", "Number of leading PCs for scree plot (k)", value = 4, min = 1, step = 1),
      numericInput("eps_total", "Total epsilon", value = 5, min = 1e-6, step = 0.1),
      numericInput("delta_total", "Total delta", value = 1e-3, min = 1e-12, max = 0.99, step = 1e-6),
      checkboxInput(
        "split_budget_between_plots",
        "Split total privacy budget between scree and score plots",
        value = FALSE
      ),
      checkboxInput("center", "Center columns", value = TRUE),
      checkboxInput("standardize", "Standardize columns", value = FALSE),
      checkboxInput("g_dppca", "Use private PC directions", value = FALSE),

      hr(),
      tabsetPanel(
        id = "settings_tab",

        tabPanel(
          "Scree controls",
          checkboxGroupInput(
            "scree_method", "DP scree method",
            choices = c("Clipped mean" = "clipped", "PMWM" = "pmwm", "Huber" = "huber"),
            selected = c("clipped", "pmwm", "huber"),
            inline = TRUE
          ),
          selectInput(
            "scree_type", "Plot type",
            choices = c("PVE" = "pve", "Raw scree" = "scree"),
            selected = "pve"
          ),
          checkboxInput("scree_mono", "Apply monotone post-processing", value = TRUE),
          helpText("Select one or more DP scree methods to display together."),

          tags$details(
            class = "control-details",
            tags$summary("Clipped control"),
            tags$div(
              class = "control-body",
              numericInput("C_clip", "C_clip", value = 3, min = 1e-6, step = 0.5)
            )
          ),

          tags$details(
            class = "control-details",
            tags$summary("PMWM control"),
            tags$div(
              class = "control-body",
              numericInput("pmwm_beta", "beta", value = 1.01, min = 1.000001, step = 0.01),
              numericInput("pmwm_a", "a", value = 0, step = 1),
              numericInput("pmwm_b", "b", value = 10, step = 1),
              numericInput("pmwm_trim_const", "trim_const", value = 10, min = 0, step = 1),
              numericInput("pmwm_eta", "eta", value = 0.01, min = 0, max = 0.49, step = 0.01),
              checkboxInput("pmwm_split_mode", "split_mode", value = TRUE)
            )
          ),

          tags$details(
            class = "control-details",
            tags$summary("Huber control"),
            tags$div(
              class = "control-body",
              numericInput("huber_mu0", "mu0", value = 0, step = 0.5),
              numericInput("huber_eta0", "eta0", value = 1, min = 1e-8, step = 0.1),
              numericInput("huber_T", "T", value = 50, min = 1, step = 1),
              numericInput("huber_M", "M", value = 20, min = 1, step = 1),
              numericInput("huber_k_min_m2", "k_min_m2", value = -40, step = 1),
              numericInput("huber_k_max_m2", "k_max_m2", value = 40, step = 1),
              numericInput("huber_m2_frac", "m2_frac", value = 0.25, min = 1e-6, max = 0.99, step = 0.05)
            )
          ),

          actionButton("run_scree", "Run DP scree plot", class = "btn-primary")
        ),

        tabPanel(
          "Score controls",
          numericInput("axis_x", "x-axis PC", value = 1, min = 1, step = 1),
          numericInput("axis_y", "y-axis PC", value = 2, min = 1, step = 1),
          checkboxGroupInput(
            "score_method", "DP score method",
            choices = c("Additive histogram" = "add", "Sparse histogram" = "sparse"),
            selected = c("add", "sparse"),
            inline = TRUE
          ),
          selectInput(
            "bin_method", "Bin recommendation method",
            choices = c("WZ", "Lei", "none"),
            selected = "WZ"
          ),
          conditionalPanel(
            condition = "input.bin_method == 'none'",
            numericInput("m_x", "m_x", value = 20, min = 1, step = 1),
            numericInput("m_y", "m_y", value = 20, min = 1, step = 1)
          ),
          checkboxInput("use_group", "Use group labels when available", value = TRUE),
          conditionalPanel(
            condition = "input.data_source == 'upload' && input.use_group",
            textInput(
              "upload_group_col",
              "group:",
              value = "",
              placeholder = "Enter group column name, e.g. color"
            ),
            helpText("For uploaded CSV files, enter the column name containing group labels. This column is excluded from PCA and used only for grouped score plots.")
          ),
          helpText("Score plots use a fixed frame inflation of 0.20 and fixed color #6A5ACD for ungrouped plots."),
          actionButton("run_score", "Run DP score plot", class = "btn-primary")
        )
      )
    ),

    mainPanel(
      width = 8,
      div(
        class = "main-plot-area",
        div(
          class = "plot-card plot-card-scree",
          h4("DP Scree Plot"),
          div(class = "plot-card-note", "PVE curves are shown in a taller, centered panel to avoid a flattened appearance."),
          plotOutput("scree_plot", height = "520px")
        ),
        div(
          class = "plot-card plot-card-score",
          h4("DP Score Plot"),
          plotOutput("score_plot", height = "520px")
        )
      )
    )
  )
)

# ------------------------------------------------------------
# Server
# ------------------------------------------------------------

server <- function(input, output, session) {

  selected_data <- reactive({
    if (input$data_source == "provided") {
      req(.supplied_X)

      # For supplied data, keep the original data for dp_score_plot_group()
      # when G is a column name such as "color".  Also prepare a numeric-only
      # feature matrix for scree/non-group score calculations.
      split <- .split_grouped_data(.supplied_X, .supplied_group)

      G_arg <- NULL
      if (!is.null(.supplied_group)) {
        G_arg <- .supplied_group
      } else {
        inferred <- .infer_group_column(.supplied_X)
        if (!is.null(inferred)) G_arg <- inferred
      }

      list(
        X_feature = split$X,
        X_score = .supplied_X,
        G_arg = G_arg,
        label = "Provided data"
      )

    } else if (input$data_source == "example") {
      ex <- .load_example_dataset(input$example_dataset)

      # For the grouped example datasets, gau_g and eur_map_g already contain
      # the group column named "color".  Pass the full data frame and G = "color"
      # to dp_score_plot_group(), matching the console usage:
      # dp_score_plot_group(X = gau_g, G = "color", ...).
      if (input$example_dataset == "gau_grouped") {
        X_raw <- .load_dppca_data("gau_g")
        return(list(
          X_feature = .split_grouped_data(X_raw, "color")$X,
          X_score = X_raw,
          G_arg = "color",
          label = "Gaussian with groups"
        ))
      }

      if (input$example_dataset == "eur_map_grouped") {
        X_raw <- .load_dppca_data("eur_map_g")
        return(list(
          X_feature = .split_grouped_data(X_raw, "color")$X,
          X_score = X_raw,
          G_arg = "color",
          label = "Europe map with groups"
        ))
      }

      list(
        X_feature = ex$X,
        X_score = ex$X,
        G_arg = NULL,
        label = ex$label
      )

    } else {
      req(input$csv_file)
      dat <- utils::read.csv(input$csv_file$datapath, header = input$csv_header, check.names = FALSE)

      # For uploaded CSV files, do not force automatic grouping.
      # If the user checks "Use group labels when available" and types a
      # group column name, pass the full uploaded data to dp_score_plot_group()
      # with G = <that column name>.  The feature matrix used for scree and
      # non-group score plots excludes that group column.
      group_col <- input$upload_group_col
      if (is.null(group_col)) group_col <- ""
      group_col <- trimws(group_col)

      if (isTRUE(input$use_group) && nzchar(group_col)) {
        validate(
          need(group_col %in% names(dat),
               paste0("The uploaded CSV does not contain a column named '", group_col, "'."))
        )

        split <- .split_grouped_data(dat, group_col)
        list(
          X_feature = split$X,
          X_score = dat,
          G_arg = group_col,
          label = "Uploaded CSV"
        )
      } else {
        list(
          X_feature = dat,
          X_score = dat,
          G_arg = NULL,
          label = "Uploaded CSV"
        )
      }
    }
  })

  X_data <- reactive({
    .to_numeric_matrix(selected_data()$X_feature)
  })

  plot_budget <- reactive({
    validate(
      need(input$eps_total > 0, "epsilon must be positive."),
      need(input$delta_total > 0 && input$delta_total < 1, "delta must be in (0, 1).")
    )

    if (isTRUE(input$split_budget_between_plots)) {
      list(
        scree_eps = input$eps_total / 2,
        scree_delta = input$delta_total / 2,
        score_eps = input$eps_total / 2,
        score_delta = input$delta_total / 2
      )
    } else {
      list(
        scree_eps = input$eps_total,
        scree_delta = input$delta_total,
        score_eps = input$eps_total,
        score_delta = input$delta_total
      )
    }
  })

  observeEvent(X_data(), {
    X <- X_data()
    n <- nrow(X)
    max_pc <- max(2L, min(n - 1L, ncol(X)))
    updateNumericInput(session, "k", max = max_pc, value = min(input$k, max_pc))
    updateNumericInput(session, "axis_x", max = max_pc, value = min(input$axis_x, max_pc))
    updateNumericInput(session, "axis_y", max = max_pc, value = min(max(input$axis_y, 2L), max_pc))
    updateNumericInput(session, "delta_total", value = .default_delta(n))
  })

  scree_event <- eventReactive(input$run_scree, {
    X <- X_data()
    max_pc <- min(nrow(X) - 1L, ncol(X))
    validate(
      need(input$k >= 1 && input$k <= max_pc, paste0("k must be between 1 and ", max_pc, ".")),
      need(length(.selected_scree_methods(input)) >= 1L, "Choose at least one DP scree method."),
      need(input$eps_total > 0, "epsilon must be positive."),
      need(input$delta_total > 0 && input$delta_total < 1, "delta must be in (0, 1)."),
      need(input$pmwm_a < input$pmwm_b, "For PMWM, a must be smaller than b."),
      need(input$huber_k_min_m2 < input$huber_k_max_m2, "For Huber, k_min_m2 must be smaller than k_max_m2.")
    )

    list(
      X = X,
      control = .make_scree_control(input),
      budget = plot_budget()
    )
  })

  output$scree_plot <- renderPlot({
    obj <- scree_event()
    tryCatch(
      {
        dp_scree_plot(
          X = obj$X,
          k = input$k,
          method = .selected_scree_methods(input),
          eps_total = obj$budget$scree_eps,
          delta_total = obj$budget$scree_delta,
          center = input$center,
          standardize = input$standardize,
          control = obj$control,
          g_dppca = input$g_dppca,
          cpp.option = FALSE,
          mono = input$scree_mono,
          type = input$scree_type
        )
      },
      error = function(e) {
        showNotification(paste("Scree plot error:", conditionMessage(e)), type = "error", duration = 10)
        plot.new()
        text(0.5, 0.5, paste("Scree plot error:\n", conditionMessage(e)), cex = 0.9)
      }
    )
  })

  score_event <- eventReactive(input$run_score, {
    X <- X_data()
    max_pc <- min(nrow(X) - 1L, ncol(X))
    axes <- c(as.integer(input$axis_x), as.integer(input$axis_y))
    validate(
      need(length(unique(axes)) == 2L, "Choose two different PCs for the score plot axes."),
      need(all(axes >= 1 & axes <= max_pc), paste0("Axes must be between 1 and ", max_pc, ".")),
      need(length(input$score_method) >= 1L, "Choose at least one DP score method."),
      need(input$eps_total > 0, "epsilon must be positive."),
      need(input$delta_total > 0 && input$delta_total < 1, "delta must be in (0, 1).")
    )

    m_x <- if (input$bin_method == "none") input$m_x else NULL
    m_y <- if (input$bin_method == "none") input$m_y else NULL

    dat <- selected_data()
    G_arg <- dat$G_arg
    budget <- plot_budget()

    tryCatch(
      {
        if (isTRUE(input$use_group) && !is.null(G_arg)) {
          dp_score_plot_group(
            X = dat$X_score,
            G = G_arg,
            center = input$center,
            standardize = input$standardize,
            g_dppca = input$g_dppca,
            cpp.option = FALSE,
            axes = axes,
            eps_total = budget$score_eps,
            delta_total = budget$score_delta,
            method = input$score_method,
            m_x = m_x,
            m_y = m_y,
            bin_method = input$bin_method
          )
        } else {
          dp_score_plot(
            X = X,
            center = input$center,
            standardize = input$standardize,
            g_dppca = input$g_dppca,
            cpp.option = FALSE,
            axes = axes,
            eps_total = budget$score_eps,
            delta_total = budget$score_delta,
            method = input$score_method,
            m_x = m_x,
            m_y = m_y,
            bin_method = input$bin_method
          )
        }
      },
      error = function(e) {
        showNotification(paste("Score plot error:", conditionMessage(e)), type = "error", duration = 10)
        structure(list(error = conditionMessage(e)), class = "dppca_app_error")
      }
    )
  })

  output$score_plot <- renderPlot({
    res <- score_event()

    if (inherits(res, "dppca_app_error")) {
      plot.new()
      text(0.5, 0.5, paste("Score plot error:\n", res$error), cex = 0.9)
    } else if (!is.null(res$plot$all)) {
      print(res$plot$all)
    } else if (!is.null(res$plot)) {
      print(res$plot)
    } else {
      print(res)
    }
  })
}

shinyApp(ui, server)

