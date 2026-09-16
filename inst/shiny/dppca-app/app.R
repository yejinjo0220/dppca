library(shiny)
library(dppca)

# Data and small application helpers -----------------------------------------

.supplied_X <- getOption("dppca.app.X", NULL)
.supplied_group <- getOption("dppca.app.group", NULL)
.has_supplied_data <- !is.null(.supplied_X)

.app_fun <- function(name) {
  fun <- get0(name, envir = asNamespace("dppca"), mode = "function")
  if (is.null(fun)) {
    stop(
      "The loaded dppca package does not contain ", name,
      ". Update the package R files, then run devtools::document() ",
      "and devtools::load_all() before opening the app.",
      call. = FALSE
    )
  }
  fun
}

.app_load_data <- function(name) {
  env <- new.env(parent = emptyenv())
  utils::data(list = name, package = "dppca", envir = env)
  if (!exists(name, envir = env, inherits = FALSE)) {
    stop("Could not load example dataset '", name, "'.", call. = FALSE)
  }
  get(name, envir = env, inherits = FALSE)
}

.app_data <- function(data, label, group = NULL) {
  data <- as.data.frame(data)
  preferred <- intersect(c("group", "groups", "label", "class"), names(data))
  group_column <- if (length(preferred)) preferred[1] else ""
  external_group <- NULL

  if (!is.null(group)) {
    if (is.character(group) && length(group) == 1L && group %in% names(data)) {
      group_column <- group
    } else {
      if (is.data.frame(group)) group <- group[[1]]
      if (is.matrix(group)) group <- group[, 1]
      if (length(group) != nrow(data)) {
        stop("Group labels must have one value per observation.", call. = FALSE)
      }
      external_group <- group
      group_column <- ".provided_group"
    }
  }

  list(
    data = data,
    label = label,
    group_column = group_column,
    external_group = external_group
  )
}

.app_matrix <- function(data) {
  numeric_columns <- vapply(data, is.numeric, logical(1))
  X <- as.matrix(data[, numeric_columns, drop = FALSE])
  storage.mode(X) <- "double"

  if (nrow(X) < 2L || ncol(X) < 2L) {
    stop("Choose data with at least two rows and two numeric variables.",
         call. = FALSE)
  }
  if (any(!is.finite(X))) {
    stop("Numeric variables must not contain missing or infinite values.",
         call. = FALSE)
  }
  if (anyDuplicated(colnames(X))) {
    stop("Variable names must be unique.", call. = FALSE)
  }
  X
}

.app_checked_data <- function(expr) {
  tryCatch(force(expr), error = function(e) {
    if (inherits(e, "shiny.silent.error")) stop(e)
    validate(need(FALSE, conditionMessage(e)))
  })
}

.app_integer <- function(value, label, lower = 1, upper = Inf) {
  if (
    length(value) != 1L || !is.numeric(value) || !is.finite(value) ||
      value != floor(value) || value < lower || value > upper
  ) {
    stop(label, " must be an integer between ", lower, " and ", upper,
         ".", call. = FALSE)
  }
  as.integer(value)
}

.app_controls <- function(input, methods) {
  if (!length(methods)) {
    stop("Choose at least one scree method.", call. = FALSE)
  }

  controls <- lapply(methods, function(method) {
    switch(
      method,
      clipped = clipped_control(C_clip = input$C_clip),
      pmwm = pmwm_control(
        beta = input$pmwm_beta,
        a = input$pmwm_a,
        b = input$pmwm_b,
        trim_const = input$pmwm_trim_const,
        eta = input$pmwm_eta,
        split_mode = input$pmwm_split_mode
      ),
      huber = huber_control(
        mu0 = input$huber_mu0,
        eta0 = input$huber_eta0,
        T = .app_integer(input$huber_T, "Huber T"),
        M = .app_integer(input$huber_M, "Huber M"),
        k_min_m2 = input$huber_k_min_m2,
        k_max_m2 = input$huber_k_max_m2,
        m2_frac = input$huber_m2_frac
      )
    )
  })
  names(controls) <- methods
  if (length(methods) == 1L) controls[[1]] else controls
}

.app_pca_split <- function(points) {
  if (
    length(points) != 2L || any(!is.finite(points)) ||
      points[1] <= 0 || points[2] >= 300 || points[1] >= points[2]
  ) {
    stop("Keep a positive share for Loading, Scree, and Score.", call. = FALSE)
  }

  # Two boundaries divide one budget into three parts. Using 300 units keeps
  # the default allocation at exactly one third for each stage.
  c(
    loading = points[1],
    scree = points[2] - points[1],
    score = 300 - points[2]
  ) / 300
}

.app_privacy <- function(
    input,
    mode,
    budget_mode = input$budget_mode
) {
  if (identical(budget_mode, "direct")) {
    budgets <- list(
      loading = c(eps = input$budget_loading_eps, delta = input$budget_loading_delta)
    )
    if (mode %in% c("pca", "scree")) {
      budgets$scree <- c(eps = input$budget_scree_eps, delta = input$budget_scree_delta)
    }
    if (mode %in% c("pca", "score")) {
      budgets$score <- c(eps = input$budget_score_eps, delta = input$budget_score_delta)
    }
    return(do.call(privacy_control, budgets))
  }

  split <- switch(
    mode,
    pca = .app_pca_split(input$pca_split_points),
    scree = c(loading = input$loading_share, scree = 1 - input$loading_share),
    score = c(loading = input$loading_share, score = 1 - input$loading_share)
  )
  privacy_control(eps = input$eps, delta = input$delta, split = split)
}

.app_budget <- function(privacy, mode, methods = 1L) {
  b <- privacy$components
  rows <- list(Loading = b$loading)
  if (!is.null(b$scree)) rows$Scree <- b$scree
  if (!is.null(b$score)) {
    rows[["Score frame"]] <- c(eps = 0.35 * b$score[["eps"]], delta = 0)
    rows[["Score histogram"]] <- c(
      eps = 0.65 * b$score[["eps"]], delta = b$score[["delta"]]
    )
  }
  rows[["Total"]] <- privacy$total

  if (methods > 1L && mode == "scree") {
    rows[["Total (selected methods)"]] <- b$loading + methods * b$scree
  }
  if (methods > 1L && mode == "score") {
    rows[["Total (selected methods)"]] <- b$loading +
      rows[["Score frame"]] + methods * rows[["Score histogram"]]
  }
  .app_budget_table(do.call(rbind, rows))
}

.app_budget_table <- function(values) {
  data.frame(
    Stage = rownames(values),
    eps = format(signif(values[, "eps"], 4), trim = TRUE),
    delta = format(values[, "delta"], scientific = TRUE, digits = 3),
    row.names = NULL,
    check.names = FALSE
  )
}

.app_capture_plot <- function(expr) {
  previous_device <- grDevices::dev.cur()
  grDevices::pdf(NULL)
  device <- grDevices::dev.cur()
  on.exit({
    if (device %in% grDevices::dev.list()) {
      grDevices::dev.off(device)
    }
    # Closing the PDF selects the next open device, which may be RStudio.
    # Explicitly return drawing to the original Shiny graphics device.
    if (previous_device %in% grDevices::dev.list()) {
      grDevices::dev.set(previous_device)
    }
  }, add = TRUE)
  force(expr)
}

.app_result_note <- function(result) {
  if (is.null(result)) {
    return(tags$p(class = "result-note", "Choose settings, then click Run analysis."))
  }
  tags$p(
    class = "result-note",
    paste0(
      result$label, " · ", result$n, " observations · ", result$p,
      " variables · ", result$time,
      ". Estimation settings take effect on the next run."
    )
  )
}

# Plot controls only redraw stored estimates. A run keeps one sampling seed.
.app_with_seed <- function(
    seed,
    expr
) {
  if (
    !is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      seed < 0 || seed > .Machine$integer.max || seed != floor(seed)
  ) {
    stop("Invalid plot seed.", call. = FALSE)
  }

  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) {
    previous_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }

  on.exit({
    if (had_seed) {
      assign(".Random.seed", previous_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)

  set.seed(seed)
  force(expr)
}

.app_pca_plots <- function(
    fit,
    seed,
    alpha = 0,
    variables = 10,
    labels = TRUE,
    loading_display = "point",
    score_display = "histogram",
    variable_display = "vector",
    circle = TRUE,
    base_size = 12
) {
  .app_with_seed(seed, .app_capture_plot(
    .app_fun("plot.dp_pca")(
      x = fit,
      type = "all",
      non_private = TRUE,
      loading = list(display = loading_display),
      score = list(display = score_display),
      variable = list(display = variable_display, circle = circle),
      biplot = list(alpha = alpha),
      variables = variables,
      labels = labels,
      base_size = base_size
    )
  ))
}

.app_biplot_plots <- function(
    fit,
    seed,
    alpha = 0,
    variables = 10,
    labels = TRUE,
    base_size = 12
) {
  # Both plot types sample private points before constructing any panels.
  # Default sampling uses the stored histogram and leaves reference points intact.
  .app_with_seed(seed, .app_capture_plot(
    .app_fun("plot.dp_pca")(
      x = fit,
      type = "biplot",
      non_private = TRUE,
      biplot = list(alpha = alpha),
      variables = variables,
      labels = labels,
      base_size = base_size
    )
  ))
}

.app_draw_pca <- function(
    result,
    view = "all",
    version = "private",
    compare = TRUE
) {
  types <- c("scree", "score", "loading", "variable", "biplot")
  view <- match.arg(view, c("all", types))
  version <- match.arg(version, c("private", "nonprivate"))

  if (view == "all") {
    plots <- result[[version]]$plots[types]
    if (length(plots) != length(types) || any(vapply(plots, is.null, logical(1)))) {
      stop("The stored PCA result does not contain all five plots.", call. = FALSE)
    }
    title <- if (version == "private") "Differentially private PCA" else "Non-private PCA"
    output <- patchwork::wrap_plots(unname(plots), ncol = 3) +
      patchwork::plot_annotation(title = title)
  } else {
    versions <- if (isTRUE(compare)) c("nonprivate", "private") else version
    plots <- lapply(versions, function(nm) result[[nm]]$plots[[view]])
    if (any(vapply(plots, is.null, logical(1)))) {
      stop("The selected PCA plot is not available in the stored result.", call. = FALSE)
    }
    output <- patchwork::wrap_plots(plots, ncol = length(plots))
  }

  print(output)
  invisible(output)
}

.app_draw_scree <- function(
    result,
    type = "pve"
) {
  type <- match.arg(type, c("pve", "scree"))
  series <- intersect(c("nonprivate", "clipped", "pmwm", "huber"), names(result))
  values <- lapply(series, function(nm) result[[nm]][[type]])
  valid <- vapply(values, function(x) {
    is.numeric(x) && length(x) > 0L && all(is.finite(x))
  }, logical(1))

  if (!length(values) || !all(valid) || length(unique(lengths(values))) != 1L) {
    stop("The stored scree result does not contain valid plot values.", call. = FALSE)
  }

  control <- .app_fun(".merge_scree_plot_control")(NULL, type)
  k <- length(values[[1L]])
  index <- seq_len(k)
  ylim <- range(unlist(values, use.names = FALSE))
  pad <- if (diff(ylim) == 0) max(0.5, abs(ylim[1L]) * 0.05) else diff(ylim) * 0.04
  ylim <- ylim + c(-pad, pad)
  if (type == "pve") {
    ylim[1L] <- max(0, ylim[1L])
  }

  # Keep the plotting panel square without changing either axis scale.
  previous_par <- graphics::par(pty = "s")
  on.exit(graphics::par(previous_par), add = TRUE)

  graphics::matplot(
    index,
    do.call(cbind, values),
    type = "b",
    col = unname(control$col[series]),
    lty = unname(control$lty[series]),
    pch = unname(control$pch[series]),
    lwd = control$lwd,
    cex = control$point_cex,
    xlab = control$xlab,
    ylab = control$ylab,
    main = control$title,
    xlim = if (k == 1L) c(0.5, 1.5) else c(1, k),
    ylim = ylim,
    xaxt = "n",
    cex.main = control$cex_main,
    cex.lab = control$cex_lab,
    cex.axis = control$cex_axis
  )
  graphics::axis(1, at = index, labels = index, cex.axis = control$cex_axis)
  graphics::legend(
    control$legend_position,
    legend = unname(control$legend_labels[series]),
    col = unname(control$col[series]),
    lty = unname(control$lty[series]),
    pch = unname(control$pch[series]),
    lwd = control$lwd,
    pt.cex = control$point_cex,
    bty = control$legend_bty,
    cex = control$cex_legend
  )

  invisible(result)
}

.app_draw_score <- function(
    result,
    display = "histogram"
) {
  choices <- c("histogram", "sample")
  if (
    !is.character(display) || !length(display) || anyNA(display) ||
      any(!display %in% choices)
  ) {
    stop("Choose at least one score display: histogram or sample.", call. = FALSE)
  }
  display <- intersect(choices, display)
  stored <- result$plot
  methods <- intersect(c("add", "sparse"), names(stored))
  plots <- list()

  if ("histogram" %in% display) {
    plots <- c(plots, list(stored$nonprivate), unname(stored[methods]))
  }
  if ("sample" %in% display) {
    private <- lapply(methods, function(method) stored$sample[[method]]$uniform)
    plots <- c(plots, list(stored$scatter), private)
  }

  if (!length(methods) || any(vapply(plots, is.null, logical(1)))) {
    stop("The selected score plots are not available in the stored result.", call. = FALSE)
  }

  output <- patchwork::wrap_plots(
    plots,
    ncol = 1L + length(methods),
    nrow = length(display),
    byrow = TRUE,
    guides = "collect"
  ) & ggplot2::theme(legend.position = "right")
  print(output)
  invisible(output)
}

.app_draw_loading <- function(
    result,
    display = "vector"
) {
  choices <- c("vector", "point", "heatmap")
  if (
    !is.character(display) || !length(display) || anyNA(display) ||
      any(!display %in% choices)
  ) {
    stop("Choose at least one loading display: vector, point, or heatmap.",
         call. = FALSE)
  }
  display <- intersect(choices, display)
  stored <- result$plot
  plots <- c(
    unname(stored$nonprivate[display]),
    unname(stored$private[display])
  )

  if (
    length(plots) != 2L * length(display) ||
      any(vapply(plots, is.null, logical(1)))
  ) {
    stop("The selected loading plots are not available in the stored result.",
         call. = FALSE)
  }

  single <- length(display) == 1L
  output <- patchwork::wrap_plots(
    plots,
    ncol = if (single) 2L else length(display),
    nrow = if (single) 1L else 2L,
    byrow = TRUE,
    guides = "collect"
  ) &
    ggplot2::theme(legend.position = "right")
  print(output)
  invisible(output)
}

# A single budget bar, with one movable divider on each side of Scree.
.app_pca_split_control <- function() {
  div(
    class = "pca-split",
    title = "Drag either divider to change the three budget shares.",
    tags$style(HTML("
      .pca-split .shiny-input-container { width: 100%; margin-bottom: 0; }
      .pca-split .irs { height: 36px; }
      .pca-split .irs-line { top: 14px; height: 10px; border: 0;
        background: linear-gradient(to right,
          #187b79 0%, #187b79 33.333333%,
          #6688b0 33.333333%, #6688b0 66.666667%,
          #bc8642 66.666667%, #bc8642 100%); }
      .pca-split .irs-bar { top: 14px; height: 10px; background: transparent; border: 0; }
      .pca-split .irs-handle { top: 8px; }
      .pca-split .irs-from, .pca-split .irs-to, .pca-split .irs-single,
      .pca-split .irs-min, .pca-split .irs-max { display: none; }
      .split-legend { display: flex; justify-content: space-between; gap: 6px;
        flex-wrap: wrap; font-size: 11px; margin: 0 0 10px; }
      .split-legend > span { white-space: nowrap; }
      .split-legend > span::before { content: ''; display: inline-block;
        width: 7px; height: 7px; border-radius: 2px; margin-right: 4px; }
      .split-loading::before { background: #187b79; }
      .split-scree::before { background: #6688b0; }
      .split-score::before { background: #bc8642; }
    ")),
    sliderInput(
      "pca_split_points", "Budget proportions",
      min = 0, max = 300, value = c(100, 200), step = 1,
      ticks = FALSE, dragRange = FALSE
    ),
    div(
      class = "split-legend",
      span(class = "split-loading", "Loading ", strong(id = "loading_percent", "33.3%")),
      span(class = "split-scree", "Scree ", strong(id = "scree_percent", "33.3%")),
      span(class = "split-score", "Score ", strong(id = "score_percent", "33.3%"))
    ),
    tags$script(HTML("
      $(document).on('shiny:connected', function() {
        const input = $('#pca_split_points');
        const slider = input.data('ionRangeSlider');
        if (!slider) return;

        // Prevent any stage from receiving a zero share.
        slider.update({ from_min: 1, to_max: 299, min_interval: 1 });

        function paintShares() {
          const left = slider.result.from / 3;
          const right = slider.result.to / 3;
          input.closest('.pca-split').find('.irs-line').css('background',
            'linear-gradient(to right, ' +
            '#187b79 0%, #187b79 ' + left + '%, ' +
            '#6688b0 ' + left + '%, #6688b0 ' + right + '%, ' +
            '#bc8642 ' + right + '%, #bc8642 100%)');
          $('#loading_percent').text(left.toFixed(1) + '%');
          $('#scree_percent').text((right - left).toFixed(1) + '%');
          $('#score_percent').text((100 - right).toFixed(1) + '%');
        }

        // Namespacing keeps one handler if the browser reconnects.
        input.off('change.budgetShares').on('change.budgetShares', paintShares);
        paintShares();
      });
    "))
  )
}

# Shared controls for the two views that estimate eigenvalues.
.app_scree_controls <- function() {
  tags$details(
    class = "app-details",
    tags$summary("Scree method settings"),
    div(
      class = "details-body",
      conditionalPanel(
        condition = "(input.analysis == 'pca' && input.pca_scree_method == 'clipped') || (input.analysis == 'scree' && input.scree_methods && input.scree_methods.indexOf('clipped') >= 0)",
        div(class = "control-caption", "Clipped"),
        numericInput("C_clip", "Clipping threshold", value = 3, min = 1e-6, step = 0.5)
      ),
      conditionalPanel(
        condition = "(input.analysis == 'pca' && input.pca_scree_method == 'pmwm') || (input.analysis == 'scree' && input.scree_methods && input.scree_methods.indexOf('pmwm') >= 0)",
        div(class = "control-caption", "PMWM"),
        numericInput("pmwm_beta", "Beta", value = 1.001, min = 1.000001, step = 0.001),
        div(
          class = "control-pair",
          numericInput("pmwm_a", "Lower bound (a)", value = 0, step = 1),
          numericInput("pmwm_b", "Upper bound (b)", value = 10, step = 1)
        ),
        numericInput("pmwm_trim_const", "Trim constant", value = 10, min = 0, step = 1),
        numericInput("pmwm_eta", "Eta", value = 0.01, min = 0, max = 0.49, step = 0.01),
        checkboxInput("pmwm_split_mode", "Split mode", value = TRUE)
      ),
      conditionalPanel(
        condition = "(input.analysis == 'pca' && input.pca_scree_method == 'huber') || (input.analysis == 'scree' && input.scree_methods && input.scree_methods.indexOf('huber') >= 0)",
        div(class = "control-caption", "Huber"),
        div(
          class = "control-pair",
          numericInput("huber_mu0", "Initial mean", value = 0, step = 0.5),
          numericInput("huber_eta0", "Initial scale", value = 1, min = 1e-8, step = 0.1)
        ),
        div(
          class = "control-pair",
          numericInput("huber_T", "Iterations (T)", value = 50, min = 1, step = 1),
          numericInput("huber_M", "M", value = 20, min = 1, step = 1)
        ),
        div(
          class = "control-pair",
          numericInput("huber_k_min_m2", "Lower search bound", value = -20, step = 1),
          numericInput("huber_k_max_m2", "Upper search bound", value = 20, step = 1)
        ),
        numericInput("huber_m2_frac", "Second-moment budget share", value = 0.25, min = 1e-6, max = 0.99, step = 0.05)
      ),
      checkboxInput("mono", "Monotone eigenvalue post-processing", value = TRUE)
    )
  )
}

ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      :root {
        --navy: #172d48;
        --teal: #187b79;
        --ink: #24364a;
        --muted: #62758a;
        --line: #dce5ec;
        --canvas: #f3f6f8;
      }
      body {
        background: var(--canvas);
        color: var(--ink);
        font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif;
        font-size: 14px;
        line-height: 1.5;
      }
      .container-fluid { max-width: 1760px; padding: 0 28px 36px; }
      .app-header {
        display: flex;
        align-items: center;
        gap: 17px;
        padding: 29px 0 25px;
      }
      .app-wordmark {
        background: var(--navy);
        color: #fff;
        font-size: 17px;
        font-weight: 700;
        letter-spacing: .04em;
        padding: 15px 16px;
        border-radius: 12px;
      }
      .app-header h1 { margin: 0 0 4px; color: var(--navy); font-size: 25px; font-weight: 650; }
      .app-header p { margin: 0; color: var(--muted); font-size: 14px; }
      .well {
        background: #fff;
        border: 1px solid var(--line);
        border-radius: 14px;
        box-shadow: 0 3px 14px rgba(23, 45, 72, .035);
        padding: 20px;
      }
      .sidebar-heading {
        margin: 0 0 13px;
        font-size: 12px;
        font-weight: 750;
        letter-spacing: .085em;
        text-transform: uppercase;
        color: var(--navy);
      }
      .control-section + .control-section { border-top: 1px solid var(--line); margin-top: 20px; padding-top: 20px; }
      .form-group { margin-bottom: 13px; }
      .control-label { font-size: 12px; font-weight: 600; }
      .form-control, .selectize-input { border-color: #d2dce5; border-radius: 6px; box-shadow: none; }
      .form-control:focus, .selectize-input.focus { border-color: var(--teal); box-shadow: 0 0 0 2px rgba(24, 123, 121, .1); }
      .checkbox label, .radio label { font-size: 12px; line-height: 1.5; }
      input[type='checkbox'], input[type='radio'] { accent-color: var(--teal); }
      .help-block { color: var(--muted); font-size: 11px; line-height: 1.55; }
      .control-pair { display: grid; grid-template-columns: minmax(0, 1fr) minmax(0, 1fr); gap: 12px; }
      .control-pair > .shiny-input-container { width: 100%; min-width: 0; }
      .control-caption { color: var(--teal); font-size: 12px; font-weight: 700; margin: 8px 0 10px; }
      .app-details { border-top: 1px solid var(--line); margin: 14px 0 0; }
      .app-details > summary {
        display: flex;
        align-items: center;
        gap: 10px;
        padding: 12px 5px;
        border-radius: 5px;
        font-size: 12px;
        font-weight: 650;
        color: var(--navy);
        cursor: pointer;
        list-style: none;
      }
      .app-details > summary::-webkit-details-marker { display: none; }
      .app-details > summary::marker { content: ''; }
      .app-details > summary::before {
        content: '';
        flex: 0 0 7px;
        width: 7px;
        height: 7px;
        border-right: 2px solid var(--teal);
        border-bottom: 2px solid var(--teal);
        transform: rotate(-45deg);
      }
      .app-details[open] > summary::before { transform: rotate(45deg); }
      .app-details > summary:hover { background: #eef5f6; }
      .app-details > summary:focus-visible {
        background: #eef5f6;
        outline: 2px solid var(--teal);
        outline-offset: 2px;
      }
      .preprocessing-details { border-top: 0; margin: 0; }
      .details-body { padding: 1px 0 7px; }
      .run-area { border-top: 1px solid var(--line); padding-top: 18px; margin-top: 20px; }
      .btn-primary { background: var(--teal); border-color: var(--teal); font-weight: 650; }
      .btn-primary:hover, .btn-primary:focus, .btn-primary:active { background: #116563; border-color: #116563; }
      #run { width: 100%; min-height: 42px; border-radius: 7px; }
      .nav-pills { background: #e7edf2; padding: 5px; border-radius: 11px; margin-bottom: 20px; display: flex; gap: 4px; }
      .nav-pills > li { float: none; flex: 1; text-align: center; margin: 0 !important; }
      .nav-pills > li > a { border-radius: 7px; color: #566a80; font-size: 13px; font-weight: 650; padding: 11px 13px; }
      .nav-pills > li.active > a, .nav-pills > li.active > a:hover, .nav-pills > li.active > a:focus {
        background: var(--navy); color: #fff; box-shadow: 0 2px 5px rgba(23, 45, 72, .12);
      }
      .nav-pills > li > a:hover { background: #dce6ee; }
      .result-card { background: #fff; border: 1px solid var(--line); border-radius: 14px; padding: 23px 24px; box-shadow: 0 3px 14px rgba(23, 45, 72, .035); margin-bottom: 20px; }
      .card-heading { margin: 0 0 5px; color: var(--navy); font-size: 21px; font-weight: 650; }
      .result-status { border-left: 3px solid #89c0bb; background: #f4f9f9; color: #4f6a73; font-size: 12px; padding: 9px 12px; margin: 13px 0 19px; border-radius: 0 6px 6px 0; }
      .plot-toolbar { display: flex; align-items: end; flex-wrap: wrap; column-gap: 24px; row-gap: 4px; border-bottom: 1px solid var(--line); padding-bottom: 4px; margin-bottom: 18px; }
      .plot-toolbar > .shiny-input-container { width: 230px; }
      .plot-toolbar .shiny-input-radiogroup { width: auto; }
      .plot-toolbar .shiny-input-radiogroup label { white-space: nowrap; }
      .plot-toolbar > .shiny-panel-conditional { max-width: 100%; }
      .toolbar-alpha { width: 210px; max-width: 100%; }
      .toolbar-alpha .shiny-input-container { width: 100%; }
      .plot-settings-grid { display: grid; grid-template-columns: repeat(3, minmax(0, 1fr)); gap: 5px 22px; }
      .plot-settings-grid .shiny-input-container { max-width: 100%; }
      .plot-area { min-height: 360px; overflow: auto; }
      .scree-plot-area { width: 100%; max-width: 720px; margin: 0 auto; }
      .plot-area .shiny-output-error-validation { padding: 50px 20px; font-size: 14px; color: var(--muted); text-align: center; }
      .table-wrap { overflow-x: auto; }
      .table { font-size: 12px; }
      .table > thead > tr > th { border-bottom: 1px solid var(--line); color: var(--navy); }
      .summary-heading { font-size: 14px; color: var(--navy); font-weight: 650; margin: 23px 0 12px; }
      .fit-info { color: var(--muted); font-size: 12px; margin-top: 18px; }
      .irs--shiny .irs-bar, .irs--shiny .irs-single { background: var(--teal); border-color: var(--teal); }
      @media (max-width: 1100px) {
        .container-fluid { padding-right: 18px; padding-left: 18px; }
        .well { padding: 16px; }
        .control-pair { grid-template-columns: 1fr; gap: 0; }
        .plot-settings-grid { grid-template-columns: repeat(2, minmax(0, 1fr)); }
        .result-card { padding: 20px; }
      }
      @media (max-width: 767px) {
        .app-header { align-items: flex-start; padding-top: 22px; gap: 12px; }
        .app-header h1 { font-size: 22px; }
        .app-header p { font-size: 12px; }
        .app-wordmark { font-size: 14px; padding: 13px 11px; }
        .control-pair, .plot-settings-grid { grid-template-columns: repeat(2, minmax(0, 1fr)); }
        .nav-pills > li > a { padding: 10px 6px; font-size: 12px; }
        .result-card { padding: 18px 15px; }
        .plot-toolbar { gap: 4px 16px; }
      }
      @media (max-width: 420px) {
        .control-pair, .plot-settings-grid { grid-template-columns: 1fr; }
      }
    "))
  ),

  div(
    class = "app-header",
    div(class = "app-wordmark", "dppca"),
    div(
      h1("DP-PCA"),
      p("Explore private principal components, from individual estimates to a complete analysis.")
    )
  ),

  sidebarLayout(
    sidebarPanel(
      width = 3,
      div(
        class = "control-section",
        h3(class = "sidebar-heading", "Data"),
        radioButtons(
          "data_source", NULL,
          choices = if (.has_supplied_data) {
            c("Provided data" = "provided", "Example data" = "example", "Upload CSV" = "upload")
          } else {
            c("Example data" = "example", "Upload CSV" = "upload")
          },
          selected = if (.has_supplied_data) "provided" else "example"
        ),
        conditionalPanel(
          condition = "input.data_source == 'provided'",
          helpText("Using the data supplied to dppca_app().")
        ),
        conditionalPanel(
          condition = "input.data_source == 'example'",
          selectInput(
            "example_dataset", "Example dataset",
            choices = c("Gaussian" = "gau", "Gaussian with groups" = "gau_grouped", "Adult" = "adult"),
            selected = "gau_grouped"
          )
        ),
        conditionalPanel(
          condition = "input.data_source == 'upload'",
          fileInput("csv_file", "CSV file", accept = ".csv"),
          checkboxInput("csv_header", "First row contains column names", value = TRUE)
        ),
        selectInput(
          "group_col", "Group column",
          choices = c("None" = ""), selected = ""
        ),
        checkboxInput(
          "use_group", "Use group labels for score plot", value = TRUE
        )
      ),

      div(
        class = "control-section",
        h3(class = "sidebar-heading", "Budget allocation"),
        conditionalPanel(
          condition = "input.analysis != 'loading'",
          radioButtons(
            "budget_mode", "Budget input",
            choices = c(
              "Total + proportions" = "total",
              "Set each component" = "direct"
            ),
            selected = "total"
          )
        ),
        conditionalPanel(
          condition = "input.analysis == 'loading' || input.budget_mode == 'total'",
          div(
            class = "control-pair",
            numericInput(
              "eps", "Total epsilon",
              value = 10, min = 1e-6, step = 0.1
            ),
            numericInput(
              "delta", "Total delta",
              value = 1e-8, min = 0, max = 0.99, step = 1e-6
            )
          )
        ),
        conditionalPanel(
          condition = "input.analysis == 'loading'",
          helpText("Loading uses the full budget to estimate all PC directions.")
        ),
        conditionalPanel(
          condition = "input.budget_mode == 'total' && (input.analysis == 'scree' || input.analysis == 'score')",
          sliderInput(
            "loading_share", "Loading share",
            min = 0.05, max = 0.95, value = 0.5, step = 0.05
          ),
          helpText("The remaining share is used for the selected scree or score method.")
        ),
        conditionalPanel(
          condition = "input.analysis == 'pca' && input.budget_mode == 'total'",
          .app_pca_split_control()
        ),
        conditionalPanel(
          condition = "input.analysis != 'loading' && input.budget_mode == 'direct'",
          div(class = "control-caption", "Loading"),
          div(
            class = "control-pair",
            numericInput(
              "budget_loading_eps", "Epsilon",
              value = 10 / 3, min = 0, step = 0.1
            ),
            numericInput(
              "budget_loading_delta", "Delta",
              value = 1e-6, min = 0, max = 0.99, step = 1e-6
            )
          ),
          conditionalPanel(
            condition = "input.analysis == 'pca' || input.analysis == 'scree'",
            div(class = "control-caption", "Scree"),
            div(
              class = "control-pair",
              numericInput(
                "budget_scree_eps", "Epsilon",
                value = 10 / 3, min = 0, step = 0.1
              ),
              numericInput(
                "budget_scree_delta", "Delta",
                value = 1e-6, min = 0, max = 0.99, step = 1e-6
              )
            )
          ),
          conditionalPanel(
            condition = "input.analysis == 'pca' || input.analysis == 'score'",
            div(class = "control-caption", "Score"),
            div(
              class = "control-pair",
              numericInput(
                "budget_score_eps", "Epsilon",
                value = 10 / 3, min = 0, step = 0.1
              ),
              numericInput(
                "budget_score_delta", "Delta",
                value = 1e-6, min = 0, max = 0.99, step = 1e-6
              )
            )
          ),
          helpText("The total budget is derived from the component budgets above.")
        ),
        div(class = "table-wrap", tableOutput("budget_preview"))
      ),

      div(
        class = "control-section",
        tags$details(
          class = "app-details preprocessing-details",
          tags$summary("Preprocessing"),
          div(
            class = "details-body",
            checkboxInput("center", "Center columns", value = TRUE),
            checkboxInput("standardize", "Standardize columns", value = FALSE),
            checkboxInput("cpp", "Use C++ implementation", value = TRUE)
          )
        )
      ),

      div(
        class = "control-section",
        h3(class = "sidebar-heading", "Analysis settings"),
        conditionalPanel(
          condition = "input.analysis == 'pca' || input.analysis == 'scree'",
          numericInput("k", "Leading PCs", value = 4, min = 1, step = 1)
        ),
        conditionalPanel(
          condition = "input.analysis != 'scree'",
          div(
            class = "control-pair",
            numericInput("axis_x", "x-axis PC", value = 1, min = 1, step = 1),
            numericInput("axis_y", "y-axis PC", value = 2, min = 1, step = 1)
          )
        ),
        conditionalPanel(
          condition = "input.analysis == 'pca'",
          selectInput(
            "pca_scree_method", "Scree method",
            choices = c("Clipped" = "clipped", "PMWM" = "pmwm", "Huber" = "huber"),
            selected = "clipped"
          )
        ),
        conditionalPanel(
          condition = "input.analysis == 'scree'",
          checkboxGroupInput(
            "scree_methods", "Scree methods",
            choices = c("Clipped" = "clipped", "PMWM" = "pmwm", "Huber" = "huber"),
            selected = "clipped"
          )
        ),
        conditionalPanel(
          condition = "input.analysis == 'pca' || input.analysis == 'scree'",
          .app_scree_controls()
        ),
        conditionalPanel(
          condition = "input.analysis == 'pca'",
          selectInput(
            "pca_score_method", "Score method",
            choices = c("Sparse histogram" = "sparse", "Additive histogram" = "add"),
            selected = "sparse"
          )
        ),
        conditionalPanel(
          condition = "input.analysis == 'score'",
          checkboxGroupInput(
            "score_methods", "Score methods",
            choices = c("Sparse histogram" = "sparse", "Additive histogram" = "add"),
            selected = "sparse"
          )
        ),
        conditionalPanel(
          condition = "input.analysis == 'pca' || input.analysis == 'score'",
          div(
            class = "control-pair",
            numericInput("bins_x", "x-axis bins", value = 10, min = 1, step = 1),
            numericInput("bins_y", "y-axis bins", value = 10, min = 1, step = 1)
          ),
          conditionalPanel(
            condition = "input.analysis == 'score'",
            checkboxGroupInput(
              "score_display", "Score displays",
              choices = c("Histogram" = "histogram", "Sampled points" = "sample"),
              selected = "histogram"
            ),
            numericInput(
              "sampling_size", "Sampled points",
              value = 500, min = 1, step = 100
            ),
            helpText("Sample size takes effect on the next run. Histogram and sample views reuse the same estimate.")
          )
        ),
        conditionalPanel(
          condition = "input.analysis == 'loading'",
          checkboxGroupInput(
            "loading_display", "Loading displays",
            choices = c("Vector" = "vector", "Point" = "point", "Heatmap" = "heatmap"),
            selected = "vector", inline = TRUE
          ),
          numericInput("variables", "Top variables", value = 10, min = 1, step = 1),
          checkboxInput("labels", "Show variable labels", value = TRUE),
          checkboxInput("align_sign", "Align non-private signs to DP directions", value = TRUE),
          conditionalPanel(
            condition = "input.loading_display && input.loading_display.indexOf('heatmap') >= 0",
            selectizeInput(
              "heatmap_components", "Heatmap PCs",
              choices = NULL, multiple = TRUE
            ),
            helpText("Leave empty to use the two selected axis PCs.")
          )
        )
      ),
      div(class = "run-area", actionButton("run", "Run PCA", class = "btn-primary"))
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "analysis", type = "pills", selected = "pca",
        tabPanel(
          "PCA", value = "pca",
          div(
            class = "result-card",
            h2(class = "card-heading", "A complete PCA analysis"),
            div(class = "result-status", uiOutput("pca_status")),
            div(
              class = "plot-toolbar",
              selectInput("pca_view", "View", choices = c("Overview" = "all", "Scree" = "scree", "Score" = "score", "Loading" = "loading", "Variable" = "variable", "Biplot" = "biplot", "Summary" = "summary"), selected = "all"),
              conditionalPanel(
                condition = "input.pca_view == 'all'",
                radioButtons("pca_version", "Overview results", choices = c("DP" = "private", "Non-private" = "nonprivate"), selected = "private", inline = TRUE)
              ),
              conditionalPanel(
                condition = "input.pca_view != 'all' && input.pca_view != 'summary'",
                checkboxInput("pca_compare", "Compare with non-private results", value = TRUE)
              ),
              conditionalPanel(
                condition = "input.pca_view == 'all' || input.pca_view == 'biplot'",
                div(
                  class = "toolbar-alpha",
                  sliderInput(
                    "alpha", "Biplot alpha",
                    min = 0, max = 1, value = 0, step = 0.01,
                    ticks = FALSE, width = "100%"
                  )
                )
              )
            ),
            conditionalPanel(
              condition = "input.pca_view != 'summary'",
              tags$details(
                class = "app-details",
                tags$summary("Plot settings"),
                div(
                  class = "details-body plot-settings-grid",
                  numericInput("pca_variables", "Top variables", value = 10, min = 1, step = 1),
                  checkboxInput("pca_labels", "Show variable labels", value = TRUE),
                  selectInput("pca_loading_display", "Loading display", choices = c("Vector" = "vector", "Point" = "point", "Heatmap" = "heatmap"), selected = "point"),
                  selectInput("pca_score_display", "Score display", choices = c("Histogram" = "histogram", "Sampled points" = "sample"), selected = "histogram"),
                  selectInput("pca_variable_display", "Variable display", choices = c("Vector" = "vector", "Point" = "point"), selected = "vector"),
                  checkboxInput("pca_circle", "Show reference circle", value = TRUE)
                )
              ),
              div(class = "plot-area", plotOutput("pca_plot", height = "auto"))
            ),
            conditionalPanel(
              condition = "input.pca_view == 'summary'",
              h3(class = "summary-heading", "Component importance"),
              div(class = "table-wrap", tableOutput("importance")),
              h3(class = "summary-heading", "Privacy allocation"),
              div(class = "table-wrap", tableOutput("pca_budget"))
            ),
            div(class = "fit-info", uiOutput("fit_info"))
          )
        ),
        tabPanel(
          "Loading", value = "loading",
          div(
            class = "result-card",
            h2(class = "card-heading", "Loading plot"),
            div(class = "result-status", uiOutput("loading_status")),
            div(class = "plot-area", plotOutput("loading_plot", height = "auto"))
          )
        ),
        tabPanel(
          "Scree", value = "scree",
          div(
            class = "result-card",
            h2(class = "card-heading", "Scree plot"),
            div(class = "result-status", uiOutput("scree_status")),
            div(
              class = "plot-toolbar",
              radioButtons(
                "scree_type", "Show",
                choices = c("PVE" = "pve", "Eigenvalues" = "scree"),
                selected = "pve", inline = TRUE
              )
            ),
            div(
              class = "plot-area scree-plot-area",
              plotOutput("scree_plot", height = "auto")
            )
          )
        ),
        tabPanel(
          "Score", value = "score",
          div(
            class = "result-card",
            h2(class = "card-heading", "Score plot"),
            div(class = "result-status", uiOutput("score_status")),
            div(class = "plot-area", plotOutput("score_plot", height = "auto"))
          )
        )
      )
    )
  )
)

# Server --------------------------------------------------------------------

server <- function(input, output, session) {
  results <- reactiveValues(pca = NULL, loading = NULL, scree = NULL, score = NULL)
  errors <- reactiveValues(pca = NULL, loading = NULL, scree = NULL, score = NULL)
  delta_sample_size <- reactiveVal(NULL)
  direct_budgets <- reactiveValues(pca = NULL, scree = NULL, score = NULL)
  active_direct_mode <- reactiveVal(NULL)

  # Clear old fits immediately, even when the new data are not yet valid.
  observeEvent(list(
    input$data_source, input$example_dataset, input$csv_file,
    input$csv_header, input$group_col
  ), {
    for (mode in c("pca", "loading", "scree", "score")) {
      results[[mode]] <- NULL
      errors[[mode]] <- NULL
    }
  }, priority = 100)

  raw_data <- reactive(.app_checked_data({
    if (identical(input$data_source, "provided")) {
      req(.supplied_X)
      return(.app_data(.supplied_X, "Provided data", .supplied_group))
    }
    if (identical(input$data_source, "upload")) {
      req(input$csv_file)
      data <- utils::read.csv(
        input$csv_file$datapath,
        header = isTRUE(input$csv_header),
        check.names = FALSE
      )
      return(.app_data(data, input$csv_file$name))
    }

    req(input$example_dataset)
    if (input$example_dataset == "gau_grouped") {
      data <- .app_load_data("gau_g")
      if (!(is.data.frame(data) || is.matrix(data))) {
        return(.app_data(.app_load_data("gau"), "Gaussian with groups", data))
      }
      return(.app_data(data, "Gaussian with groups"))
    }
    label <- if (input$example_dataset == "gau") "Gaussian" else "Adult"
    .app_data(.app_load_data(input$example_dataset), label)
  }))

  observeEvent(raw_data(), {
    data <- raw_data()
    choices <- c("None" = "", stats::setNames(names(data$data), names(data$data)))
    if (!is.null(data$external_group)) {
      choices <- c(choices, "Provided group labels" = ".provided_group")
    }
    updateSelectInput(
      session, "group_col", choices = choices, selected = data$group_column
    )
  })

  prepared_data <- reactive(.app_checked_data({
    data <- raw_data()
    group <- NULL
    selected <- input$group_col
    if (is.null(selected)) selected <- ""
    features <- data$data

    if (identical(selected, ".provided_group")) {
      group <- data$external_group
    } else if (nzchar(selected)) {
      validate(need(selected %in% names(features), "Choose a group column."))
      group <- features[[selected]]
      features <- features[, names(features) != selected, drop = FALSE]
    }

    X <- .app_matrix(features)
    list(
      X = X,
      group = group,
      label = data$label
    )
  }))

  observeEvent(prepared_data(), {
    data <- prepared_data()
    p <- ncol(data$X)
    n <- nrow(data$X)
    updateNumericInput(session, "k", max = p, value = min(4L, p))
    updateNumericInput(session, "axis_x", max = p, value = 1)
    updateNumericInput(session, "axis_y", max = p, value = 2)
    updateNumericInput(session, "variables", max = p, value = min(10L, p))
    updateNumericInput(session, "pca_variables", max = p, value = min(10L, p))
    updateSelectizeInput(
      session, "heatmap_components",
      choices = stats::setNames(seq_len(p), paste0("PC", seq_len(p))),
      selected = c(1, 2)
    )
    if (!identical(n, delta_sample_size())) {
      updateNumericInput(
        session, "delta",
        value = 10^(-ceiling(log10(as.numeric(n)^2)))
      )
      delta_sample_size(n)
    }
  })

  show_direct_budget <- function(budgets) {
    for (component in names(budgets)) {
      budget <- budgets[[component]]
      updateNumericInput(
        session, paste0("budget_", component, "_eps"), value = budget[["eps"]]
      )
      updateNumericInput(
        session, paste0("budget_", component, "_delta"), value = budget[["delta"]]
      )
    }
    invisible(NULL)
  }

  initialize_direct_budget <- function(mode) {
    if (is.null(mode) || mode == "loading") return(FALSE)
    tryCatch({
      privacy <- .app_privacy(input, mode, budget_mode = "total")
      direct_budgets[[mode]] <- privacy$components
      show_direct_budget(privacy$components)
      TRUE
    }, error = function(e) {
      showNotification(conditionMessage(e), type = "error", duration = 10)
      FALSE
    })
  }

  # Each analysis keeps its own direct entries. First visits start from that
  # analysis's total allocation, so hidden placeholder fields are never reused.
  observeEvent(list(input$budget_mode, input$analysis), {
    previous <- active_direct_mode()
    if (!is.null(previous)) {
      components <- switch(
        previous,
        pca = c("loading", "scree", "score"),
        scree = c("loading", "scree"),
        score = c("loading", "score")
      )
      saved <- lapply(components, function(component) {
        values <- lapply(c("eps", "delta"), function(parameter) {
          value <- input[[paste0("budget_", component, "_", parameter)]]
          if (is.null(value) || !length(value)) NA_real_ else value
        })
        stats::setNames(unlist(values), c("eps", "delta"))
      })
      names(saved) <- components
      direct_budgets[[previous]] <- saved
    }

    active_direct_mode(NULL)
    mode <- input$analysis
    if (
      identical(input$budget_mode, "direct") &&
        !is.null(mode) && mode != "loading"
    ) {
      if (is.null(direct_budgets[[mode]])) {
        if (!initialize_direct_budget(mode)) {
          components <- switch(
            mode,
            pca = c("loading", "scree", "score"),
            scree = c("loading", "scree"),
            score = c("loading", "score")
          )
          empty <- stats::setNames(
            lapply(components, function(component) c(eps = NA_real_, delta = NA_real_)),
            components
          )
          direct_budgets[[mode]] <- empty
          show_direct_budget(empty)
        }
      } else {
        show_direct_budget(direct_budgets[[mode]])
      }
      active_direct_mode(mode)
    }
  })

  observeEvent(input$analysis, {
    label <- switch(
      input$analysis,
      pca = "Run PCA",
      loading = "Run Loading",
      scree = "Run Scree",
      score = "Run Score"
    )
    updateActionButton(session, "run", label = label)
  })

  output$budget_preview <- renderTable({
    req(input$analysis)
    if (input$analysis == "loading") {
      validate(
        need(is.finite(input$eps) && input$eps > 0, "eps must be positive."),
        need(
          is.finite(input$delta) && input$delta > 0 && input$delta < 1,
          "delta must be between 0 and 1."
        )
      )
      return(.app_budget_table(rbind(Loading = c(eps = input$eps, delta = input$delta))))
    }
    privacy <- tryCatch(.app_privacy(input, input$analysis), error = identity)
    validate(need(!inherits(privacy, "error"),
                  if (inherits(privacy, "error")) conditionMessage(privacy) else ""))
    methods <- switch(
      input$analysis,
      scree = length(input$scree_methods),
      score = length(input$score_methods),
      1L
    )
    .app_budget(privacy, input$analysis, methods)
  }, striped = TRUE, spacing = "xs", align = "lrr")

  # Estimation is only triggered by this button. Renderers use saved results.
  observeEvent(input$run, {
    mode <- input$analysis
    errors[[mode]] <- NULL

    computed <- tryCatch({
      if (identical(input$data_source, "upload") && is.null(input$csv_file)) {
        stop("Upload a CSV file before running the analysis.", call. = FALSE)
      }
      data <- prepared_data()
      X <- data$X
      p <- ncol(X)
      privacy <- if (mode == "loading") NULL else .app_privacy(input, mode)
      axes <- NULL
      bins <- NULL
      k <- NULL

      if (mode != "scree") {
        axes <- c(
          .app_integer(input$axis_x, "Horizontal PC", upper = p),
          .app_integer(input$axis_y, "Vertical PC", upper = p)
        )
        if (axes[1] == axes[2]) stop("Choose two different PC axes.", call. = FALSE)
      }
      if (mode %in% c("pca", "scree")) {
        k <- .app_integer(input$k, "Number of scree components", upper = p)
        if (mode == "pca" && (k < 2L || max(axes) > k)) {
          stop("For PCA, k must include both selected axes and be at least 2.",
               call. = FALSE)
        }
      }
      if (mode %in% c("pca", "score")) {
        bins <- c(
          .app_integer(input$bins_x, "Horizontal bins"),
          .app_integer(input$bins_y, "Vertical bins")
        )
      }
      if (mode == "loading" && !length(input$loading_display)) {
        stop("Choose at least one loading display.", call. = FALSE)
      }
      if (mode == "score" && !length(input$score_display)) {
        stop("Choose at least one score display.", call. = FALSE)
      }

      value <- withProgress(message = paste("Computing", toupper(mode)), value = 0, {
        incProgress(0.1, detail = "Estimating private quantities")
        if (mode == "pca") {
          value <- dp_pca(
            X = X,
            privacy = privacy,
            scree = list(
              k = k,
              method = input$pca_scree_method,
              control = .app_controls(input, input$pca_scree_method),
              mono = isTRUE(input$mono)
            ),
            score = list(method = input$pca_score_method, bins = bins),
            axes = axes,
            center = isTRUE(input$center),
            standardize = isTRUE(input$standardize),
            cpp.option = isTRUE(input$cpp),
            non_private = TRUE
          )
        } else if (mode == "loading") {
          components <- as.integer(input$heatmap_components)
          if (!length(components)) components <- axes
          value <- dp_loading_plot(
            X = X,
            eps = input$eps,
            delta = input$delta,
            axes = axes,
            display = c("vector", "point", "heatmap"),
            heatmap_components = components,
            center = isTRUE(input$center),
            standardize = isTRUE(input$standardize),
            cpp.option = isTRUE(input$cpp),
            variables = .app_integer(input$variables, "Variables", upper = p),
            labels = isTRUE(input$labels),
            align_sign = isTRUE(input$align_sign)
          )
        } else if (mode == "scree") {
          methods <- unique(input$scree_methods)
          value <- dp_scree(
            X = X,
            k = k,
            privacy = privacy,
            method = methods,
            control = .app_controls(input, methods),
            center = isTRUE(input$center),
            standardize = isTRUE(input$standardize),
            cpp.option = isTRUE(input$cpp),
            mono = isTRUE(input$mono)
          )
        } else {
          methods <- unique(input$score_methods)
          if (!length(methods)) stop("Choose at least one score method.", call. = FALSE)
          args <- list(
            X = X,
            privacy = privacy,
            bins = bins,
            method = methods,
            private_plot = c("histogram", "sample"),
            sampling_control = sampling_control(
              method = "uniform",
              sample_size = .app_integer(input$sampling_size, "Sample points")
            ),
            center = isTRUE(input$center),
            standardize = isTRUE(input$standardize),
            cpp.option = isTRUE(input$cpp),
            axes = axes
          )
          if (isTRUE(input$use_group) && !is.null(data$group)) {
            args$group <- data$group
            value <- do.call(dp_score_plot_group, args)
          } else {
            value <- do.call(dp_score_plot, args)
          }
        }
        incProgress(0.9, detail = "Saving results for display")
        value
      })

      list(
        value = value,
        seed = sample.int(.Machine$integer.max, 1L),
        label = data$label,
        n = nrow(X),
        p = p,
        time = format(Sys.time(), "%H:%M:%S"),
        privacy = privacy
      )
    }, error = function(e) {
      errors[[mode]] <- conditionMessage(e)
      showNotification(conditionMessage(e), type = "error", duration = 10)
      NULL
    })

    results[[mode]] <- computed
  }, ignoreInit = TRUE)

  for (mode in c("pca", "loading", "scree", "score")) {
    local({
      current_mode <- mode
      output[[paste0(current_mode, "_status")]] <- renderUI({
        error <- errors[[current_mode]]
        if (!is.null(error)) {
          return(tags$p(class = "result-error", error))
        }
        .app_result_note(results[[current_mode]])
      })
    })
  }

  # Same stored seed = same synthetic points; appearance changes spend no
  # additional estimation budget and do not change the fitted quantities.
  pca_base_plots <- reactive({
    result <- results$pca
    req(result)
    variables <- .app_integer(input$pca_variables, "Variables", upper = result$p)
    .app_pca_plots(
      fit = result$value,
      seed = result$seed,
      alpha = 0,
      variables = variables,
      labels = isTRUE(input$pca_labels),
      loading_display = input$pca_loading_display,
      score_display = input$pca_score_display,
      variable_display = input$pca_variable_display,
      circle = isTRUE(input$pca_circle)
    )
  })

  biplot_plots <- reactive({
    result <- results$pca
    req(result)
    .app_biplot_plots(
      fit = result$value,
      seed = result$seed,
      alpha = input$alpha,
      variables = .app_integer(input$pca_variables, "Variables", upper = result$p),
      labels = isTRUE(input$pca_labels)
    )
  })

  output$pca_plot <- renderPlot({
    validate(need(!is.null(results$pca), "Run PCA to see the plots."))
    req(input$pca_view != "summary")
    if (input$pca_view == "biplot") {
      plots <- biplot_plots()
    } else {
      plots <- pca_base_plots()
      if (input$pca_view == "all") {
        biplots <- biplot_plots()
        for (version in c("nonprivate", "private")) {
          plots[[version]]$plots$biplot <- biplots[[version]]$plots$biplot
        }
      }
    }
    .app_draw_pca(
      plots,
      view = input$pca_view,
      version = if (input$pca_view == "all") input$pca_version else "private",
      compare = isTRUE(input$pca_compare)
    )
  }, height = function() {
    if (identical(input$pca_view, "all")) 840 else 590
  }, res = 110)

  output$fit_info <- renderUI({
    req(results$pca)
    fit <- results$pca$value
    tags$p(
      class = "result-note",
      paste0(
        "Scree: ", fit$methods$scree, " · Score: ", fit$methods$score,
        " · Axes: PC", fit$axes[1], ", PC", fit$axes[2],
        " · PVE is normalized over ", fit$k, " estimated components."
      )
    )
  })

  pca_summary <- reactive({
    req(results$pca)
    .app_fun("summary.dp_pca")(results$pca$value)
  })

  output$importance <- renderTable({
    values <- pca_summary()$importance
    data.frame(
      Component = colnames(values),
      Eigenvalue = signif(values["Eigenvalue", ], 5),
      `PVE (%)` = round(100 * values["PVE", ], 2),
      `Cumulative PVE (%)` = round(100 * values["Cumulative PVE", ], 2),
      row.names = NULL,
      check.names = FALSE
    )
  }, striped = TRUE, spacing = "s", digits = 5)

  output$pca_budget <- renderTable({
    .app_budget_table(pca_summary()$budget)
  }, striped = TRUE, spacing = "s", align = "lrr")

  output$loading_plot <- renderPlot({
    validate(need(!is.null(results$loading), "Run Loading to compare the loadings."))
    validate(need(length(input$loading_display) > 0, "Choose at least one loading display."))
    .app_draw_loading(results$loading$value, display = input$loading_display)
  }, height = function() {
    displays <- input$loading_display
    if (length(displays) > 1L) {
      if (length(displays) == 2L) 1100 else 900
    } else if ("heatmap" %in% displays) {
      690
    } else {
      590
    }
  }, res = 110)

  output$scree_plot <- renderPlot({
    validate(need(!is.null(results$scree), "Run Scree to compare the methods."))
    .app_draw_scree(results$scree$value, type = input$scree_type)
  }, height = function() {
    width <- session$clientData$output_scree_plot_width
    if (length(width) != 1L || !is.finite(width) || width <= 0) return(640)
    max(360, min(720, round(width)))
  }, res = 110)

  output$score_plot <- renderPlot({
    validate(need(!is.null(results$score), "Run Score to see the distributions."))
    validate(need(length(input$score_display) > 0, "Choose at least one score display."))
    .app_draw_score(results$score$value, display = input$score_display)
  }, height = function() {
    if (length(input$score_display) > 1L) 1050 else 590
  }, res = 110)
}

shinyApp(ui, server)
