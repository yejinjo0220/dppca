library(shiny)
library(dppca)

# Requires the S3-enabled dppca sources supplied with this update.
# Estimate on Run; plot the stored fits through their S3 methods.

# Data and small application helpers -----------------------------------------

.supplied_X <- getOption("dppca.app.X", NULL)
.supplied_group <- getOption("dppca.app.group", NULL)
.has_supplied_data <- !is.null(.supplied_X)

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
  # Numeric preview for UI limits; estimators receive the original data frame
  # so their common preprocessing can report excluded columns.
  numeric_columns <- vapply(data, function(x) {
    is.numeric(x) && !is.complex(x) && is.null(dim(x))
  }, logical(1))
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
      points[1] <= 0 || points[2] >= 100 || points[1] >= points[2]
  ) {
    stop("Keep a positive share for Loading, Scree, and Score.", call. = FALSE)
  }

  # The two slider handles divide 100% into three strictly positive shares.
  c(
    loading = points[1],
    scree = points[2] - points[1],
    score = 100 - points[2]
  ) / 100
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
    scree = c(loading = input$scree_loading_share / 100,
              scree = 1 - input$scree_loading_share / 100),
    score = c(loading = input$score_loading_share / 100,
             score = 1 - input$score_loading_share / 100)
  )
  privacy_control(eps = input$eps, delta = input$delta, split = split)
}

.app_budget <- function(privacy, mode, methods = 1L) {
  b <- privacy$budgets
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
  tagList(
    tags$p(
      class = "result-note",
      paste0(
        result$label, " · ", result$n, " observations · ", result$p,
        " numeric variables"
      )
    ),
    if (length(result$warnings)) tags$div(
      class = "result-note", style = "white-space: pre-line;",
      paste0("Warning:\n", paste(result$warnings, collapse = "\n\n"))
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

# Capture the S3 method's return value without drawing on another device.
# Reusing a run's seed keeps synthetic points stable when switching views.
.app_stored_plots <- function(result, center = NULL) {
  plot_args <- result$plot_args
  if (!is.null(center)) plot_args$center <- center
  .app_with_seed(result$seed, .app_capture_plot(
    do.call(plot, c(list(x = result$value), plot_args))
  ))
}

.app_pca_plots <- function(
    fit,
    seed,
    variables = 10,
    biplot_labels = TRUE,
    loading_display = "vector",
    score_display = "histogram",
    center = TRUE,
    sampling = NULL,
    base_size = 12
) {
  .app_with_seed(seed, .app_capture_plot(
    plot(
      x = fit,
      display = "all",
      loading = list(type = loading_display),
      score = list(
        type = score_display,
        plot_control = score_plot_control(show_group_legend = TRUE)
      ),
      biplot = list(alpha = 0, labels = biplot_labels),
      center = center,
      sampling = sampling,
      variables = variables,
      base_size = base_size
    )
  ))
}

.app_biplot_plots <- function(
    fit,
    seed,
    alpha = 0,
    variables = 10,
    biplot_labels = TRUE,
    center = TRUE,
    sampling = NULL,
    base_size = 12
) {
  # Both plot types sample private points before constructing any panels.
  # Default sampling uses the stored histogram and leaves reference points intact.
  .app_with_seed(seed, .app_capture_plot(
    plot(
      x = fit,
      display = "biplot",
      score = list(plot_control = score_plot_control(show_group_legend = TRUE)),
      biplot = list(alpha = alpha, labels = biplot_labels),
      center = center,
      sampling = sampling,
      variables = variables,
      base_size = base_size
    )
  ))
}

# Extract one compact, multi-column guide without altering a stored plot.
.app_right_legend <- function(plot, title = NULL) {
  if (!inherits(plot, "ggplot")) return(NULL)
  group <- plot$data$group
  if (is.null(group) || identical(plot$theme$legend.position, "none")) return(NULL)
  guide_plot <- plot + ggplot2::theme(
    legend.position = "right",
    legend.text = ggplot2::element_text(size = 9),
    legend.title = ggplot2::element_text(size = 10, face = "bold"),
    legend.key.height = grid::unit(3, "mm"),
    legend.key.width = grid::unit(3, "mm")
  ) + ggplot2::guides(
    colour = ggplot2::guide_legend(ncol = 1L),
    fill = ggplot2::guide_legend(ncol = 1L)
  )
  if (!is.null(title)) {
    guide_plot <- guide_plot + ggplot2::labs(colour = title, fill = title)
  }
  g <- ggplot2::ggplotGrob(guide_plot)
  guides <- g$grobs[grepl("^guide-box", g$layout$name)]
  guides <- Filter(function(x) !inherits(x, "zeroGrob"), guides)
  if (length(guides)) guides[[1L]] else NULL
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
    if (is.null(result[[version]])) version <- "private"
    plots <- result[[version]]$plots[types]
    if (length(plots) != length(types) || any(vapply(plots, is.null, logical(1)))) {
      stop("The stored PCA result does not contain all five plots.", call. = FALSE)
    }
    # Keep one vertical group legend immediately to the right of the biplot.
    # Hide only the overview's score legend; extracted plots stay unchanged.
    if (!is.null(plots$biplot$data$group)) {
      plots$score <- plots$score + ggplot2::theme(legend.position = "none")
    }
    top <- patchwork::wrap_plots(
      plots[c("scree", "score", "loading")], ncol = 3L
    )
    bottom <- patchwork::wrap_plots(
      plots[c("variable", "biplot")], ncol = 2L,
      widths = c(1, 1.35)
    )
    output <- patchwork::wrap_plots(
      top, bottom, ncol = 1L, heights = c(1, 1)
    )
  } else {
    if (is.null(result[[version]])) version <- "private"
    versions <- if (isTRUE(compare) && !is.null(result$nonprivate)) {
      c("nonprivate", "private")
    } else {
      version
    }
    plots <- lapply(versions, function(nm) result[[nm]]$plots[[view]])
    if (any(vapply(plots, is.null, logical(1)))) {
      stop("The selected PCA plot is not available in the stored result.", call. = FALSE)
    }
    output <- patchwork::wrap_plots(plots, ncol = length(plots))
  }

  print(output)
  invisible(output)
}

# Preserve the app's square panel while using the package S3 scree method.
.app_draw_scree <- function(
    result,
    type = "pve"
) {
  previous_par <- graphics::par(pty = "s")
  on.exit(graphics::par(previous_par), add = TRUE)
  plot(result, type = type)
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
  has_reference <- !is.null(stored$nonprivate)
  if (!length(methods)) {
    stop("The selected score plots are not available in the stored result.",
         call. = FALSE)
  }

  # With one method and no reference, display histogram and sample side
  # by side. Otherwise, stack each method's panels vertically, and place
  # that method's single group legend beside its sampled-point panel.
  keys <- c(if (has_reference) "reference", methods)
  single_pair <- length(keys) == 1L &&
    identical(display, c("histogram", "sample"))
  columns <- lapply(keys, function(key) {
    method_plots <- lapply(display, function(tp) {
      if (key == "reference") {
        if (tp == "histogram") return(stored$nonprivate)
        sample <- stored$reference_samples$uniform
        if (is.null(sample)) stored$scatter else sample
      } else if (tp == "histogram") {
        stored[[key]]
      } else {
        stored$sample[[key]]$uniform
      }
    })
    if (any(vapply(method_plots, is.null, logical(1)))) {
      stop("The selected score plots are not available in the stored result.",
           call. = FALSE)
    }
    grouped <- !is.null(result$group_colors)
    displayed <- lapply(seq_along(method_plots), function(i) {
      panel <- method_plots[[i]]
      if (grouped && i < length(method_plots)) {
        panel <- panel + ggplot2::theme(legend.position = "none")
      }
      panel
    })
    patchwork::wrap_plots(displayed, ncol = if (single_pair) 2L else 1L)
  })
  output <- patchwork::wrap_plots(columns, ncol = length(columns))
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
  has_reference <- !is.null(stored$nonprivate) &&
    length(stored$nonprivate) > 0L
  plots <- c(
    if (has_reference) unname(stored$nonprivate[display]),
    unname(stored$private[display])
  )

  if (
    length(plots) != (1L + as.integer(has_reference)) * length(display) ||
      any(vapply(plots, is.null, logical(1)))
  ) {
    stop("The selected loading plots are not available in the stored result.",
         call. = FALSE)
  }

  single <- length(display) == 1L
  output <- patchwork::wrap_plots(
    plots,
    ncol = if (single) 1L + as.integer(has_reference) else length(display),
    nrow = if (single) 1L else 1L + as.integer(has_reference),
    byrow = TRUE,
    guides = "collect"
  ) &
    ggplot2::theme(legend.position = "right")
  print(output)
  invisible(output)
}

# A single budget bar, with one movable divider on each side of Scree.
# One draggable allocation slider with three directly editable percentages.
# The slider is the single Shiny input used to compute the actual budget.
.app_pca_split_control <- function() {
  div(
    class = "pca-split",
    tags$style(HTML("
      .pca-split .shiny-input-container { width: 100%; margin-bottom: 0; }
      .pca-split .irs { height: 34px; }
      .pca-split .irs-line { top: 14px; height: 10px; border: 0; }
      .pca-split .irs-bar { top: 14px; height: 10px; background: transparent; border: 0; }
      .pca-split .irs-handle { top: 8px; }
      .pca-split .irs-from, .pca-split .irs-to, .pca-split .irs-single,
      .pca-split .irs-min, .pca-split .irs-max { display: none; }
      .pca-shares { display: grid; grid-template-columns: repeat(3, minmax(0, 1fr));
        gap: 9px; margin: 4px 0 8px; }
      .pca-share label { display: block; font-size: 11px; font-weight: 600;
        color: var(--ink); margin-bottom: 4px; white-space: nowrap; }
      .pca-share label::before { content: ''; display: inline-block; width: 7px;
        height: 7px; border-radius: 2px; margin-right: 5px; }
      .pca-share.loading label::before { background: #187b79; }
      .pca-share.scree label::before { background: #6688b0; }
      .pca-share.score label::before { background: #bc8642; }
      .pca-share .pct-wrap { display: flex; align-items: center; border: 1px solid #d2dce5;
        border-radius: 6px; background: white; padding-right: 7px; }
      .pca-share .pct-wrap:focus-within { border-color: #187b79; }
      .pca-share input { min-width: 0; width: 100%; border: none; border-radius: 6px;
        height: 32px; padding: 5px 3px 5px 8px; font-size: 12px; outline: none; }
      .pca-share .pct-wrap span { font-size: 11px; color: #62758a; }
    ")),
    sliderInput(
      "pca_split_points", "Budget proportions",
      min = 0, max = 100, value = c(33.3, 66.7), step = 0.1,
      ticks = FALSE, dragRange = FALSE
    ),
    div(class = "pca-shares",
      div(class = "pca-share loading", tags$label(`for` = "pca_loading_pct", "Loading"),
        div(class = "pct-wrap", tags$input(id = "pca_loading_pct", type = "number",
          min = "0.1", max = "99.8", step = "0.1", value = "33.3"), span("%"))),
      div(class = "pca-share scree", tags$label(`for` = "pca_scree_pct", "Scree"),
        div(class = "pct-wrap", tags$input(id = "pca_scree_pct", type = "number",
          min = "0.1", max = "99.8", step = "0.1", value = "33.4"), span("%"))),
      div(class = "pca-share score", tags$label(`for` = "pca_score_pct", "Score"),
        div(class = "pct-wrap", tags$input(id = "pca_score_pct", type = "number",
          min = "0.1", max = "99.8", step = "0.1", value = "33.3"), span("%")))
    ),
    tags$script(HTML("
      function initBudgetShares() {
        var input = $('#pca_split_points');
        var slider = input.data('ionRangeSlider');
        if (!slider || input.data('budgetSharesReady')) return;
        input.data('budgetSharesReady', true);
        var ids = ['#pca_loading_pct', '#pca_scree_pct', '#pca_score_pct'];
        var colors = ['#187b79', '#6688b0', '#bc8642'];
        slider.update({ from_min: 0.1, to_max: 99.9, min_interval: 0.1 });

        function currentShares() {
          return [slider.result.from, slider.result.to - slider.result.from,
                  100 - slider.result.to];
        }
        function syncFromSlider() {
          var values = currentShares();
          for (var i = 0; i < 3; i++) $(ids[i]).val(values[i].toFixed(1));
          var left = slider.result.from;
          var right = slider.result.to;
          input.closest('.pca-split').find('.irs-line').css('background',
            'linear-gradient(to right, ' +
            colors[0] + ' 0%, ' + colors[0] + ' ' + left + '%, ' +
            colors[1] + ' ' + left + '%, ' + colors[1] + ' ' + right + '%, ' +
            colors[2] + ' ' + right + '%, ' + colors[2] + ' 100%)');
        }
        input.off('change.budgetShares').on('change.budgetShares', syncFromSlider);
        ids.forEach(function(id, changed) {
          $(id).off('change.budgetShares').on('change.budgetShares', function() {
            var entered = Number($(id).val());
            var shares = currentShares();
            if (!Number.isFinite(entered) || entered < 0.1 || entered > 99.8) {
              syncFromSlider(); return;
            }
            entered = Math.round(entered * 10) / 10;
            if (changed === 0) {
              shares[0] = entered;
              shares[2] = 100 - shares[0] - shares[1];
            } else if (changed === 1) {
              shares[1] = entered;
              shares[2] = 100 - shares[0] - shares[1];
            } else {
              shares[2] = entered;
              shares[1] = 100 - shares[0] - shares[2];
            }
            if (shares.some(function(x) { return x < 0.099999; })) {
              syncFromSlider(); return;
            }
            slider.update({ from: Math.round(shares[0] * 10) / 10,
              to: Math.round((shares[0] + shares[1]) * 10) / 10 });
            // Ion.RangeSlider's programmatic update does not always emit a
            // Shiny change event. Notify the binding explicitly.
            input.trigger('change');
            syncFromSlider();
          });
        });
        syncFromSlider();
      }
      $(initBudgetShares);
      $(document).on('shiny:connected.budgetShares', initBudgetShares);
    "))
  )
}


# Two-part budget control for standalone Scree and Score analyses.
# Shiny uses only the slider value; editable percentages update that slider.
.app_two_split_control <- function(mode) {
  stopifnot(mode %in% c("scree", "score"))
  label <- if (mode == "scree") "Scree" else "Score"
  slider_id <- paste0(mode, "_loading_share")
  loading_id <- paste0(mode, "_loading_pct")
  method_id <- paste0(mode, "_method_pct")
  script <- "
    function initializeTwoShare_@MODE@() {
      var input = $('#@MODE@_loading_share');
      var slider = input.data('ionRangeSlider');
      if (!slider || input.data('twoSharesReady')) return;
      input.data('twoSharesReady', true);
      var loading = $('#@MODE@_loading_pct');
      var method = $('#@MODE@_method_pct');
      var container = input.closest('.two-split');
      var color = '@COLOR@';
      function sync() {
        var left = slider.result.from;
        loading.val(left.toFixed(1));
        method.val((100 - left).toFixed(1));
        container.find('.irs-line').css('background',
          'linear-gradient(to right, #187b79 0%, #187b79 ' + left + '%, ' +
          color + ' ' + left + '%, ' + color + ' 100%)');
      }
      input.off('change.twoShares').on('change.twoShares', sync);
      function update(changed) {
        var entered = Number(changed.val());
        if (!Number.isFinite(entered) || entered < 0.1 || entered > 99.9) {
          sync(); return;
        }
        entered = Math.round(entered * 10) / 10;
        var next = changed.is(loading) ? entered : 100 - entered;
        slider.update({ from: Math.round(next * 10) / 10 });
        input.trigger('change');
        sync();
      }
      loading.off('change.twoShares').on('change.twoShares', function(){ update(loading); });
      method.off('change.twoShares').on('change.twoShares', function(){ update(method); });
      sync();
    }
    $(initializeTwoShare_@MODE@);
    $(document).on('shiny:connected.twoShares', initializeTwoShare_@MODE@);
  "
  script <- gsub("@MODE@", mode, script, fixed = TRUE)
  script <- gsub("@COLOR@", if (mode == "scree") "#6688b0" else "#bc8642", script, fixed = TRUE)
  div(
    class = "pca-split two-split",
    sliderInput(slider_id, "Budget proportions", min = 0.1, max = 99.9,
                value = 50, step = 0.1, ticks = FALSE),
    div(class = "pca-shares two-shares",
      div(class = "pca-share loading",
        tags$label(`for` = loading_id, "Loading"),
        div(class = "pct-wrap", tags$input(id = loading_id, type = "number",
          min = "0.1", max = "99.9", step = "0.1", value = "50.0"), span("%"))),
      div(class = paste("pca-share", mode),
        tags$label(`for` = method_id, label),
        div(class = "pct-wrap", tags$input(id = method_id, type = "number",
          min = "0.1", max = "99.9", step = "0.1", value = "50.0"), span("%")))
    ),
    tags$script(HTML(script))
  )
}

# Shared controls for the two views that estimate eigenvalues.
.app_scree_controls <- function() {
  div(
    class = "scree-method-fields",
    conditionalPanel(
      condition = "input.analysis == 'scree'",
      tags$hr(class = "scree-method-divider")
    ),
    div(
      conditionalPanel(
        condition = "(input.analysis == 'pca' && input.pca_scree_method == 'clipped') || (input.analysis == 'scree' && input.scree_methods && input.scree_methods.indexOf('clipped') >= 0)",
        conditionalPanel(
          condition = "input.analysis == 'scree'",
          div(class = "control-caption scree-method-name", "Clipped")
        ),
        numericInput("C_clip", "C_clip (Clipping threshold)", value = 3, min = 1e-6, step = 0.5)
      ),
      conditionalPanel(
        condition = "(input.analysis == 'pca' && input.pca_scree_method == 'pmwm') || (input.analysis == 'scree' && input.scree_methods && input.scree_methods.indexOf('pmwm') >= 0)",
        conditionalPanel(
          condition = "input.analysis == 'scree'",
          div(class = "control-caption scree-method-name", "PMWM")
        ),
        numericInput("pmwm_beta", "beta (Shape parameter)", value = 1.001, min = 1.000001, step = 0.001),
        div(
          class = "control-pair",
          numericInput("pmwm_a", "a (Lower bound)", value = 0, step = 1),
          numericInput("pmwm_b", "b (Upper bound)", value = 10, step = 1)
        ),
        numericInput("pmwm_trim_const", "trim_const (Trim constant)", value = 10, min = 0, step = 1),
        numericInput("pmwm_eta", "eta (Trimming fraction)", value = 0.01, min = 0, max = 0.49, step = 0.01),
        checkboxInput("pmwm_split_mode", "split_mode (Split mode)", value = TRUE)
      ),
      conditionalPanel(
        condition = "(input.analysis == 'pca' && input.pca_scree_method == 'huber') || (input.analysis == 'scree' && input.scree_methods && input.scree_methods.indexOf('huber') >= 0)",
        conditionalPanel(
          condition = "input.analysis == 'scree'",
          div(class = "control-caption scree-method-name", "Huber")
        ),
        div(
          class = "control-pair",
          numericInput("huber_mu0", "mu0 (Initial mean)", value = 0, step = 0.5),
          numericInput("huber_eta0", "eta0 (Initial scale)", value = 1, min = 1e-8, step = 0.1)
        ),
        div(
          class = "control-pair",
          numericInput("huber_T", "T (Iterations)", value = 50, min = 1, step = 1),
          numericInput("huber_M", "M (Huber parameter)", value = 20, min = 1, step = 1)
        ),
        div(
          class = "control-pair",
          numericInput("huber_k_min_m2", "k_min_m2 (Lower search bound)", value = -20, step = 1),
          numericInput("huber_k_max_m2", "k_max_m2 (Upper search bound)", value = 20, step = 1)
        ),
        numericInput("huber_m2_frac", "m2_frac (Second-moment budget share)", value = 0.25, min = 1e-6, max = 0.99, step = 0.05)
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
      .app-layout { display: grid; grid-template-columns: minmax(0, 2fr) minmax(0, 3fr); gap: 24px; align-items: start; }
      .app-controls { display: grid; grid-template-columns: repeat(2, minmax(0, 1fr)); min-width: 0; margin-bottom: 0; }
      .control-column, .app-results { min-width: 0; }
      .control-column:first-child { padding-right: 20px; }
      .control-column + .control-column { border-left: 1px solid var(--line); padding-left: 20px; }
      .app-controls .shiny-input-container { width: 100%; max-width: 100%; min-width: 0; }
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
      .plot-toolbar { display: flex; align-items: end; flex-wrap: wrap; column-gap: 20px; row-gap: 4px; border-bottom: 1px solid var(--line); padding-bottom: 6px; margin-bottom: 10px; }
      .plot-toolbar > .shiny-input-container { width: 230px; }
      .plot-toolbar .shiny-input-radiogroup { width: auto; }
      .plot-toolbar .shiny-input-radiogroup label { white-space: nowrap; }
      .plot-toolbar > .shiny-panel-conditional { max-width: 100%; }
      .plot-toolbar > .shiny-panel-conditional:last-child { margin-left: auto; }
      .toolbar-alpha { width: 210px; max-width: 100%; }
      .toolbar-alpha .shiny-input-container { width: 100%; }
      .plot-settings-grid { display: grid; grid-template-columns: repeat(3, minmax(0, 1fr)); gap: 5px 22px; }
      .plot-settings-grid .shiny-input-container { max-width: 100%; }
      .plot-area { min-height: 360px; overflow: auto; }
      .pca-plot-area { min-height: 400px; overflow: visible; margin-top: 0; }
      .two-split .pca-shares { grid-template-columns: repeat(2, minmax(0, 1fr)); }
      .two-split .irs--shiny .irs-bar { background: transparent; border: 0; }
      .plot-toolbar { position: relative; z-index: 3; }
      .toolbar-plot-settings {
        position: relative; margin: 0 0 12px; border: 0;
        min-width: 145px;
      }
      .toolbar-plot-settings > summary {
        padding: 9px 13px; border: 1px solid var(--line);
        background: #f6f9fa; border-radius: 7px;
      }
      .toolbar-plot-settings[open] > summary {
        background: #e9f3f3; border-color: #a1c6c2;
      }
      .toolbar-plot-settings .details-body {
        position: absolute; z-index: 30; top: calc(100% + 6px); right: 0;
        width: min(510px, calc(100vw - 65px));
        padding: 18px 18px 10px; border: 1px solid var(--line);
        border-radius: 10px; background: #fff;
        box-shadow: 0 10px 32px rgba(23,45,72,.16);
      }
      .plot-settings-grid { align-items: start; }
      .plot-settings-grid .settings-subheading {
        grid-column: 1 / -1; border-top: 1px solid var(--line);
        padding-top: 9px; color: var(--teal); font-size: 12px;
        font-weight: 700; margin-top: 2px;
      }
      .toolbar-plot-settings .plot-settings-grid { grid-template-columns: repeat(2, minmax(0,1fr)); }
      .toolbar-plot-settings .shiny-input-container { width: 100%; min-width: 0; }
      .toolbar-plot-settings .irs { margin-top: 3px; }
      .pca-plot-area .shiny-plot-output { display: block; }
      .cpp-setting { border-top: 1px solid var(--line); padding-top: 11px;
        margin-top: 8px; }
      .scree-method-fields { border-top: 0; padding-top: 0;
        margin: 0 0 12px; }
      .scree-method-divider { border: 0; border-top: 1px solid var(--line);
        margin: 10px 0 14px; }
      .scree-method-fields .scree-method-name { margin: 12px 0 8px; }
      .scree-method-fields > div { min-width: 0; }
      .scree-plot-area { width: 100%; max-width: 720px; margin: 0 auto; }
      .plot-area .shiny-output-error-validation { padding: 50px 20px; font-size: 14px; color: var(--muted); text-align: center; }
      .table-wrap { overflow-x: auto; }
      .table { font-size: 12px; }
      .table > thead > tr > th { border-bottom: 1px solid var(--line); color: var(--navy); }
      .summary-heading { font-size: 14px; color: var(--navy); font-weight: 650; margin: 23px 0 12px; }
      .fit-info { color: var(--muted); font-size: 12px; margin-top: 18px; }
      .irs--shiny .irs-bar, .irs--shiny .irs-single { background: var(--teal); border-color: var(--teal); }
      @media (max-width: 1600px) {
        .budget-pair { grid-template-columns: 1fr; gap: 0; }
      }
      @media (max-width: 1400px) {
        .control-pair { grid-template-columns: 1fr; gap: 0; }
      }
      @media (max-width: 1100px) {
        .container-fluid { padding-right: 18px; padding-left: 18px; }
        .app-layout { grid-template-columns: 1fr; }
        .well { padding: 16px; }
        .control-pair { grid-template-columns: repeat(2, minmax(0, 1fr)); gap: 12px; }
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
      @media (max-width: 600px) {
        .app-controls { grid-template-columns: 1fr; }
        .control-column:first-child { padding-right: 0; }
        .control-column + .control-column { border-left: 0; border-top: 1px solid var(--line); padding: 20px 0 0; margin-top: 20px; }
      }
      @media (max-width: 420px) {
        .control-pair, .plot-settings-grid { grid-template-columns: 1fr; gap: 0; }
      }
    "))
  ),

  div(
    class = "app-header",
    div(class = "app-wordmark", "dppca"),
    div(
      h1("DP-PCA"),
      p("Estimate and visualize private PC directions, variation, and score distributions.")
    )
  ),

  div(
    class = "app-layout",
    tags$form(
      class = "well app-controls", role = "complementary",
      div(
        class = "control-column",
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
              choices = c(
                "Gaussian" = "gau",
                "Gaussian with groups" = "gau_grouped",
                "Vowel" = "vowel"
              ),
              selected = "gau_grouped"
            )
          ),
          conditionalPanel(
            condition = "input.data_source == 'upload'",
            fileInput("csv_file", "CSV file", accept = ".csv"),
            checkboxInput("csv_header", "First row contains column names", value = TRUE)
          ),
          conditionalPanel(
            condition = "input.analysis == 'pca' || input.analysis == 'score'",
            selectInput(
              "group_col", "Group column",
              choices = c("None" = ""), selected = ""
            ),
            checkboxInput(
              "use_group", "Use group labels for PCA and Score", value = TRUE
            )
          ),
          checkboxInput(
            "include_nonprivate", "Include non-private comparison", value = FALSE
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
                "Set each budget" = "direct"
              ),
              selected = "total"
            )
          ),
          conditionalPanel(
            condition = "input.analysis == 'loading' || input.budget_mode == 'total'",
            div(
              class = "control-pair budget-pair",
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
            condition = "input.analysis == 'scree' && input.budget_mode == 'total'",
            .app_two_split_control("scree")
          ),
          conditionalPanel(
            condition = "input.analysis == 'score' && input.budget_mode == 'total'",
            .app_two_split_control("score")
          ),
          conditionalPanel(
            condition = "input.analysis == 'pca' && input.budget_mode == 'total'",
            .app_pca_split_control()
          ),
          conditionalPanel(
            condition = "input.analysis != 'loading' && input.budget_mode == 'direct'",
            div(class = "control-caption", "Loading"),
            div(
              class = "control-pair budget-pair",
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
                class = "control-pair budget-pair",
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
                class = "control-pair budget-pair",
                numericInput(
                  "budget_score_eps", "Epsilon",
                  value = 10 / 3, min = 0, step = 0.1
                ),
                numericInput(
                  "budget_score_delta", "Delta",
                  value = 1e-6, min = 0, max = 0.99, step = 1e-6
                )
              )
            )
          ),
          div(class = "table-wrap", tableOutput("budget_preview")),
          div(class = "cpp-setting",
            checkboxInput("cpp", "Use C++ implementation", value = TRUE)
          )
        )
      ),
      div(
        class = "control-column",
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
              checkboxInput("score_center", "Center DP score display", value = TRUE),
              numericInput(
                "sampling_size", "Sampled points",
                value = 500, min = 1, step = 100
              ),
              helpText("Centering only changes the display. Sample size takes effect on the next run; histogram and sample views reuse the same estimate.")
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
            conditionalPanel(
              condition = "input.loading_display && input.loading_display.indexOf('heatmap') >= 0",
              selectizeInput(
                "heatmap_components", "Heatmap PCs",
                choices = NULL, multiple = TRUE
              ),
              NULL
            )
          )
        ),
        div(class = "run-area", actionButton("run", "Run PCA", class = "btn-primary"))
      )
    ),

    tags$main(
      class = "app-results",
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
              uiOutput("pca_reference_controls"),
              conditionalPanel(
                condition = "input.pca_view != 'summary'",
                tags$details(
                  class = "app-details toolbar-plot-settings",
                  tags$summary("Plot settings"),
                  div(
                    class = "details-body plot-settings-grid",
                    numericInput("pca_variables", "Top variables", value = 10, min = 1, step = 1),
                    selectInput("pca_loading_display", "Loading display", choices = c("Vector" = "vector", "Point" = "point", "Heatmap" = "heatmap"), selected = "vector"),
                    selectInput("pca_score_display", "Score display", choices = c("Histogram" = "histogram", "Sampled points" = "sample"), selected = "histogram"),
                    checkboxInput("pca_center", "Center DP score and biplot displays", value = TRUE),
                    div(class = "settings-subheading", "Biplot"),
                    conditionalPanel(
                      condition = "input.pca_view == 'all' || input.pca_view == 'biplot'",
                      sliderInput("alpha", "Biplot alpha", min = 0, max = 1,
                                  value = 0, step = 0.01, ticks = FALSE, width = "100%")
                    ),
                    checkboxInput("pca_labels", "Show variable labels", value = TRUE),
                    div(class = "settings-subheading", "Synthetic sampling"),
                    checkboxInput("pca_custom_sample", "Custom sample size", value = FALSE),
                    conditionalPanel(
                      condition = "input.pca_custom_sample",
                      numericInput("pca_sampling_size", "Sample points (total)",
                        value = 1000, min = 1, step = 100)
                    )
                  )
                )
              )
            ),
            conditionalPanel(
              condition = "input.pca_view != 'summary'",
              div(class = "plot-area pca-plot-area", plotOutput("pca_plot", height = "auto"))
            ),
            conditionalPanel(
              condition = "input.pca_view == 'summary'",
              h3(class = "summary-heading", "Component importance"),
              div(class = "table-wrap", tableOutput("importance")),
              h3(class = "summary-heading", "Privacy allocation"),
              div(class = "table-wrap", tableOutput("pca_budget"))
            )
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
  preview_variables <- reactiveVal(NULL)
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
    label <- switch(
      input$example_dataset,
      gau = "Gaussian",
      vowel = "Vowel"
    )
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
    group_arg <- NULL
    selected <- input$group_col
    if (is.null(selected)) selected <- ""
    features <- data$data

    if (identical(selected, ".provided_group")) {
      group <- data$external_group
      group_arg <- group
    } else if (nzchar(selected)) {
      validate(need(selected %in% names(features), "Choose a group column."))
      group <- features[[selected]]
      group_arg <- selected
      if (isTRUE(input$use_group) && input$analysis %in% c("pca", "score")) {
        features <- features[, names(features) != selected, drop = FALSE]
      }
    }

    list(
      X = data$data,
      numeric_X = .app_matrix(features),
      group = group,
      group_arg = group_arg,
      label = data$label
    )
  }))

  observeEvent(prepared_data(), {
    data <- prepared_data()
    p <- ncol(data$numeric_X)
    n <- nrow(data$numeric_X)
    if (!identical(colnames(data$numeric_X), preview_variables())) {
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
      preview_variables(colnames(data$numeric_X))
    }
    if (!identical(n, delta_sample_size())) {
      updateNumericInput(session, "pca_sampling_size", value = as.integer(n))
      updateNumericInput(session, "sampling_size", value = as.integer(n))
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
      direct_budgets[[mode]] <- privacy$budgets
      show_direct_budget(privacy$budgets)
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
    warning_messages <- character()

    computed <- tryCatch(withCallingHandlers({
      if (identical(input$data_source, "upload") && is.null(input$csv_file)) {
        stop("Upload a CSV file before running the analysis.", call. = FALSE)
      }
      data <- prepared_data()
      X <- data$X
      p <- ncol(data$numeric_X)
      group_arg <- if (isTRUE(input$use_group)) data$group_arg else NULL
      privacy <- if (mode == "loading") NULL else .app_privacy(input, mode)
      axes <- NULL
      bins <- NULL
      k <- NULL
      plot_args <- list()

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
            cpp.option = isTRUE(input$cpp),
            non_private = isTRUE(input$include_nonprivate),
            group = group_arg
          )
        } else if (mode == "loading") {
          components <- as.integer(input$heatmap_components)
          if (!length(components)) components <- axes
          plot_args <- list(
            axes = axes,
            type = c("vector", "point", "heatmap"),
            heatmap_components = components,
            variables = .app_integer(input$variables, "Variables", upper = p)
          )
          value <- dp_loading(
            X = X,
            privacy = privacy_control(
              loading = c(eps = input$eps, delta = input$delta)
            ),
            cpp.option = isTRUE(input$cpp),
            non_private = isTRUE(input$include_nonprivate)
          )
        } else if (mode == "scree") {
          methods <- unique(input$scree_methods)
          value <- dp_scree(
            X = X,
            k = k,
            privacy = privacy,
            method = methods,
            control = .app_controls(input, methods),
            cpp.option = isTRUE(input$cpp),
            mono = isTRUE(input$mono),
            non_private = isTRUE(input$include_nonprivate)
          )
        } else {
          methods <- unique(input$score_methods)
          if (!length(methods)) stop("Choose at least one score method.", call. = FALSE)
          plot_args <- list(
            type = c("histogram", "sample"),
            plot_control = score_plot_control(show_group_legend = TRUE),
            sampling_control = sampling_control(
              method = "uniform",
              sample_size = .app_integer(input$sampling_size, "Sample points")
            )
          )
          args <- list(
            X = X,
            privacy = privacy,
            bins = bins,
            method = methods,
            cpp.option = isTRUE(input$cpp),
            axes = axes,
            non_private = isTRUE(input$include_nonprivate),
            group = group_arg
          )
          value <- do.call(dp_score, args)
        }
        incProgress(0.9, detail = "Saving results for display")
        value
      })

      list(
        value = value,
        plot_args = plot_args,
        seed = sample.int(.Machine$integer.max, 1L),
        label = data$label,
        n = nrow(X),
        p = p,
        warnings = warning_messages,
        privacy = privacy
      )
    }, warning = function(w) {
      warning_messages <<- unique(c(warning_messages, conditionMessage(w)))
      invokeRestart("muffleWarning")
    }), error = function(e) {
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

  # Only offer comparison controls when the completed fit contains a reference.
  output$pca_reference_controls <- renderUI({
    result <- results$pca
    if (is.null(result) || is.null(result$value$nonprivate)) return(NULL)
    if (identical(input$pca_view, "all")) {
      radioButtons("pca_version", "Overview results",
        choices = c("DP" = "private", "Non-private" = "nonprivate"),
        selected = "private", inline = TRUE)
    } else if (!identical(input$pca_view, "summary")) {
      checkboxInput("pca_compare", "Compare with non-private results", value = TRUE)
    }
  })

  # NULL preserves the package default: n private synthetic points and
  # observed (rather than resampled) non-private comparison points.
  pca_sampling <- reactive({
    if (!isTRUE(input$pca_custom_sample)) return(NULL)
    sampling_control(
      method = "uniform",
      sample_size = .app_integer(input$pca_sampling_size, "Sample points")
    )
  })

  # These reactives depend only on saved runs, not the current form inputs.
  # View changes select panels from the same S3 plots and synthetic samples.
  loading_plots <- reactive({
    result <- results$loading
    req(result)
    .app_stored_plots(result)
  })

  score_plots <- reactive({
    result <- results$score
    req(result)
    .app_stored_plots(result, center = isTRUE(input$score_center))
  })

  # Same stored seed = same synthetic points; appearance changes spend no
  # additional estimation budget and do not change the fitted quantities.
  pca_base_plots <- reactive({
    result <- results$pca
    req(result)
    variables <- .app_integer(input$pca_variables, "Variables", upper = result$p)
    .app_pca_plots(
      fit = result$value,
      seed = result$seed,
      variables = variables,
      biplot_labels = isTRUE(input$pca_labels),
      loading_display = input$pca_loading_display,
      score_display = input$pca_score_display,
      center = isTRUE(input$pca_center),
      sampling = pca_sampling()
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
      biplot_labels = isTRUE(input$pca_labels),
      center = isTRUE(input$pca_center),
      sampling = pca_sampling()
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
        for (version in intersect(names(plots), names(biplots))) {
          plots[[version]]$plots$biplot <- biplots[[version]]$plots$biplot
        }
      }
    }
    .app_draw_pca(
      plots,
      view = input$pca_view,
      version = if (input$pca_view == "all" &&
        !is.null(input$pca_version)) input$pca_version else "private",
      compare = isTRUE(input$pca_compare)
    )
  }, height = function() {
    width <- session$clientData$output_pca_plot_width
    if (length(width) != 1L || !is.finite(width) || width <= 0) width <- 900
    if (identical(input$pca_view, "all")) {
      max(420, min(650, round(width * 0.65)))
    } else {
      max(360, min(570, round(width * 0.57)))
    }
  }, res = 110)


  pca_summary <- reactive({
    req(results$pca)
    summary(results$pca$value)
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
    .app_draw_loading(loading_plots(), display = input$loading_display)
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
    .app_draw_score(score_plots(), display = input$score_display)
  }, height = function() {
    fit <- results$score
    n_methods <- if (!is.null(fit)) length(fit$value$metadata$methods) else
      length(input$score_methods)
    has_reference <- !is.null(fit) && !is.null(fit$value$nonprivate)
    n_displays <- max(1L, length(input$score_display))
    cols <- max(1L, n_methods + as.integer(has_reference))
    single_pair <- cols == 1L &&
      identical(sort(input$score_display), c("histogram", "sample"))
    width <- session$clientData$output_score_plot_width
    if (length(width) != 1L || !is.finite(width) || width <= 0) width <- 900
    plot_cols <- if (single_pair) 2L else cols
    plot_rows <- if (single_pair) 1L else n_displays
    per_panel <- max(280, min(600, round(width / plot_cols * 1.15)))
    as.integer(plot_rows * per_panel)
  }, res = 110)
}

shinyApp(ui, server)
