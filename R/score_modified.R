# `code/MW-score_plot_modified.R` 파일안의 dp_score_ver02 함수가 스크립트 내부에 작성된 Internel helper functions를 기준으로 잘 동작하는지 체크해줘. 간단한 시뮬레이션을 돌려도 좋아. 이 과정에서 맞거나 옳게 수정해야하는 부분이 있다면 수정해주고, 대신에 기존 버전의 파일도 log에 저장해줘. 해당 파일을 수정할 때는 너 스타일로 코드를 너무 바꾸거나 수정하지 말고, 지금 있는 스타일과 구조를 그대로 살리면서 오류가 될 부분만 수정해줘.

# Ver 02 DP score implementation
#
# Minwoo Kim, 2026/08/10(Mon) 
#
# Main changes
# - Durfee's DP quantile estimates are used instead of smooth sensitivity approach
# - privacy budget allocation is changed
# - dp_score provides only one of "add" and "sparse" histogram, not both.

dp_score_ver02 <- function(
    X,
    eps,
    delta,
    bins,
    method = c("add", "sparse"),
    center = TRUE,
    standardize = FALSE,
    g_dppca = FALSE,
    cpp.option = FALSE,
    axes = c(1, 2),
    estimate_center = TRUE
) {
  X <- validate_score_matrix(X)
  validate_score_common(
    X = X,
    eps = eps,
    delta = delta,
    bins = bins,
    center = center,
    standardize = standardize,
    g_dppca = g_dppca,
    cpp.option = cpp.option,
    axes = axes,
    estimate_center = estimate_center
  )
  
  axes <- as.integer(axes)
  bins <- as.integer(bins)
  m_x <- bins[1]
  m_y <- bins[2]
  
  method <- match.arg(method)
  
  budget <- split_score_privacy_budget_ver02(
    eps = eps,
    delta = delta,
    g_dppca = g_dppca
  )
  
  score_res <- compute_score_coordinates(
    X = X,
    axes = axes,
    center = center,
    standardize = standardize,
    g_dppca = g_dppca,
    cpp.option = cpp.option,
    eps_pc = budget$eps_pc,
    delta_pc = budget$delta_pc
  )
  
  frame_out <- dp_frame_ver02(
    X = score_res$score,
    eps_frame = budget$eps_frame,
    inflate = 0.20,
    center = center,
    estimate_center = estimate_center
  )
  
  hist <- score_histograms_ver02(
    X_score = score_res$score,
    xlim = frame_out$xlim,
    ylim = frame_out$ylim,
    bins = bins,
    eps_hist = budget$eps_hist,
    delta_hist = budget$delta_hist,
    method = method
  )
  
  list(
    score = score_res$score,
    frame = frame_out,
    hist_nonpriv = hist$non_dp,
    hist_dp = hist$dp,
    method = method
  )
}

# Description:
# - make plots that indicates by `method` argument
dp_score_plot_ver02 <- function(
    X,
    eps,
    delta,
    bins,
    method = c("add", "sparse"),
    center = TRUE,
    standardize = FALSE,
    g_dppca = FALSE,
    cpp.option = FALSE,
    axes = c(1, 2),
    estimate_center = TRUE,
    opt.kernel = FALSE  # also provide kernel plot if TRUE / under-developed...
) {
  method <- match.arg(method, choices = c("add", "sparse"), several.ok = TRUE)
  color <- "#6A5ACD"
  
  p_add <- NULL
  p_sparse <- NULL
  p_none <- NULL 
  p_scatter <- NULL
  
  if("add" %in% method){
    # calculate score and private histograms
    score_res <- dp_score_ver02(
      X = X,
      eps = eps,
      delta = delta,
      bins = bins,
      center = center,
      standardize = standardize,
      g_dppca = g_dppca,
      cpp.option = cpp.option,
      axes = axes,
      method = "add",
      estimate_center = estimate_center
    )
    xlim <- score_res$frame$xlim
    ylim <- score_res$frame$ylim
    pc_names <- colnames(score_res$score)
    
    p_add <- dppca:::make_hist_plot_dp(score_res$hist_dp, xlim, ylim, color, "Add DP Hist",
                                       xlab = pc_names[1], ylab = pc_names[2])
  }
  
  if("sparse" %in% method){
    # calculate score and private histograms
    score_res <- dp_score_ver02(
      X = X,
      eps = eps,
      delta = delta,
      bins = bins,
      center = center,
      standardize = standardize,
      g_dppca = g_dppca,
      cpp.option = cpp.option,
      axes = axes,
      method = "sparse",
      estimate_center = estimate_center
    )
    xlim <- score_res$frame$xlim
    ylim <- score_res$frame$ylim
    pc_names <- colnames(score_res$score)
    
    p_sparse <- dppca:::make_hist_plot_dp(score_res$hist_dp, xlim, ylim, color, "Sparse DP Hist",
                                          xlab = pc_names[1], ylab = pc_names[2])
  }
  
  # get non-private scatter plot
  X_score <- as.data.frame(score_res$score)
  colnames(X_score) <- c("pc_x", "pc_y")
  p_scatter <- ggplot2::ggplot(
    X_score,
    ggplot2::aes(x = .data$pc_x, y = .data$pc_y)
  ) +
    ggplot2::geom_point(alpha = 0.6, size = 1.8, color = color) +
    ggplot2::coord_fixed(xlim = xlim, ylim = ylim) +
    ggplot2::scale_x_continuous(expand = c(0, 0), breaks = pretty(xlim, n = 5)) +
    ggplot2::scale_y_continuous(expand = c(0, 0), breaks = pretty(ylim, n = 5)) +
    dppca:::theme_dp_base() +
    ggplot2::labs(x = pc_names[1], y = pc_names[2])
  p_scatter <- dppca:::add_title_dp(p_scatter, "Original Scatter")
  
  # get non-private histogram plot
  p_none <- dppca:::make_hist_plot_dp(score_res$hist_nonpriv, xlim, ylim, color, "Original Hist",
                                      xlab = pc_names[1], ylab = pc_names[2])
  
  # Aggregated plots
  plot_panels <- list(p_scatter, p_none, p_add, p_sparse)
  plot_panels <- Filter(Negate(is.null), plot_panels)
  p_all <- patchwork::wrap_plots(plot_panels, nrow = 2)
  
  # Final output
  list(
    score = score_res,
    plot = list(
      scatter = p_scatter,
      none = p_none,
      add = p_add,
      sparse = p_sparse,
      all = p_all
    )
  )
}



# Internal helpers ------------------------------------------------------------
# some functions are modified from the original package implementations
# the Durfee quantile functions are newly added

split_score_privacy_budget_ver02 <- function(eps, delta, g_dppca) {
  # unequal splitting of privacy budget with more weight on histogram and Durfee's approach
  if (isTRUE(g_dppca)) {
    list(
      eps_pc = 0.2 * eps,
      eps_frame = 0.2 * eps,
      eps_hist = 0.6 * eps,
      delta_pc = 0.2 * delta,
      delta_frame = NULL,
      delta_hist = 0.8 * delta
    )
  } else {
    list(
      eps_pc = NULL,
      eps_frame = 0.2 * eps,
      eps_hist = 0.8 * eps,
      delta_pc = NULL,
      delta_frame = NULL,
      delta_hist = delta
    )
  }
}

#' original version as in the dppca package implementation
compute_score_coordinates <- function(
  X,
  axes,
  center,
  standardize,
  g_dppca,
  cpp.option,
  eps_pc,
  delta_pc
) {
  k_max <- max(axes)
  
  X_proc <- dppca:::prep_matrix_for_pca(
    X = X,
    center = center,
    standardize = standardize
  )
  
  V_all <- dppca::dp_pc_dir(
    X = X,
    k = k_max,
    center = center,
    standardize = standardize,
    g_dppca = g_dppca,
    eps = eps_pc,
    delta = delta_pc,
    cpp.option = cpp.option
  )
  
  V <- V_all[, axes, drop = FALSE]
  X_score <- as.matrix(X_proc %*% V)
  colnames(X_score) <- paste0("PC", axes)
  
  list(score = X_score, directions = V)
}




# modified version...
# estimates radius for each coordinate
dp_frame_ver02 <- function(
    X,
    eps_frame,
    inflate = 0.20,
    center = TRUE,
    estimate_center = TRUE
) {
  X <- as.matrix(X)
  
  if (!is.numeric(X) || ncol(X) != 2L) {
    stop("`X` must be a numeric matrix with exactly two columns.", call. = FALSE)
  }
  if (nrow(X) < 2L) {
    stop("`X` must have at least two rows.", call. = FALSE)
  }
  if (anyNA(X) || any(!is.finite(X))) {
    stop("`X` must contain only finite values.", call. = FALSE)
  }
  if (
    !is.numeric(eps_frame) || length(eps_frame) != 1L ||
    !is.finite(eps_frame) || eps_frame <= 0
  ) {
    stop("`eps_frame` must be a positive number.", call. = FALSE)
  }
  if (
    !is.numeric(inflate) || length(inflate) != 1L ||
    !is.finite(inflate) || inflate < 0
  ) {
    stop("`inflate` must be a nonnegative number.", call. = FALSE)
  }
  validate_logical_value(center, "center")
  validate_logical_value(estimate_center, "estimate_center")
  if (!estimate_center && !center) {
    stop(
      "`estimate_center = FALSE` requires `center = TRUE`.",
      call. = FALSE
    )
  }
  
  if (estimate_center) {
    eps_each <- eps_frame / 4
    
    center_x <- dp_quantile_durfee_signed(
      X[, 1], q = 0.5, epsilon = eps_each
    )$value
    
    center_y <- dp_quantile_durfee_signed(
      X[, 2], q = 0.5, epsilon = eps_each
    )$value
    
  } else {
    eps_each <- eps_frame / 2
    
    center_x <- 0
    center_y <- 0
  }
  
  radius_values_x <- abs(X[, 1] - center_x)
  r_x_max <- dp_quantile_durfee(
    radius_values_x, q = 0.99, epsilon = eps_each
  )$value
  
  radius_values_y <- abs(X[, 2] - center_y)
  r_y_max <- dp_quantile_durfee(
    radius_values_y, q = 0.99, epsilon = eps_each
  )$value
  
  if ((!is.finite(r_x_max) || r_x_max <= 0) || (!is.finite(r_y_max) || r_y_max <= 0)){
    stop(
      "The private frame radius must be positive and finite. ",
      "Try a larger privacy budget.",
      call. = FALSE
    )
  }
  
  inflated_r_x_max <- (1 + inflate) * r_x_max
  inflated_r_y_max <- (1 + inflate) * r_y_max
  
  list(
    xlim = c(center_x - inflated_r_x_max, center_x + inflated_r_x_max),
    ylim = c(center_y - inflated_r_y_max, center_y + inflated_r_y_max)
  )
}

# using Durfee's DP quantile estimation. Here, only epsilon budget is used.
dp_frame_ver02_old <- function(
    X,
    eps_frame,
    inflate = 0.20,
    center = TRUE,
    estimate_center = TRUE
) {
  X <- as.matrix(X)
  
  if (!is.numeric(X) || ncol(X) != 2L) {
    stop("`X` must be a numeric matrix with exactly two columns.", call. = FALSE)
  }
  if (nrow(X) < 2L) {
    stop("`X` must have at least two rows.", call. = FALSE)
  }
  if (anyNA(X) || any(!is.finite(X))) {
    stop("`X` must contain only finite values.", call. = FALSE)
  }
  if (
    !is.numeric(eps_frame) || length(eps_frame) != 1L ||
    !is.finite(eps_frame) || eps_frame <= 0
  ) {
    stop("`eps_frame` must be a positive number.", call. = FALSE)
  }
  if (
    !is.numeric(inflate) || length(inflate) != 1L ||
    !is.finite(inflate) || inflate < 0
  ) {
    stop("`inflate` must be a nonnegative number.", call. = FALSE)
  }
  validate_logical_value(center, "center")
  validate_logical_value(estimate_center, "estimate_center")
  if (!estimate_center && !center) {
    stop(
      "`estimate_center = FALSE` requires `center = TRUE`.",
      call. = FALSE
    )
  }
  
  q_radius <- 0.99
  
  if (estimate_center) {
    eps_each <- eps_frame / 3
    center_x <- dp_quantile_durfee_signed(
      X[, 1], q = 0.5, epsilon = eps_each
    )$value
    center_y <- dp_quantile_durfee_signed(
      X[, 2], q = 0.5, epsilon = eps_each
    )$value
  } else {
    eps_each <- eps_frame
    center_x <- 0
    center_y <- 0
  }
  
  radius_values <- sqrt((X[, 1] - center_x)^2 + (X[, 2] - center_y)^2)
  radius_out <- dp_quantile_durfee(
    radius_values, q = q_radius, epsilon = eps_each
  )
  radius <- radius_out$value
  
  if (!is.finite(radius) || radius <= 0) {
    stop(
      "The private frame radius is not positive. ",
      "Try a larger privacy budget.",
      call. = FALSE
    )
  }
  
  inflated_radius <- (1 + inflate) * radius
  
  list(
    xlim = c(center_x - inflated_radius, center_x + inflated_radius),
    ylim = c(center_y - inflated_radius, center_y + inflated_radius)
  )
}

score_histogram_grid <- function(xlim, ylim, m_x, m_y) {
  x_breaks <- seq(xlim[1], xlim[2], length.out = m_x + 1L)
  y_breaks <- seq(ylim[1], ylim[2], length.out = m_y + 1L)
  
  base_coord <- do.call(
    rbind,
    lapply(seq_len(m_y), function(j) {
      data.frame(
        xmin = x_breaks[-length(x_breaks)],
        xmax = x_breaks[-1],
        ymin = y_breaks[j],
        ymax = y_breaks[j + 1L]
      )
    })
  )
  
  list(
    x_breaks = x_breaks,
    y_breaks = y_breaks,
    base_coord = base_coord,
    m_x = m_x,
    m_y = m_y,
    m = m_x * m_y
  )
}


# Modified
# - privacy budget splitting; output only one version
score_histograms_ver02 <- function(
    X_score,
    xlim,
    ylim,
    bins,
    eps_hist,
    delta_hist,
    method=c("add", "sparse")
) {
  method <- match.arg(method)
  
  validate_bins(bins)
  bins <- as.integer(bins)
  grid <- score_histogram_grid(
    xlim = xlim,
    ylim = ylim,
    m_x = bins[1],
    m_y = bins[2]
  )
  
  score_histograms_from_grid_ver02(
    X_score = X_score,
    grid = grid,
    eps_hist = eps_hist,
    delta_hist = delta_hist,
    method = method,
    group_name = NULL
  )
}

score_histograms_from_grid_ver02 <- function(
    X_score,
    grid,
    eps_hist,
    delta_hist,
    method = c("add", "sparse"),  # add or sparse
    group_name = NULL
) {
  method <- match.arg(method)
  
  n <- nrow(X_score)
  if (n < 1L) {
    stop("Each histogram must contain at least one observation.", call. = FALSE)
  }
  
  bx <- cut(
    X_score[, 1],
    breaks = grid$x_breaks,
    include.lowest = TRUE,
    labels = FALSE
  )
  by <- cut(
    X_score[, 2],
    breaks = grid$y_breaks,
    include.lowest = TRUE,
    labels = FALSE
  )
  
  if (anyNA(bx)) {
    idx <- which(is.na(bx))
    bx[idx] <- findInterval(X_score[idx, 1], grid$x_breaks, all.inside = TRUE)
  }
  if (anyNA(by)) {
    idx <- which(is.na(by))
    by[idx] <- findInterval(X_score[idx, 2], grid$y_breaks, all.inside = TRUE)
  }
  
  bidx <- (by - 1L) * grid$m_x + bx
  counts <- as.numeric(table(factor(bidx, levels = seq_len(grid$m))))
  p_hat <- counts / n
  
  hist_none <- grid$base_coord
  hist_none$prob <- p_hat
  
  hist <- NULL
  
  if (method == "add") {
    sigma <- sqrt(2) * sqrt(2 * log(1.25 / delta_hist)) / eps_hist
    c_tilde <- pmax(counts + stats::rnorm(grid$m, mean = 0, sd = sigma), 0)
    
    if (sum(c_tilde) <= 0) {
      prefix <- if (is.null(group_name)) "" else paste0("Group `", group_name, "`: ")
      stop(
        prefix,
        "all privatized bin counts are zero after Gaussian noise and clipping. ",
        "Try a larger `eps` or fewer bins.",
        call. = FALSE
      )
    }
    
    hist <- grid$base_coord
    hist$prob <- c_tilde / sum(c_tilde)
    
  } else if (method == "sparse") {
    q_sparse <- numeric(grid$m)
    scale_lap <- 2 / (eps_hist * n)
    threshold <- (2 * log(2 / delta_hist)) / (eps_hist * n) + 1 / n
    
    for (kk in seq_len(grid$m)) {
      if (p_hat[kk] == 0) {
        q_sparse[kk] <- 0
      } else {
        z_k <- VGAM::rlaplace(1, location = 0, scale = scale_lap)
        q_k <- max(p_hat[kk] + z_k, 0)
        q_sparse[kk] <- if (q_k < threshold) 0 else q_k
      }
    }
    
    if (sum(q_sparse) > 0) {
      q_sparse <- q_sparse / sum(q_sparse)
    }
    
    hist <- grid$base_coord
    hist$prob <- q_sparse
  }
  
  list(
    non_dp = hist_none,
    dp = hist
  )
}

# the following function is written by ChatGPT
dp_quantile_durfee <- function(x, q, epsilon, lower = 0,
                               beta = 1.01, max_steps = 5000L) {
  x <- as.numeric(x)
  if (length(x) < 1L || anyNA(x) || any(!is.finite(x))) {
    stop("Durfee quantile input must contain finite values.", call. = FALSE)
  }
  if (any(x < lower)) {
    stop("Durfee's one-sided estimator requires a valid public lower bound.",
         call. = FALSE)
  }
  if (!is.numeric(q) || length(q) != 1L || q <= 0 || q >= 1) {
    stop("`q` must lie in (0, 1).", call. = FALSE)
  }
  if (!is.numeric(epsilon) || length(epsilon) != 1L || epsilon <= 0) {
    stop("Durfee quantile epsilon must be positive.", call. = FALSE)
  }
  if (!is.numeric(beta) || length(beta) != 1L || beta <= 1) {
    stop("`durfee_beta` must be greater than 1.", call. = FALSE)
  }
  max_steps <- as.integer(max_steps)
  if (!is.finite(max_steps) || max_steps < 1L) {
    stop("`durfee_max_steps` must be a positive integer.", call. = FALSE)
  }
  
  # Monotone AboveThreshold with sensitivity one. Exponential noise and
  # epsilon_1 = epsilon_2 = epsilon / 2 give pure epsilon-DP under
  # fixed-size replacement adjacency.
  eps_threshold <- epsilon / 2
  eps_query <- epsilon / 2
  noisy_threshold <- q * length(x) + stats::rexp(1L, rate = eps_threshold)
  sorted_x <- sort(x)
  estimate <- lower
  halted <- FALSE
  
  for (i in seq_len(max_steps)) {
    estimate <- beta^i + lower - 1
    if (!is.finite(estimate)) {
      estimate <- .Machine$double.xmax
    }
    count_below <- findInterval(estimate, sorted_x, left.open = TRUE)
    noisy_count <- count_below + stats::rexp(1L, rate = eps_query)
    if (noisy_count >= noisy_threshold) {
      halted <- TRUE
      break
    }
  }
  if (!halted) {
    warning(
      "Durfee quantile search reached `durfee_max_steps`; consider a smaller `durfee_beta` only with a larger step limit, or increase the step limit.",
      call. = FALSE
    )
  }
  
  list(value = estimate, halted = halted, steps = i)
}

# fully unbounded extension in Algorithm 4 of Durfee (2023)
dp_quantile_durfee_signed <- function(x, q, epsilon,
                                      beta = 1.01, max_steps = 5000L) {
  x <- as.numeric(x)
  if (length(x) < 1L || anyNA(x) || any(!is.finite(x))) {
    stop("Durfee quantile input must contain finite values.", call. = FALSE)
  }
  if (
    !is.numeric(q) || length(q) != 1L ||
    !is.finite(q) || q <= 0 || q >= 1
  ) {
    stop("`q` must lie in (0, 1).", call. = FALSE)
  }
  if (
    !is.numeric(epsilon) || length(epsilon) != 1L ||
    !is.finite(epsilon) || epsilon <= 0
  ) {
    stop("Durfee quantile epsilon must be positive.", call. = FALSE)
  }
  if (
    !is.numeric(beta) || length(beta) != 1L ||
    !is.finite(beta) || beta <= 1
  ) {
    stop("`durfee_beta` must be greater than 1.", call. = FALSE)
  }
  max_steps <- as.integer(max_steps)
  if (!is.finite(max_steps) || max_steps < 1L) {
    stop("`durfee_max_steps` must be a positive integer.", call. = FALSE)
  }
  
  # Algorithm 4 uses two independent AboveThreshold calls. Each call
  # receives half of epsilon and splits it equally between threshold
  # and query noise.
  run_search <- function(values, quantile) {
    eps_threshold <- epsilon / 4
    eps_query <- epsilon / 4
    noisy_threshold <- quantile * length(values) +
      stats::rexp(1L, rate = eps_threshold)
    sorted_values <- sort(values)
    estimate <- 0
    halted <- FALSE
    step <- 0L
    
    for (i in seq.int(0L, max_steps - 1L)) {
      estimate <- beta^i - 1
      if (!is.finite(estimate)) {
        estimate <- .Machine$double.xmax
      }
      count_below <- findInterval(
        estimate, sorted_values, left.open = TRUE
      )
      noisy_count <- count_below + stats::rexp(1L, rate = eps_query)
      step <- i
      if (noisy_count >= noisy_threshold) {
        halted <- TRUE
        break
      }
    }
    
    list(value = estimate, halted = halted, steps = step)
  }
  
  positive_out <- run_search(x, q)
  negative_out <- run_search(-x, 1 - q)
  
  estimate <- if (positive_out$halted && positive_out$steps > 0L) {
    positive_out$value
  } else if (negative_out$halted && negative_out$steps > 0L) {
    -negative_out$value
  } else {
    0
  }
  
  if (!positive_out$halted || !negative_out$halted) {
    warning(
      "Durfee signed quantile search reached `durfee_max_steps`; consider a smaller `durfee_beta` only with a larger step limit, or increase the step limit.",
      call. = FALSE
    )
  }
  
  list(
    value = estimate,
    halted = positive_out$halted && negative_out$halted,
    steps = c(
      positive = positive_out$steps,
      negative = negative_out$steps
    )
  )
}

validate_score_matrix <- function(X) {
  X <- as.matrix(X)
  
  if (!is.numeric(X)) {
    stop("`X` must be numeric or coercible to a numeric matrix.", call. = FALSE)
  }
  if (nrow(X) < 2L) {
    stop("`X` must have at least two rows.", call. = FALSE)
  }
  if (ncol(X) < 2L) {
    stop("`X` must have at least two columns.", call. = FALSE)
  }
  if (anyNA(X) || any(!is.finite(X))) {
    stop("`X` must contain only finite values.", call. = FALSE)
  }
  
  X
}

validate_score_common <- function(
    X,
    eps,
    delta,
    bins,
    center,
    standardize,
    g_dppca,
    cpp.option,
    axes,
    estimate_center
) {
  validate_logical_value(center, "center")
  validate_logical_value(standardize, "standardize")
  validate_logical_value(g_dppca, "g_dppca")
  validate_logical_value(cpp.option, "cpp.option")
  validate_logical_value(estimate_center, "estimate_center")
  if (!estimate_center && !center) {
    stop(
      "`estimate_center = FALSE` requires `center = TRUE`.",
      call. = FALSE
    )
  }
  
  if (!is.numeric(eps) || length(eps) != 1L || !is.finite(eps) || eps <= 0) {
    stop("`eps` must be a positive number.", call. = FALSE)
  }
  if (
    !is.numeric(delta) || length(delta) != 1L ||
    !is.finite(delta) || delta <= 0 || delta >= 1
  ) {
    stop("`delta` must be a number in `(0, 1)`.", call. = FALSE)
  }
  
  validate_bins(bins)
  
  if (length(axes) != 2L || !is.numeric(axes)) {
    stop("`axes` must be an integer vector of length 2.", call. = FALSE)
  }
  if (anyNA(axes) || any(!is.finite(axes)) || any(axes <= 0) ||
      any(axes != as.integer(axes))) {
    stop("`axes` must contain positive integers.", call. = FALSE)
  }
  if (max(axes) > ncol(X)) {
    stop("The largest value in `axes` must not exceed `ncol(X)`.", call. = FALSE)
  }
  
  invisible(TRUE)
}

validate_logical_value <- function(x, arg) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop("`", arg, "` must be `TRUE` or `FALSE`.", call. = FALSE)
  }
  
  invisible(TRUE)
}

validate_positive_integer <- function(x, arg) {
  if (
    !is.numeric(x) || length(x) != 1L || !is.finite(x) ||
    x < 1 || x != as.integer(x)
  ) {
    stop("`", arg, "` must be a positive integer.", call. = FALSE)
  }
  
  invisible(TRUE)
}

validate_bins <- function(bins) {
  if (
    !is.numeric(bins) || length(bins) != 2L || anyNA(bins) ||
    any(!is.finite(bins)) || any(bins < 1) || any(bins != as.integer(bins))
  ) {
    stop("`bins` must be an integer vector of length 2 with positive values.", call. = FALSE)
  }
  
  invisible(TRUE)
}
