#' Plot differentially private scree curves
#'
#' User-facing function that computes and plots nonprivate and/or differentially
#' private scree curves using one or more DP eigenvalue estimators.
#'
#' @param X Numeric data matrix with observations in rows.
#' @param k Number of leading components.
#' @param dp_lambda Which DP estimator(s) to plot. One of
#'   \code{"both"}, \code{"clipped"}, \code{"pmwm"}, or \code{"huber"}.
#' @param eps_total Total privacy budget epsilon.
#' @param delta_total Total privacy budget delta.
#' @param center Logical; whether to center columns before PCA.
#' @param standardize Logical; whether to standardize columns before PCA.
#' @param B_mult Clipping multiplier for the clipped estimator.
#' @param beta Log-binning base for PMWM.
#' @param a,b Support bounds for PMWM.
#' @param trim_const,eta PMWM trimming parameters.
#' @param split_mode Logical; whether PMWM splits the sample.
#' @param mu0,eta0,T,M Huber estimator parameters.
#' @param k_min_star,k_max_star,include_zero_bin Huber scale-estimation parameters.
#' @param m2_frac Huber privacy split parameter.
#' @param step_schedule Huber GD step schedule.
#' @param dp_v_flag Logical; whether to privatize the PCA basis.
#' @param cpp.option Logical passed to \code{mech_tau_sph()}.
#' @param mono Logical; whether to apply monotone postprocessing.
#' @param type Either \code{"lambda"} or \code{"evr"}.
#' @param show_nonprivate Logical; whether to overlay nonprivate scree.
#' @param show_points Logical; whether to add points.
#' @param col_map Named color vector. Defaults are used when \code{NULL}.
#' @param lty_map Named line-type vector. Defaults are used when \code{NULL}.
#' @param pch_map Named point-shape vector. Defaults are used when \code{NULL}.
#' @param main,xlab,ylab,ylim Plot controls.
#' @param legend_pos Legend position.
#' @param ... Additional graphical arguments passed to \code{plot()}.
#'
#' @return Invisibly returns a list containing nonprivate results and
#'   method-specific \code{dp_scree()} outputs.
#' @export
plot_dp_scree <- function(
    X,
    k,
    dp_lambda = c("both", "clipped", "pmwm", "huber"),
    eps_total,
    delta_total,
    center = TRUE,
    standardize = TRUE,
    B_mult = 2,
    beta = 1.01,
    a = NULL,
    b = NULL,
    trim_const = 10,
    eta = 0.01,
    split_mode = TRUE,
    mu0 = 0,
    eta0 = 1,
    T = NULL,
    M = NULL,
    k_min_star = -40L,
    k_max_star = 40L,
    include_zero_bin = TRUE,
    m2_frac = 1 / 4,
    step_schedule = c("fixed", "1/sqrt(t)"),
    dp_v_flag = FALSE,
    cpp.option = FALSE,
    mono = TRUE,
    type = c("evr", "lambda"),
    show_nonprivate = TRUE,
    show_points = TRUE,
    col_map = NULL,
    lty_map = NULL,
    pch_map = NULL,
    main = NULL,
    xlab = "Component",
    ylab = NULL,
    ylim = NULL,
    legend_pos = "topright",
    ...
) {
  dp_lambda <- match.arg(dp_lambda)
  type <- match.arg(type)
  step_schedule <- match.arg(step_schedule)

  X <- as.matrix(X)
  validate_dp_pca_inputs(X, k, eps_total, delta_total)

  need_pmwm <- dp_lambda %in% c("both", "pmwm")
  if (need_pmwm && (is.null(a) || is.null(b))) {
    stop("For dp_lambda = 'pmwm' or 'both', you must provide 'a' and 'b'.")
  }

  if (is.null(col_map)) {
    col_map <- c(
      nonprivate = "black",
      clipped = "red",
      pmwm = "green",
      huber = "blue"
    )
  }
  if (is.null(lty_map)) {
    lty_map <- c(
      nonprivate = 1,
      clipped = 1,
      pmwm = 1,
      huber = 1
    )
  }
  if (is.null(pch_map)) {
    pch_map <- c(
      nonprivate = 16,
      clipped = 17,
      pmwm = 15,
      huber = 18
    )
  }

  results <- list()

  if (dp_lambda %in% c("both", "clipped")) {
    results$clipped <- dp_scree(
      X = X,
      k = k,
      method = "clipped",
      eps_total = eps_total,
      delta_total = delta_total,
      center = center,
      standardize = standardize,
      dp_v_flag = dp_v_flag,
      cpp.option = cpp.option,
      mono = mono,
      B_mult = B_mult
    )
  }

  if (dp_lambda %in% c("both", "pmwm")) {
    results$pmwm <- dp_scree(
      X = X,
      k = k,
      method = "pmwm",
      eps_total = eps_total,
      delta_total = delta_total,
      center = center,
      standardize = standardize,
      dp_v_flag = dp_v_flag,
      cpp.option = cpp.option,
      mono = mono,
      beta = beta,
      a = a,
      b = b,
      trim_const = trim_const,
      eta = eta,
      split_mode = split_mode
    )
  }

  if (dp_lambda %in% c("both", "huber")) {
    results$huber <- dp_scree(
      X = X,
      k = k,
      method = "huber",
      eps_total = eps_total,
      delta_total = delta_total,
      center = center,
      standardize = standardize,
      dp_v_flag = dp_v_flag,
      cpp.option = cpp.option,
      mono = mono,
      mu0 = mu0,
      eta0 = eta0,
      T = T,
      M = M,
      k_min_star = k_min_star,
      k_max_star = k_max_star,
      include_zero_bin = include_zero_bin,
      m2_frac = m2_frac,
      step_schedule = step_schedule
    )
  }

  # nonprivate reference: any computed result already contains it
  ref_obj <- if (length(results) > 0) results[[1]] else NULL

  y_list <- list()

  if (show_nonprivate) {
    if (is.null(ref_obj)) {
      stop("Nothing to plot: no DP result was computed.")
    }
    y_list$nonprivate <- if (type == "lambda") ref_obj$lambda_np else ref_obj$evr_np
  }

  for (nm in names(results)) {
    y_list[[nm]] <- if (type == "lambda") results[[nm]]$lambda_dp else results[[nm]]$evr
  }

  if (length(y_list) == 0) {
    stop("Nothing to plot. Check 'dp_lambda' and 'show_nonprivate'.")
  }

  if (is.null(ylab)) {
    ylab <- if (type == "lambda") "Eigenvalue" else "Explained Variance Ratio"
  }

  if (is.null(main)) {
    main <- if (type == "lambda") "DP Scree Plot" else "DP Explained Variance Ratio Plot"
  }

  if (is.null(ylim)) {
    y_all <- unlist(y_list, use.names = FALSE)
    y_all <- y_all[is.finite(y_all)]
    if (length(y_all) == 0) {
      ylim <- c(0, 1)
    } else {
      ylim <- range(y_all)
      if (diff(ylim) == 0) {
        ylim <- ylim + c(-0.5, 0.5)
      }
    }
  }

  idx <- seq_len(k)
  series_names <- names(y_list)

  nm1 <- series_names[1]
  y1 <- y_list[[nm1]]

  graphics::plot(
    idx, y1,
    type = if (show_points) "b" else "l",
    col = unname(col_map[nm1]),
    lty = unname(lty_map[nm1]),
    pch = unname(pch_map[nm1]),
    xlab = xlab,
    ylab = ylab,
    main = main,
    ylim = ylim,
    ...
  )

  if (length(series_names) >= 2) {
    for (nm in series_names[-1]) {
      graphics::lines(
        idx, y_list[[nm]],
        type = if (show_points) "b" else "l",
        col = unname(col_map[nm]),
        lty = unname(lty_map[nm]),
        pch = unname(pch_map[nm])
      )
    }
  }

  graphics::legend(
    legend_pos,
    legend = series_names,
    col = unname(col_map[series_names]),
    lty = unname(lty_map[series_names]),
    pch = unname(pch_map[series_names]),
    bty = "n"
  )

  invisible(list(
    nonprivate = if (is.null(ref_obj)) NULL else list(
      lambda = ref_obj$lambda_np,
      evr = ref_obj$evr_np
    ),
    results = results
  ))
}
