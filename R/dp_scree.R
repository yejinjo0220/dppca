# ============================================================
# dp_scree_final.R
# User-facing functions for DP scree estimation
# ============================================================

#' Compute differentially private scree estimates
#'
#' @description
#' User-facing function for differentially private scree estimation.
#' It computes the leading non-private scree values together with one selected
#' DP scree estimator.
#'
#' Supported methods are:
#' \itemize{
#'   \item \code{"clipped"}: clipped mean-based scree estimator,
#'   \item \code{"pmw"}: PMW-style scree estimator based on private quantiles and winsorized means,
#'   \item \code{"huber"}: Huber-type robust private mean-based scree estimator.
#' }
#'
#' @param X Numeric data matrix with observations in rows.
#' @param k Integer number of leading principal components.
#' @param method One of \code{"clipped"}, \code{"pmw"}, or \code{"huber"}.
#' @param eps_total Total privacy epsilon allocated to the full scree routine.
#' @param delta_total Total privacy delta allocated to the full scree routine.
#' @param center Logical; whether to center columns before PCA.
#' @param standardize Logical; whether to standardize columns before PCA.
#' @param g_dppca Logical; whether to privatize the PCA direction matrix.
#' @param cpp.option Logical passed to \code{mech_tau_sph()} when private directions are computed.
#' @param mono Logical; whether to apply monotone post-processing to the final DP scree vector.
#' @param C_clip Positive clipping threshold used by the clipped estimator.
#' @param beta Log-binning base used by the PMW quantile estimator.
#' @param a,b Finite support bounds supplied to the PMW quantile routine.
#' @param trim_const,eta PMW practical clipping parameters.
#' @param split_mode Logical; whether the PMW estimator splits the sample.
#' @param mu0 Initial value used by the Huber noisy gradient descent.
#' @param eta0 Fixed step size used by the Huber noisy gradient descent.
#' @param T Optional integer number of Huber gradient descent iterations.
#' @param M Optional integer number of blocks used in \code{dp_m2()}.
#' @param k_min_m2 Integer lower bound for dyadic histogram bins used in the Huber scale-proxy step.
#' @param k_max_m2 Integer upper bound for dyadic histogram bins used in the Huber scale-proxy step.
#' @param m2_frac Fraction of the Huber scree privacy budget allocated to the private scale-proxy step.
#'
#' @return A list containing \code{method}, \code{scree_np}, \code{evr_np},
#'   \code{scree}, and \code{evr}.
#'
#' @export
dp_scree <- function(
    X, k, method = c("clipped", "pmw", "huber"),
    eps_total, delta_total,
    center = TRUE, standardize = FALSE,
    g_dppca = FALSE, cpp.option = FALSE, mono = TRUE,
    C_clip = 3,
    beta = 1.01, a = NULL, b = NULL,
    trim_const = 10, eta = 0.01, split_mode = TRUE,
    mu0 = 0, eta0 = 1, T = NULL, M = NULL,
    k_min_m2 = -40, k_max_m2 = 40,
    m2_frac = 1/4
) {
  method <- match.arg(method)

  X <- as.matrix(X)
  validate_scree_inputs(
    X = X,
    k = k,
    eps_total = eps_total,
    delta_total = delta_total
  )

  if (!requireNamespace("rARPACK", quietly = TRUE)) {
    stop("Package 'rARPACK' is required.")
  }

  if (method == "clipped" && is.null(C_clip)) {
    stop("For method = 'clipped', you must provide 'C_clip'.")
  }

  if (method == "pmw" && (is.null(a) || is.null(b))) {
    stop("For method = 'pmw', you must provide 'a' and 'b'.")
  }

  X_proc <- prep_matrix_for_pca(
    X = X,
    center = center,
    standardize = standardize
  )

  S_np <- stats::cov(X_proc)
  eig_np <- rARPACK::eigs_sym(S_np, k = k)

  scree_np <- as.numeric(eig_np$values)
  evr_np <- scree_to_evr(scree_np)

  result <- switch(
    method,
    clipped = dp_scree_clipped(
      X = X, k = k,
      eps_total = eps_total, delta_total = delta_total,
      center = center, standardize = standardize,
      C_clip = C_clip,
      g_dppca = g_dppca, cpp.option = cpp.option,
      mono = mono
    ),
    pmw = dp_scree_pmw(
      X = X, k = k,
      eps_total = eps_total, delta_total = delta_total,
      g_dppca = g_dppca, cpp.option = cpp.option,
      split_mode = split_mode,
      center = center, standardize = standardize,
      beta = beta, a = a, b = b,
      trim_const = trim_const, eta = eta,
      mono = mono
    ),
    huber = dp_scree_huber(
      X = X, k = k,
      eps_total = eps_total, delta_total = delta_total,
      g_dppca = g_dppca, cpp.option = cpp.option,
      center = center, standardize = standardize,
      mu0 = mu0, eta0 = eta0, T = T, M = M,
      k_min_m2 = k_min_m2, k_max_m2 = k_max_m2,
      m2_frac = m2_frac, mono = mono
    )
  )

  list(
    method   = method,
    scree_np = scree_np,
    evr_np   = evr_np,
    scree    = result$scree,
    evr      = result$evr
  )
}

#' Plot differentially private scree curves
#'
#' @description
#' User-facing function that computes and plots non-private and/or
#' differentially private scree curves using one or more DP scree estimators.
#'
#' By default, explained variance ratios (EVR) are plotted. When
#' \code{type = "scree"}, raw scree values are plotted instead.
#'
#' @param X Numeric data matrix with observations in rows.
#' @param k Integer number of leading principal components.
#' @param dp_scree_method Which DP estimator(s) to plot. One of
#'   \code{"all"}, \code{"clipped"}, \code{"pmw"}, or \code{"huber"}.
#' @param eps_total Total privacy epsilon allocated to the full scree routine.
#' @param delta_total Total privacy delta allocated to the full scree routine.
#' @param center Logical; whether to center columns before PCA.
#' @param standardize Logical; whether to standardize columns before PCA.
#' @param C_clip Positive clipping threshold used by the clipped estimator.
#' @param beta Log-binning base used by the PMW quantile estimator.
#' @param a,b Finite support bounds supplied to the PMW quantile routine.
#' @param trim_const,eta PMW practical clipping parameters.
#' @param split_mode Logical; whether the PMW estimator splits the sample into quantile and mean subsets.
#' @param mu0 Initial value used by the Huber noisy gradient descent.
#' @param eta0 Fixed step size used by the Huber noisy gradient descent.
#' @param T Optional integer number of Huber gradient descent iterations.
#' @param M Optional integer number of blocks used in \code{dp_m2()}.
#' @param k_min_m2 Integer lower bound for dyadic histogram bins used in the Huber scale-proxy step.
#' @param k_max_m2 Integer upper bound for dyadic histogram bins used in the Huber scale-proxy step.
#' @param m2_frac Fraction of the Huber scree privacy budget allocated to the private scale-proxy step.
#' @param g_dppca Logical; whether to privatize the PCA direction matrix.
#' @param cpp.option Logical passed to \code{mech_tau_sph()} when private directions are computed.
#' @param mono Logical; whether to apply monotone post-processing.
#' @param type Either \code{"evr"} or \code{"scree"}.
#' @param show_nonprivate Logical; whether to overlay the non-private curve.
#' @param show_points Logical; whether to draw points on the plotted curves.
#' @param col_map Optional named color vector.
#' @param lty_map Optional named line-type vector.
#' @param pch_map Optional named point-shape vector.
#' @param main,xlab,ylab,ylim Plot controls.
#' @param legend_pos Legend position.
#' @param ... Additional graphical arguments passed to \code{plot()}.
#'
#' @return Invisibly returns a list containing non-private results and the
#'   method-specific \code{dp_scree()} outputs that were plotted.
#'
#' @export
dp_scree_plot <- function(
    X, k, dp_scree_method = c("all", "clipped", "pmw", "huber"),
    eps_total, delta_total,
    center = TRUE, standardize = TRUE,
    C_clip = 3,
    beta = 1.01, a = NULL, b = NULL,
    trim_const = 10, eta = 0.01, split_mode = TRUE,
    mu0 = 0, eta0 = 1, T = NULL, M = NULL,
    k_min_m2 = -40, k_max_m2 = 40,
    m2_frac = 1/4,
    g_dppca = FALSE, cpp.option = FALSE, mono = TRUE,
    type = c("evr", "scree"), show_nonprivate = TRUE, show_points = TRUE,
    col_map = NULL, lty_map = NULL, pch_map = NULL,
    main = NULL, xlab = "Component", ylab = NULL, ylim = NULL,
    legend_pos = "topright", ...
) {
  dp_scree_method <- match.arg(dp_scree_method)
  type <- match.arg(type)

  X <- as.matrix(X)
  validate_scree_inputs(
    X = X,
    k = k,
    eps_total = eps_total,
    delta_total = delta_total
  )

  need_clipped <- dp_scree_method %in% c("all", "clipped")
  need_pmw <- dp_scree_method %in% c("all", "pmw")
  need_huber <- dp_scree_method %in% c("all", "huber")

  if (need_clipped && is.null(C_clip)) {
    stop("For dp_scree_method = 'clipped' or 'all', you must provide 'C_clip'.")
  }
  if (need_pmw && (is.null(a) || is.null(b))) {
    stop("For dp_scree_method = 'pmw' or 'all', you must provide 'a' and 'b'.")
  }

  if (is.null(col_map)) col_map <- c(nonprivate = "black", clipped = "red", pmw = "forestgreen", huber = "blue")
  if (is.null(lty_map)) lty_map <- c(nonprivate = 1, clipped = 1, pmw = 1, huber = 1)
  if (is.null(pch_map)) pch_map <- c(nonprivate = 16, clipped = 17, pmw = 15, huber = 18)

  results <- list()

  if (need_clipped) {
    results$clipped <- dp_scree(
      X = X, k = k, method = "clipped",
      eps_total = eps_total, delta_total = delta_total,
      center = center, standardize = standardize,
      g_dppca = g_dppca, cpp.option = cpp.option, mono = mono,
      C_clip = C_clip
    )
  }

  if (need_pmw) {
    results$pmw <- dp_scree(
      X = X, k = k, method = "pmw",
      eps_total = eps_total, delta_total = delta_total,
      center = center, standardize = standardize,
      g_dppca = g_dppca, cpp.option = cpp.option, mono = mono,
      beta = beta, a = a, b = b,
      trim_const = trim_const, eta = eta, split_mode = split_mode
    )
  }

  if (need_huber) {
    results$huber <- dp_scree(
      X = X, k = k, method = "huber",
      eps_total = eps_total, delta_total = delta_total,
      center = center, standardize = standardize,
      g_dppca = g_dppca, cpp.option = cpp.option, mono = mono,
      mu0 = mu0, eta0 = eta0, T = T, M = M,
      k_min_m2 = k_min_m2, k_max_m2 = k_max_m2,
      m2_frac = m2_frac
    )
  }

  ref_obj <- if (length(results) > 0) results[[1]] else NULL
  if (is.null(ref_obj)) stop("Nothing to plot: no DP result was computed.")

  y_list <- list()
  if (show_nonprivate) y_list$nonprivate <- if (type == "scree") ref_obj$scree_np else ref_obj$evr_np
  for (nm in names(results)) y_list[[nm]] <- if (type == "scree") results[[nm]]$scree else results[[nm]]$evr

  if (is.null(ylab)) ylab <- if (type == "scree") "Scree Value" else "Explained Variance Ratio"
  if (is.null(main)) main <- if (type == "scree") "DP Scree Plot" else "DP Explained Variance Ratio Plot"

  if (is.null(ylim)) {
    y_all <- unlist(y_list, use.names = FALSE)
    y_all <- y_all[is.finite(y_all)]
    if (length(y_all) == 0) {
      ylim <- c(0, 1)
    } else {
      ylim <- range(y_all)
      if (diff(ylim) == 0) ylim <- ylim + c(-0.5, 0.5)
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
    xlab = xlab, ylab = ylab, main = main, ylim = ylim, ...
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
    nonprivate = list(scree = ref_obj$scree_np, evr = ref_obj$evr_np),
    results = results
  ))
}
