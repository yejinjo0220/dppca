# ============================================================
# dp_scree.R
# User-facing functions for differentially private scree estimation
# ============================================================

#' Control parameters for the clipped DP scree estimator
#'
#' @description
#' Creates a control list for \code{method = "clipped"} in \code{dp_scree()}
#' and \code{dp_scree_plot()}.
#'
#' @param C_clip Positive clipping threshold used by the clipped estimator.
#'
#' @return A list of clipped-estimator tuning parameters.
#'
#' @export
clipped_control <- function(C_clip = 3) {
  list(C_clip = C_clip)
}

#' Control parameters for the PMWM DP scree estimator
#'
#' @description
#' Creates a control list for \code{method = "pmwm"} in \code{dp_scree()}
#' and \code{dp_scree_plot()}.
#'
#' @param beta Log-binning base used by the PMWM private quantile estimator.
#' @param a,b Finite support bounds supplied to the PMWM private quantile routine.
#' @param trim_const,eta Practical clipping parameters used by the PMWM estimator.
#' @param split_mode Logical; whether the PMWM estimator splits the sample into
#'   quantile-estimation and mean-estimation subsets.
#'
#' @return A list of PMWM-estimator tuning parameters.
#'
#' @export
pmwm_control <- function(
    beta = 1.01, a = NULL, b = NULL,
    trim_const = 10, eta = 0.01, split_mode = TRUE
) {
  list(
    beta = beta,
    a = a,
    b = b,
    trim_const = trim_const,
    eta = eta,
    split_mode = split_mode
  )
}

#' Control parameters for the Huber DP scree estimator
#'
#' @description
#' Creates a control list for \code{method = "huber"} in \code{dp_scree()}
#' and \code{dp_scree_plot()}.
#'
#' @param mu0 Initial value used by the Huber noisy gradient descent.
#' @param eta0 Fixed step size used by the Huber noisy gradient descent.
#' @param T Optional integer number of Huber gradient descent iterations.
#' @param M Optional integer number of blocks used in \code{dp_m2()}.
#' @param k_min_m2 Integer lower bound for dyadic histogram bins used in the
#'   Huber scale-proxy step.
#' @param k_max_m2 Integer upper bound for dyadic histogram bins used in the
#'   Huber scale-proxy step.
#' @param m2_frac Fraction of the Huber scree privacy budget allocated to the
#'   private scale-proxy step.
#'
#' @return A list of Huber-estimator tuning parameters.
#'
#' @export
huber_control <- function(
    mu0 = 0, eta0 = 1, T = NULL, M = NULL,
    k_min_m2 = -40, k_max_m2 = 40,
    m2_frac = 1/4
) {
  list(
    mu0 = mu0,
    eta0 = eta0,
    T = T,
    M = M,
    k_min_m2 = k_min_m2,
    k_max_m2 = k_max_m2,
    m2_frac = m2_frac
  )
}

.default_scree_control <- function(method) {
  switch(
    method,
    clipped = clipped_control(),
    pmwm = pmwm_control(),
    huber = huber_control()
  )
}

.merge_scree_control <- function(method, control) {
  default <- .default_scree_control(method)

  if (is.null(control)) {
    control <- default
  } else {
    if (!is.list(control)) {
      stop("'control' must be a control list created by clipped_control(), pmwm_control(), or huber_control().")
    }
    control <- utils::modifyList(default, control, keep.null = TRUE)
  }

  if (method == "clipped" && is.null(control$C_clip)) {
    stop("For method = 'clipped', provide control = clipped_control(C_clip = ...).")
  }

  if (method == "pmwm" && (is.null(control$a) || is.null(control$b))) {
    stop("For method = 'pmwm', provide finite bounds using control = pmwm_control(a = ..., b = ...).")
  }

  control
}

.extract_plot_control <- function(control, method) {
  if (is.null(control)) {
    return(NULL)
  }

  if (!is.list(control)) {
    stop("'control' must be a control list, or a named list of control lists for dp_scree_plot().")
  }

  if (!is.null(control[[method]]) && is.list(control[[method]])) {
    return(control[[method]])
  }

  control
}

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
#'   \item \code{"pmwm"}: PMWM-style scree estimator based on private quantiles and winsorized means,
#'   \item \code{"huber"}: Huber-type robust private mean-based scree estimator.
#' }
#'
#' Method-specific tuning parameters are supplied through \code{control}:
#' \itemize{
#'   \item \code{control = clipped_control(C_clip = ...)} for \code{method = "clipped"},
#'   \item \code{control = pmwm_control(beta = ..., a = ..., b = ..., trim_const = ..., eta = ..., split_mode = ...)} for \code{method = "pmwm"},
#'   \item \code{control = huber_control(mu0 = ..., eta0 = ..., T = ..., M = ..., k_min_m2 = ..., k_max_m2 = ..., m2_frac = ...)} for \code{method = "huber"}.
#' }
#'
#' @param X Numeric data matrix with observations in rows.
#' @param k Integer number of leading principal components.
#' @param method One of \code{"clipped"}, \code{"pmwm"}, or \code{"huber"}.
#' @param eps_total Total privacy epsilon allocated to the full scree routine.
#' @param delta_total Total privacy delta allocated to the full scree routine.
#' @param center Logical; whether to center columns before PCA.
#' @param standardize Logical; whether to standardize columns before PCA.
#' @param g_dppca Logical; whether to privatize the PCA direction matrix.
#' @param cpp.option Logical passed to \code{mech_tau_sph()} when private directions are computed.
#' @param mono Logical; whether to apply monotone post-processing to the final DP scree vector.
#' @param control Method-specific control list created by \code{clipped_control()},
#'   \code{pmwm_control()}, or \code{huber_control()}.
#'
#' @return A list containing \code{method}, \code{scree_np}, \code{evr_np},
#'   \code{scree}, and \code{evr}.
#'
#' @examples
#' \dontrun{
#' dp_scree(X, k = 3, method = "clipped", eps_total = 1, delta_total = 1e-6,
#'          control = clipped_control(C_clip = 3))
#'
#' dp_scree(X, k = 3, method = "pmwm", eps_total = 1, delta_total = 1e-6,
#'          control = pmwm_control(a = 0, b = 10))
#'
#' dp_scree(X, k = 3, method = "huber", eps_total = 1, delta_total = 1e-6,
#'          control = huber_control(T = 50, M = 20))
#' }
#'
#' @export
dp_scree <- function(
    X, k, method = c("clipped", "pmwm", "huber"),
    eps_total, delta_total,
    center = TRUE, standardize = FALSE,
    g_dppca = FALSE, cpp.option = FALSE, mono = TRUE,
    control = NULL
) {
  method <- match.arg(method)
  control <- .merge_scree_control(method, control)

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
      C_clip = control$C_clip,
      g_dppca = g_dppca, cpp.option = cpp.option,
      mono = mono
    ),
    pmwm = dp_scree_pmwm(
      X = X, k = k,
      eps_total = eps_total, delta_total = delta_total,
      g_dppca = g_dppca, cpp.option = cpp.option,
      split_mode = control$split_mode,
      center = center, standardize = standardize,
      beta = control$beta, a = control$a, b = control$b,
      trim_const = control$trim_const, eta = control$eta,
      mono = mono
    ),
    huber = dp_scree_huber(
      X = X, k = k,
      eps_total = eps_total, delta_total = delta_total,
      g_dppca = g_dppca, cpp.option = cpp.option,
      center = center, standardize = standardize,
      mu0 = control$mu0, eta0 = control$eta0,
      T = control$T, M = control$M,
      k_min_m2 = control$k_min_m2,
      k_max_m2 = control$k_max_m2,
      m2_frac = control$m2_frac,
      mono = mono
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
#' By default, proportions of variance explained (PVE) are plotted. When
#' \code{type = "scree"}, raw scree values are plotted instead.
#'
#' Method-specific tuning parameters are supplied through \code{control}.
#' For a single method, use a single control object, for example
#' \code{control = pmwm_control(a = ..., b = ...)}. For
#' \code{dp_scree_method = "all"}, use a named list, for example
#' \code{control = list(clipped = clipped_control(...), pmwm = pmwm_control(...), huber = huber_control(...))}.
#'
#' @param X Numeric data matrix with observations in rows.
#' @param k Integer number of leading principal components.
#' @param dp_scree_method Which DP estimator(s) to plot. One of
#'   \code{"all"}, \code{"clipped"}, \code{"pmwm"}, or \code{"huber"}.
#' @param eps_total Total privacy epsilon allocated to the full scree routine.
#' @param delta_total Total privacy delta allocated to the full scree routine.
#' @param center Logical; whether to center columns before PCA.
#' @param standardize Logical; whether to standardize columns before PCA.
#' @param control Method-specific control list, or a named list of control lists
#'   when \code{dp_scree_method = "all"}.
#' @param g_dppca Logical; whether to privatize the PCA direction matrix.
#' @param cpp.option Logical passed to \code{mech_tau_sph()} when private directions are computed.
#' @param mono Logical; whether to apply monotone post-processing.
#' @param type Either \code{"pve"} or \code{"scree"}.
#'
#' @details
#' Plot appearance is handled internally. The non-private curve is overlaid,
#' points are shown on each curve, and colors, line types, point shapes,
#' labels, limits, and legend position are set automatically.
#'
#' @return Invisibly returns a list containing non-private results and the
#'   method-specific \code{dp_scree()} outputs that were plotted.
#'
#' @examples
#' \dontrun{
#' dp_scree_plot(
#'   X, k = 3, dp_scree_method = "all",
#'   eps_total = 1, delta_total = 1e-6,
#'   control = list(
#'     clipped = clipped_control(C_clip = 3),
#'     pmwm = pmwm_control(a = 0, b = 10),
#'     huber = huber_control(T = 50, M = 20)
#'   )
#' )
#' }
#'
#' @export
dp_scree_plot <- function(
    X, k, dp_scree_method = c("clipped", "pmwm", "huber", "all"),
    eps_total, delta_total,
    center = TRUE, standardize = FALSE,
    control = NULL,
    g_dppca = FALSE, cpp.option = FALSE, mono = TRUE,
    type = c("pve", "scree")
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
  need_pmwm <- dp_scree_method %in% c("all", "pmwm")
  need_huber <- dp_scree_method %in% c("all", "huber")

  col_map <- c(
    nonprivate = "black",
    clipped = "red",
    pmwm = "forestgreen",
    huber = "blue"
  )
  lty_map <- c(
    nonprivate = 1,
    clipped = 1,
    pmwm = 1,
    huber = 1
  )
  pch_map <- c(
    nonprivate = 16,
    clipped = 17,
    pmwm = 15,
    huber = 18
  )
  legend_pos <- "topright"

  results <- list()

  if (need_clipped) {
    results$clipped <- dp_scree(
      X = X, k = k, method = "clipped",
      eps_total = eps_total, delta_total = delta_total,
      center = center, standardize = standardize,
      g_dppca = g_dppca, cpp.option = cpp.option, mono = mono,
      control = .extract_plot_control(control, "clipped")
    )
  }

  if (need_pmwm) {
    results$pmwm <- dp_scree(
      X = X, k = k, method = "pmwm",
      eps_total = eps_total, delta_total = delta_total,
      center = center, standardize = standardize,
      g_dppca = g_dppca, cpp.option = cpp.option, mono = mono,
      control = .extract_plot_control(control, "pmwm")
    )
  }

  if (need_huber) {
    results$huber <- dp_scree(
      X = X, k = k, method = "huber",
      eps_total = eps_total, delta_total = delta_total,
      center = center, standardize = standardize,
      g_dppca = g_dppca, cpp.option = cpp.option, mono = mono,
      control = .extract_plot_control(control, "huber")
    )
  }

  ref_obj <- if (length(results) > 0) results[[1]] else NULL
  if (is.null(ref_obj)) stop("Nothing to plot: no DP result was computed.")

  y_list <- list()
  y_list$nonprivate <- if (type == "scree") ref_obj$scree_np else ref_obj$evr_np
  for (nm in names(results)) {
    y_list[[nm]] <- if (type == "scree") results[[nm]]$scree else results[[nm]]$evr
  }

  ylab <- if (type == "scree") "Scree Value" else "Proportion of Variance Explained"
  main <- if (type == "scree") "DP Scree Plot" else "DP PVE Plot"
  xlab <- "Component"

  y_all <- unlist(y_list, use.names = FALSE)
  y_all <- y_all[is.finite(y_all)]
  if (length(y_all) == 0) {
    ylim <- c(0, 1)
  } else {
    ylim <- range(y_all)
    if (diff(ylim) == 0) ylim <- ylim + c(-0.5, 0.5)
  }

  idx <- seq_len(k)
  series_names <- names(y_list)
  nm1 <- series_names[1]
  y1 <- y_list[[nm1]]

  graphics::plot(
    idx, y1,
    type = "b",
    col = unname(col_map[nm1]),
    lty = unname(lty_map[nm1]),
    pch = unname(pch_map[nm1]),
    xlab = xlab, ylab = ylab, main = main, ylim = ylim
  )

  if (length(series_names) >= 2) {
    for (nm in series_names[-1]) {
      graphics::lines(
        idx, y_list[[nm]],
        type = "b",
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
    nonprivate = list(scree = ref_obj$scree_np, pve = ref_obj$evr_np),
    results = results
  ))
}

