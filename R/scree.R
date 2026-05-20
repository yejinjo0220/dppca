# ============================================================-
# scree.R
# functions for differentially private scree estimation
# ============================================================-

#' Differentially private scree values
#'
#' This function computes estimates of scree values, eigenvalues of covariance
#' matrix, for principal component analysis,
#' including both the usual non-private estimates and differentially private
#' estimates.
#' The private estimates are computed as the private mean of the squared
#' principal component scores.
#' See Details for the estimating equations and method-specific construction.
#'
#' @param X A numeric matrix or data frame. Rows correspond to observations and
#'  columns correspond to variables.
#' @param k Positive integer defining the number of leading principal components
#'   to estimate. Must be an integer between `1` and the number of columns in `X`.
#' @param method Scree value estimation method. One of `"clipped"`, `"pmwm"`, or
#'  `"huber"`.
#' @param control Optional method-specific control list created by
#'   [clipped_control()], [pmwm_control()], or [huber_control()].
#' @param eps Positive number defining the total `epsilon` privacy parameter.
#'   If `g_dppca = TRUE`, it is split between private direction estimation and
#'   private scree estimation.
#' @param delta Number in `(0, 1)` defining the total `delta` privacy parameter.
#'   If `g_dppca = TRUE`, it is split between private direction estimation and
#'   private scree estimation.
#' @param center A logical value indicating whether to center the columns of `X`
#'   before computing principal component directions. The default is `TRUE`.
#' @param standardize A logical value indicating whether to scale the columns of
#'   `X` by their sample standard deviations after optional centering. The
#'   default is `FALSE`.
#' @param g_dppca A logical value indicating whether to use private principal
#'   component directions for scree estimation. The default is `FALSE`. See
#'   [dp_pc_dir()] for details.
#' @param cpp.option A logical value passed to [dp_pc_dir()] when
#'   `g_dppca = TRUE`. The default is `FALSE`.
#' @param mono A logical value indicating whether to apply monotone
#'   post-processing to the vector of private scree values. The default is `TRUE`.
#'
#' @details
#' Let \eqn{X} denote the preprocessed data matrix and let \eqn{v_l} be the \eqn{l}th
#' principal component direction. The \eqn{l}th score vector is
#' \eqn{z_l = X v_l}. The corresponding sample scree value can be written as
#' \deqn{
#'   \hat{\lambda}_l
#'   = v_l^\top \widehat{\Sigma} v_l
#'   = \frac{1}{n - 1}\sum_{i = 1}^n z_{il}^2
#'   = \frac{n}{n - 1}\left(\frac{1}{n}\sum_{i = 1}^n w_{il}\right),
#'   \qquad w_{il} = z_{il}^2.
#' }
#' Therefore, each scree value is estimated by privately estimating the mean of
#' \eqn{w_{1l}, \ldots, w_{nl}} and multiplying by \eqn{n/(n - 1)}.
#'
#' The supported methods differ in how this private mean is estimated:
#' \itemize{
#'   \item `"clipped"` clips the squared scores \eqn{w_{i\ell}} at `C_clip` and
#'   then applies the Gaussian mechanism \insertCite{dwork2014algorithmic}{dppca}.
#'   This is the simplest option but depends directly on the clipping threshold.
#'   \item `"pmwm"` uses the private modified winsorized mean approach of
#'   \insertCite{ramsay2025pmw;textual}{dppca}, adapted from the accompanying
#'   Python implementation into R. It privately estimates tail cutoffs,
#'   winsorizes the squared scores \eqn{w_{i\ell}}, and releases a noisy
#'   winsorized mean.
#'   \item `"huber"` uses a Huber-type private robust mean estimator based on
#'   noisy gradient descent, following \insertCite{yu2024gaussian;textual}{dppca}.
#' }
#'
#' The argument `g_dppca` controls how the principal component directions are
#' obtained. If `g_dppca = FALSE`, the directions are computed non-privately as
#' an eigenvector of sample covariance, and the full privacy parameters
#' `eps` and `delta` are used for private scree value estimation. If `g_dppca = TRUE`,
#' the directions are computed privately using [dp_pc_dir()].
#' In that case, the privacy parameters are split equally:
#' `eps / 2` and `delta / 2` are used for private direction estimation, and the
#' remaining `eps / 2` and `delta / 2` are used for private scree value estimation.
#' When `mono = TRUE`, the final monotone adjustment is a post-processing step and
#' does not change the privacy guarantee.
#'
#' For a detailed procedure and mathematical formulations,
#' refer \url{https://yejinjo0220.github.io/dppca/articles/dp_scree}.
#'
#' @return A list with components:
#' \itemize{
#'   \item `method`: scree value estimation method.
#'   \item `scree_np`: non-private scree estimates.
#'   \item `pve_np`: non-private proportions of variance explained.
#'   \item `scree`: differentially private scree value estimates.
#'   \item `pve`: differentially private proportions of variance explained.
#' }
#'
#' @seealso
#' [dp_pc_dir()] for principal component direction estimation.
#' [clipped_control()], [pmwm_control()], and [huber_control()] for
#' method-specific tuning parameters.
#'
#' @references
#' \insertRef{dwork2014algorithmic}{dppca}
#'
#' \insertRef{ramsay2025pmw}{dppca}
#'
#' \insertRef{yu2024gaussian}{dppca}
#'
#' \insertRef{kim2025robustdppca}{dppca}
#'
#' @examples
#' data(gau, package = "dppca")
#'
#' # Use a small subset to keep the example fast.
#' X <- gau[1:100, ]
#'
#' # Estimate the private scree values using the clipped mean method.
#' set.seed(123)
#' dp_scree(
#'   X,
#'   k = 2,
#'   method = "clipped",
#'   control = clipped_control(C_clip = 3),
#'   eps = 2,
#'   delta = 1e-3
#' )
#'
#' # Other scree methods can be used by changing `method` and `control`, e.g.,
#' # method = "pmwm",
#' # control = pmwm_control(a = 0, b = 50, trim_const = 10, eta = 0.01)
#' #
#' # method = "huber",
#' # control = huber_control(k_min_m2 = -10, k_max_m2 = 10, m2_frac = 1 / 4)
#'
#' @export
dp_scree <- function(
    X, k, method = c("clipped", "pmwm", "huber"),
    control = NULL,
    eps, delta,
    center = TRUE, standardize = FALSE,
    g_dppca = FALSE, cpp.option = FALSE, mono = TRUE
) {
  method <- match.arg(
    method,
    choices = c("clipped", "pmwm", "huber"),
    several.ok = TRUE
  )
  method <- unique(method)
  control <- .merge_scree_control(method, control)

  X <- as.matrix(X)
  validate_scree_inputs(
    X = X,
    k = k,
    eps = eps,
    delta = delta
  )

  if (!requireNamespace("rARPACK", quietly = TRUE)) {
    stop("Package `rARPACK` is required.", call. = FALSE)
  }

  X_proc <- prep_matrix_for_pca(
    X = X,
    center = center,
    standardize = standardize
  )

  S_np <- stats::cov(X_proc)
  eig_np <- rARPACK::eigs_sym(S_np, k = k)

  scree_np <- as.numeric(eig_np$values)
  pve_np <- scree_to_pve(scree_np)

  result <- switch(
    method,
    clipped = dp_scree_clipped(
      X = X, k = k,
      eps = eps, delta = delta,
      center = center, standardize = standardize,
      C_clip = control$C_clip,
      g_dppca = g_dppca, cpp.option = cpp.option,
      mono = mono
    ),
    pmwm = dp_scree_pmwm(
      X = X, k = k,
      eps = eps, delta = delta,
      g_dppca = g_dppca, cpp.option = cpp.option,
      split_mode = control$split_mode,
      center = center, standardize = standardize,
      beta = control$beta, a = control$a, b = control$b,
      trim_const = control$trim_const, eta = control$eta,
      mono = mono
    ),
    huber = dp_scree_huber(
      X = X, k = k,
      eps = eps, delta = delta,
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
    pve_np   = pve_np,
    scree    = result$scree,
    pve      = result$pve
  )
}

#' Plot differentially private scree estimates
#'
#' @description
#' This function computes and visualizes scree curves for principal component
#' analysis, including the usual non-private curve and one or more
#' differentially private estimates. It is a plotting wrapper around
#' [dp_scree()] and returns a `ggplot` object.
#'
#' @param X A numeric matrix or data frame. Rows correspond to observations and
#'  columns correspond to variables.
#' @param k Positive integer defining the number of leading principal components
#'   to estimate. Must be an integer between `1` and the number of columns in `X`.
#' @param method Scree estimation method or methods to plot. One or more of
#'   `"clipped"`, `"pmwm"`, or `"huber"`.
#' @param control Optional method-specific control list, or a named list of
#'   control lists when multiple methods are requested. Use [clipped_control()],
#'   [pmwm_control()], and [huber_control()].
#' @param eps Positive number defining the total `epsilon` privacy parameter.
#'   If `g_dppca = TRUE`, it is split between private direction estimation and
#'   private scree estimation.
#' @param delta Number in `(0, 1)` defining the total `delta` privacy parameter.
#'   If `g_dppca = TRUE`, it is split between private direction estimation and
#'   private scree estimation.
#' @param center A logical value indicating whether to center the columns of `X`
#'   before computing principal component directions. The default is `TRUE`.
#' @param standardize A logical value indicating whether to scale the columns of
#'   `X` by their sample standard deviations after optional centering. The
#'   default is `FALSE`.
#' @param g_dppca A logical value indicating whether to use private principal
#'   component directions for scree estimation. The default is `FALSE`. See
#'   [dp_pc_dir()] for details.
#' @param cpp.option A logical value passed to [dp_pc_dir()] when
#'   `g_dppca = TRUE`. The default is `FALSE`.
#' @param mono A logical value indicating whether to apply monotone
#'   post-processing to the private scree vector. The default is `TRUE`.
#' @param type Quantity to plot. Use `"pve"` to plot proportions of variance
#'   explained and `"scree"` to plot raw scree values. The default is `"pve"`.
#'
#' @details
#' This function is a plotting wrapper around [dp_scree()]. For each requested
#' method, it computes a private scree estimate and overlays it with the
#' corresponding non-private curve. When `type = "pve"`, the plotted quantity is
#' the proportion of variance explained (PVE); when `type = "scree"`, the raw
#' scree values are shown.
#'
#' To plot multiple methods, pass a character vector to `method`. If a method
#' requires tuning parameters, pass `control` as a named list, for example
#' `control = list(clipped = clipped_control(), pmwm = pmwm_control(),
#' huber = huber_control())`.
#'
#' For the estimating equations, privacy-budget allocation, and method-specific
#' construction, see [dp_scree()].
#'
#' @return Invisibly returns a list with components:
#' \itemize{
#'   \item `nonprivate`: non-private scree and PVE values.
#'   \item `results`: method-specific [dp_scree()] outputs used in the plot.
#' }
#'
#' @seealso
#' [dp_pc_dir()] for principal component direction estimation.
#' [dp_scree()] for computing non-private and differentially private scree
#' estimates.
#' [clipped_control()], [pmwm_control()], and [huber_control()] for
#' method-specific tuning parameters.
#'
#' @references
#' \insertRef{dwork2014algorithmic}{dppca}
#'
#' \insertRef{ramsay2025pmw}{dppca}
#'
#' \insertRef{yu2024gaussian}{dppca}
#'
#' \insertRef{kim2025robustdppca}{dppca}
#'
#' @examples
#' data(gau, package = "dppca")
#'
#' # Use a small subset to keep the example fast.
#' X <- gau[1:200, ]
#'
#' # Draw a private scree plot using the clipped mean method.
#' dp_scree_plot(
#'   X,
#'   k = 3,
#'   method = "clipped",
#'   control = clipped_control(C_clip = 3),
#'   eps = 3,
#'   delta = 1e-3
#' )

#' # Multiple scree methods can be overlaid by passing a vector to `method`
#' # and a named list to `control`, for example:
#' #
#' # dp_scree_plot(
#' #   X,
#' #   k = 2,
#' #   method = c("clipped", "pmwm", "huber"),
#' #   control = list(
#' #     clipped = clipped_control(C_clip = 3),
#' #     pmwm = pmwm_control(a = 0, b = 50, trim_const = 10, eta = 0.01),
#' #     huber = huber_control(k_min_m2 = -10, k_max_m2 = 10, m2_frac = 1 / 4)
#' #   ),
#' #   eps = 3,
#' #   delta = 1e-3
#' # )
#' @export
dp_scree_plot <- function(
    X,
    k,
    method = c("clipped", "pmwm", "huber"),
    control = NULL,
    eps, delta,
    center = TRUE, standardize = FALSE,
    g_dppca = FALSE, cpp.option = FALSE, mono = TRUE,
    type = c("pve", "scree")
) {
  method <- match.arg(method)
  type <- match.arg(type)

  X <- as.matrix(X)
  validate_scree_inputs(
    X = X,
    k = k,
    eps = eps,
    delta = delta
  )

  need_clipped <- "clipped" %in% method
  need_pmwm <- "pmwm" %in% method
  need_huber <- "huber" %in% method

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
      eps = eps, delta = delta,
      center = center, standardize = standardize,
      g_dppca = g_dppca, cpp.option = cpp.option, mono = mono,
      control = .extract_plot_control(control, "clipped")
    )
  }

  if (need_pmwm) {
    results$pmwm <- dp_scree(
      X = X, k = k, method = "pmwm",
      eps = eps, delta = delta,
      center = center, standardize = standardize,
      g_dppca = g_dppca, cpp.option = cpp.option, mono = mono,
      control = .extract_plot_control(control, "pmwm")
    )
  }

  if (need_huber) {
    results$huber <- dp_scree(
      X = X, k = k, method = "huber",
      eps = eps, delta = delta,
      center = center, standardize = standardize,
      g_dppca = g_dppca, cpp.option = cpp.option, mono = mono,
      control = .extract_plot_control(control, "huber")
    )
  }

  ref_obj <- if (length(results) > 0) results[[1]] else NULL
  if (is.null(ref_obj)) stop("Nothing to plot: no DP result was computed.")

  y_list <- list()
  y_list$nonprivate <- if (type == "scree") ref_obj$scree_np else ref_obj$pve_np
  for (nm in names(results)) {
    y_list[[nm]] <- if (type == "scree") results[[nm]]$scree else results[[nm]]$pve
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
    nonprivate = list(scree = ref_obj$scree_np, pve = ref_obj$pve_np),
    results = results
  ))
}
