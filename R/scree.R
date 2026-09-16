# ============================================================-
# scree.R
# functions for differentially private scree estimation
# ============================================================-

#' Differentially private scree values
#'
#' This function computes estimates of scree values, eigenvalues of the
#' covariance matrix, for principal component analysis. It returns the usual
#' non-private estimates once, together with one or more differentially private
#' estimates.
#'
#' The private estimates are computed as private estimates of the mean of the
#' squared principal component scores. See Details for the estimating equations
#' and method-specific construction.
#'
#' @param X A numeric matrix or data frame. Rows correspond to observations and
#'   columns correspond to variables.
#' @param k Positive integer defining the number of leading principal components
#'   to estimate. Must be an integer between `1` and the number of columns in `X`.
#' @param privacy A [privacy_control()] object with exactly the `loading` and
#'   `scree` components. For example, use
#'   `privacy_control(eps = 3, delta = 1e-4, split = c(loading = 0.5, scree = 0.5))`.
#'   The loading budget is used once; each requested scree method receives
#'   the full `scree` component budget, without division across methods.
#' @param method Scree value estimation method or methods. One or more of
#'   `"clipped"`, `"pmwm"`, or `"huber"`. If omitted, `"clipped"` is used.
#' @param control Optional method-specific control list created by
#'   [clipped_control()], [pmwm_control()], or [huber_control()]. When multiple
#'   methods are requested, use a named list with method names.
#' @param center A logical value indicating whether to center the columns of `X`
#'   before computing principal component directions. The default is `TRUE`.
#' @param standardize A logical value indicating whether to scale the columns of
#'   `X` by their sample standard deviations after optional centering. The
#'   default is `FALSE`.
#' @param cpp.option A logical value passed to [dp_loading()] for private
#'   loading estimation. The default is `FALSE`.
#' @param mono A logical value indicating whether to apply monotone
#'   post-processing to the vector of private scree values. The default is
#'   `TRUE`.
#'
#' @details
#' Let \eqn{X} denote the preprocessed data matrix and let \eqn{v_l} be the
#' \eqn{l}th principal component direction. The \eqn{l}th score vector is
#' \eqn{z_l = X v_l}, with sample mean \eqn{\bar z_l}. The corresponding sample
#' scree value can be written as
#' \deqn{
#'   \hat{\lambda}_l
#'   = v_l^\top \widehat{\Sigma} v_l
#'   = \frac{1}{n - 1}\sum_{i = 1}^n (z_{il} - \bar z_l)^2
#'   = \frac{n}{n - 1}\left(\frac{1}{n}\sum_{i = 1}^n w_{il}\right),
#'   \qquad w_{il} = (z_{il} - \bar z_l)^2.
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
#'   Python implementation into R. The lower and upper tail cutoffs are
#'   estimated with the pure-DP fully unbounded quantile mechanism of
#'   \insertCite{durfee2023unbounded;textual}{dppca} and truncated to the public
#'   bounds supplied through [pmwm_control()].
#'   The squared scores \eqn{w_{i\ell}} are then winsorized to those cutoffs,
#'   and the final winsorized mean is released with a Gaussian mechanism
#'   calibrated to its share of the `scree` component budget.
#'   \item `"huber"` uses a Huber-type private robust mean estimator based on
#'   noisy gradient descent, following \insertCite{yu2024gaussian;textual}{dppca}.
#' }
#'
#' Principal component directions are always estimated privately using
#' [dp_loading()]. The data are preprocessed once, and one loading estimate
#' is computed using the `loading` component of `privacy`. The first `k`
#' columns of its `private` matrix produce one score matrix shared by all
#' requested methods. Each method uses the `scree` component directly; there
#' is no additional loading allocation within a scree method.
#'
#' When multiple methods are requested, each receives the same full `scree`
#' epsilon and delta budgets for method comparison. These budgets are not
#' divided across methods. For a joint release, composition counts the shared
#' loading allocation once, plus the privacy costs of all requested methods.
#' The total recorded in `privacy` covers loading plus one scree method, not
#' a joint release of several scree methods.
#'
#' For `method = "pmwm"`, the private quantile step is pure DP and does not
#' consume `delta`; the scree component's delta budget is used by the Gaussian
#' winsorized-mean release. For `method = "huber"`, the private scale-proxy step
#' is pure DP. Noisy gradient descent receives the `(1 - m2_frac)` share of
#' the scree component's delta budget; the remaining `m2_frac` share is unused.
#'
#' The `nonprivate` component is provided only as a non-private reference and is
#' not itself differentially private.
#'
#' Proportions of variance explained are normalized over the `k` estimated
#' scree values. Thus, for each returned result, the `k` PVE values sum to one
#' whenever the corresponding scree values have a positive finite sum.
#'
#' When `mono = TRUE`, the final monotone adjustment is a post-processing step
#' and does not change the privacy guarantee.
#'
#' For a detailed procedure and mathematical formulations,
#' refer \url{https://yejinjo0220.github.io/dppca/articles/dp_scree}.
#'
#' @return A named list. The first component, `nonprivate`, contains:
#' \itemize{
#'   \item `scree`: the usual non-private scree estimates.
#'   \item `pve`: the corresponding non-private proportions of variance
#'   explained, normalized over the `k` returned components.
#' }
#' Each requested private method is returned as an additional named component
#' (`clipped`, `pmwm`, and/or `huber`), each containing:
#' \itemize{
#'   \item `scree`: the differentially private scree estimates.
#'   \item `pve`: the corresponding private proportions of variance explained,
#'   normalized over the `k` returned components.
#' }
#' All scree and PVE vectors are named `PC1`, ..., `PCk`. The return structure is
#' the same whether one or multiple private methods are requested.
#'
#' @seealso
#' [dp_loading()] for private principal component direction estimation.
#' [clipped_control()], [pmwm_control()], and [huber_control()] for
#' method-specific tuning parameters.
#'
#' @references
#' \insertRef{dwork2014algorithmic}{dppca}
#'
#' \insertRef{ramsay2025pmw}{dppca}
#'
#' \insertRef{durfee2023unbounded}{dppca}
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
#' privacy <- privacy_control(
#'   eps = 3,
#'   delta = 1e-4,
#'   split = c(loading = 0.5, scree = 0.5)
#' )
#'
#' # Estimate private scree values using the clipped mean method.
#' set.seed(123)
#' out <- dp_scree(
#'   X,
#'   k = 2,
#'   privacy = privacy,
#'   method = "clipped",
#'   control = clipped_control(C_clip = 3)
#' )
#'
#' out$nonprivate
#' out$clipped
#'
#' # Multiple methods can be requested together by using a named control list.
#' # Loading is estimated once; each method receives the full scree budget.
#' # Change the split in privacy_control() to choose a different allocation.
#'
#' @export
dp_scree <- function(
    X,
    k,
    privacy,
    method = c("clipped", "pmwm", "huber"),
    control = NULL,
    center = TRUE,
    standardize = FALSE,
    cpp.option = FALSE,
    mono = TRUE
) {
  if (missing(method)) {
    method <- "clipped"
  } else {
    method <- match.arg(
      method,
      choices = c("clipped", "pmwm", "huber"),
      several.ok = TRUE
    )
    method <- unique(method)
  }

  privacy <- .validate_privacy_control(privacy, c("loading", "scree"))
  b <- privacy$components
  X <- as.matrix(X)

  validate_scree_inputs(
    X = X,
    k = k,
    eps = privacy$total[["eps"]],
    delta = privacy$total[["delta"]]
  )

  k <- as.integer(k)

  # For multiple methods, controls must be supplied as a named list so that
  # method-specific tuning parameters are unambiguous.
  if (length(method) > 1L && !is.null(control)) {
    if (
      !is.list(control) ||
      is.null(names(control)) ||
      any(names(control) == "")
    ) {
      stop(
        "When multiple methods are requested, `control` must be a named list, ",
        "for example list(clipped = clipped_control(...), ",
        "pmwm = pmwm_control(...), huber = huber_control(...)).",
        call. = FALSE
      )
    }
  }

  # Preprocess once, then share one private loading and score matrix.
  X_proc <- prep_matrix_for_pca(
    X = X,
    center = center,
    standardize = standardize
  )
  V_used <- dp_loading(
    X = X_proc,
    eps = b$loading[["eps"]],
    delta = b$loading[["delta"]],
    center = FALSE,
    standardize = FALSE,
    cpp.option = cpp.option
  )$private[, seq_len(k), drop = FALSE]
  Y <- X_proc %*% V_used

  # Compute the ordinary non-private PCA scree values only once.
  S_np <- stats::cov(X_proc)
  eig_np <- eigen(S_np, symmetric = TRUE, only.values = TRUE)

  scree_np <- as.numeric(eig_np$values[seq_len(k)])
  pve_np <- scree_to_pve(scree_np)

  pc_names <- paste0("PC", seq_len(k))
  names(scree_np) <- pc_names
  names(pve_np) <- pc_names

  out <- list(
    nonprivate = list(
      scree = scree_np,
      pve = pve_np
    )
  )

  # All methods use the same scores and the full scree component budget.
  # The budget is not divided across methods, supporting method comparison.
  for (m in method) {
    control_m <- if (length(method) == 1L) {
      control
    } else if (!is.null(control) && m %in% names(control)) {
      control[[m]]
    } else {
      NULL
    }

    control_m <- .merge_scree_control(m, control_m)

    result <- switch(
      m,
      clipped = dp_scree_clipped(
        Y = Y,
        eps = b$scree[["eps"]],
        delta = b$scree[["delta"]],
        C_clip = control_m$C_clip,
        mono = mono
      ),
      pmwm = dp_scree_pmwm(
        Y = Y,
        eps = b$scree[["eps"]],
        delta = b$scree[["delta"]],
        split_mode = control_m$split_mode,
        beta = control_m$beta,
        a = control_m$a,
        b = control_m$b,
        trim_const = control_m$trim_const,
        eta = control_m$eta,
        mono = mono
      ),
      huber = dp_scree_huber(
        Y = Y,
        eps = b$scree[["eps"]],
        delta = b$scree[["delta"]],
        mu0 = control_m$mu0,
        eta0 = control_m$eta0,
        T = control_m$T,
        M = control_m$M,
        k_min_m2 = control_m$k_min_m2,
        k_max_m2 = control_m$k_max_m2,
        m2_frac = control_m$m2_frac,
        mono = mono
      )
    )

    scree_m <- as.numeric(result$scree)
    pve_m <- as.numeric(result$pve)

    names(scree_m) <- pc_names
    names(pve_m) <- pc_names

    out[[m]] <- list(
      scree = scree_m,
      pve = pve_m
    )
  }

  out
}


#' Plot differentially private scree estimates
#'
#' @description
#' This function computes and visualizes scree curves for principal component
#' analysis using base R graphics. It overlays the usual non-private curve with
#' one or more differentially private estimates and acts as a plotting wrapper
#' around [dp_scree()].
#'
#' @param X A numeric matrix or data frame. Rows correspond to observations and
#'   columns correspond to variables.
#' @param k Positive integer defining the number of leading principal components
#'   to estimate. Must be an integer between `1` and the number of columns in `X`.
#' @param privacy A [privacy_control()] object with exactly the `loading` and
#'   `scree` components. For example, use
#'   `privacy_control(eps = 3, delta = 1e-4, split = c(loading = 0.5, scree = 0.5))`.
#'   The loading budget is used once; each requested scree method receives
#'   the full `scree` component budget, without division across methods.
#' @param method Scree estimation method or methods to plot. One or more of
#'   `"clipped"`, `"pmwm"`, or `"huber"`. If omitted, `"clipped"` is used.
#' @param control Optional method-specific control list, or a named list of
#'   control lists when multiple methods are requested. Use [clipped_control()],
#'   [pmwm_control()], and [huber_control()].
#' @param center A logical value indicating whether to center the columns of `X`
#'   before computing principal component directions. The default is `TRUE`.
#' @param standardize A logical value indicating whether to scale the columns of
#'   `X` by their sample standard deviations after optional centering. The
#'   default is `FALSE`.
#' @param cpp.option A logical value passed to [dp_loading()] for private
#'   loading estimation. The default is `FALSE`.
#' @param mono A logical value indicating whether to apply monotone
#'   post-processing to the private scree vector. The default is `TRUE`.
#' @param type Quantity to plot. Use `"pve"` to plot proportions of variance
#'   explained and `"scree"` to plot raw scree values. The default is `"pve"`.
#' @param plot_control Optional plotting control list created by
#'   [scree_plot_control()]. If `NULL`, default plotting settings are used.
#'
#' @details
#' This function calls [dp_scree()] once and plots its returned `nonprivate`
#' result together with each requested private method using base R graphics.
#' All private methods share the same estimated loadings and projected scores.
#' The `loading` component of `privacy` is used once; each method receives the
#' full `scree` component budget. Releasing several methods together requires
#' composition; the recorded total covers loading plus one scree method.
#'
#' The default legend labels are `"Non-private"`, `"Clipped"`, `"PMWM"`, and
#' `"Huber"`. Plot appearance, including the title, axis labels, legend
#' position, colors, line types, point symbols, and text sizes, can be changed
#' with [scree_plot_control()].
#'
#' Because the plot includes the `nonprivate` reference curve, the complete plot
#' is intended for comparison and is not itself a differentially private release.
#'
#' @return Invisibly returns the named list produced by [dp_scree()].
#'
#' @seealso
#' [dp_loading()] for private principal component direction estimation.
#' [dp_scree()] for computing non-private and differentially private scree
#' estimates.
#' [scree_plot_control()] for plot appearance.
#' [clipped_control()], [pmwm_control()], and [huber_control()] for
#' method-specific tuning parameters.
#'
#' @references
#' \insertRef{dwork2014algorithmic}{dppca}
#'
#' \insertRef{ramsay2025pmw}{dppca}
#'
#' \insertRef{durfee2023unbounded}{dppca}
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
#' privacy <- privacy_control(
#'   eps = 3,
#'   delta = 1e-4,
#'   split = c(loading = 0.5, scree = 0.5)
#' )
#'
#' # Draw a private PVE plot using the clipped mean method.
#' set.seed(123)
#' dp_scree_plot(
#'   X,
#'   k = 5,
#'   privacy = privacy,
#'   method = "clipped",
#'   control = clipped_control(C_clip = 3)
#' )
#'
#' # Customize the plot using a separate plotting control.
#' # dp_scree_plot(
#' #   X,
#' #   k = 5,
#' #   privacy = privacy,
#' #   method = "clipped",
#' #   control = clipped_control(C_clip = 3),
#' #   plot_control = scree_plot_control(
#' #     title = "PVE comparison",
#' #     xlab = "Principal Component",
#' #     legend_position = "topright"
#' #   )
#' # )
#'
#' @export
dp_scree_plot <- function(
    X,
    k,
    privacy,
    method = c("clipped", "pmwm", "huber"),
    control = NULL,
    center = TRUE,
    standardize = FALSE,
    cpp.option = FALSE,
    mono = TRUE,
    type = c("pve", "scree"),
    plot_control = NULL
) {
  if (missing(method)) {
    method <- "clipped"
  } else {
    method <- match.arg(
      method,
      choices = c("clipped", "pmwm", "huber"),
      several.ok = TRUE
    )
    method <- unique(method)
  }

  type <- match.arg(type)

  plot_control <- .merge_scree_plot_control(
    plot_control = plot_control,
    type = type
  )

  results <- dp_scree(
    X = X,
    k = k,
    privacy = privacy,
    method = method,
    control = control,
    center = center,
    standardize = standardize,
    cpp.option = cpp.option,
    mono = mono
  )

  series_names <- c("nonprivate", method)

  y_list <- lapply(series_names, function(nm) {
    if (type == "scree") {
      results[[nm]]$scree
    } else {
      results[[nm]]$pve
    }
  })
  names(y_list) <- series_names

  y_all <- unlist(y_list, use.names = FALSE)
  y_all <- y_all[is.finite(y_all)]

  if (length(y_all) == 0L) {
    ylim <- c(0, 1)
  } else {
    ylim <- range(y_all)

    if (diff(ylim) == 0) {
      pad <- if (ylim[1] == 0) 0.5 else 0.05 * abs(ylim[1])
      if (!is.finite(pad) || pad <= 0) {
        pad <- 0.5
      }
      ylim <- ylim + c(-pad, pad)
    } else {
      pad <- 0.04 * diff(ylim)
      ylim <- ylim + c(-pad, pad)
    }

    if (type == "pve") {
      ylim[1] <- max(0, ylim[1])
    }
  }

  idx <- seq_len(k)

  nm1 <- series_names[1L]
  graphics::plot(
    idx,
    y_list[[nm1]],
    type = "b",
    col = unname(plot_control$col[nm1]),
    lty = unname(plot_control$lty[nm1]),
    lwd = plot_control$lwd,
    pch = unname(plot_control$pch[nm1]),
    cex = plot_control$point_cex,
    xlab = plot_control$xlab,
    ylab = plot_control$ylab,
    main = plot_control$title,
    ylim = ylim,
    xaxt = "n",
    cex.main = plot_control$cex_main,
    cex.lab = plot_control$cex_lab,
    cex.axis = plot_control$cex_axis
  )

  graphics::axis(
    1,
    at = idx,
    labels = idx,
    cex.axis = plot_control$cex_axis
  )

  if (length(series_names) > 1L) {
    for (nm in series_names[-1L]) {
      graphics::lines(
        idx,
        y_list[[nm]],
        type = "b",
        col = unname(plot_control$col[nm]),
        lty = unname(plot_control$lty[nm]),
        lwd = plot_control$lwd,
        pch = unname(plot_control$pch[nm]),
        cex = plot_control$point_cex
      )
    }
  }

  graphics::legend(
    plot_control$legend_position,
    legend = unname(plot_control$legend_labels[series_names]),
    col = unname(plot_control$col[series_names]),
    lty = unname(plot_control$lty[series_names]),
    lwd = plot_control$lwd,
    pch = unname(plot_control$pch[series_names]),
    pt.cex = plot_control$point_cex,
    bty = plot_control$legend_bty,
    cex = plot_control$cex_legend
  )

  invisible(results)
}
