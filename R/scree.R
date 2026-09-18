# ============================================================-
# scree.R
# functions for differentially private scree estimation
# ============================================================-

#' Differentially private scree values
#'
#' Estimates scree values along private principal component directions using
#' means of squared differences from disjoint score pairs. An ordinary PCA
#' reference is computed only when `non_private = TRUE`.
#'
#' @param X A numeric matrix or data frame with observations in rows and
#'   variables in columns.
#' @param k Number of leading components, an integer from `1` to `ncol(X)`.
#' @param privacy A [privacy_control()] object with exactly `loading` and
#'   `scree` components. Without `V_dp`, the loading allocation is used once
#'   to estimate directions. With `V_dp`, it records the budget already used
#'   to obtain those directions and is not spent again or redistributed.
#'   Each requested method receives the full `scree` allocation.
#' @param method One or more of `"clipped"`, `"pmwm"`, and `"huber"`.
#'   If omitted, `"clipped"` is used.
#' @param control Method-specific list created by [clipped_control()],
#'   [pmwm_control()], or [huber_control()]. For several methods, supply a
#'   named list of these controls.
#' @param center Whether to center columns of `X` before projection.
#'   The default is `TRUE`.
#' @param standardize Whether to divide columns by their sample standard
#'   deviations after optional centering. The default is `FALSE`; see Details
#'   for the preprocessing assumptions of the scree calibration.
#' @param cpp.option Logical option passed to [dp_loading()] when directions
#'   must be estimated. The default is `FALSE`; ignored when `V_dp` is supplied.
#' @param mono Whether to enforce a nonnegative, nonincreasing private scree
#'   sequence by post-processing. The default is `TRUE`.
#' @param V_dp Optional full `p` by `p` private direction matrix, where
#'   `p = ncol(X)`. Its columns must be orthonormal and its rows must match
#'   the variables of `X` in order. If both row and variable names are present,
#'   they must match. Directions must correspond to the same preprocessing.
#'   Supplying this matrix skips direction estimation; its privacy provenance
#'   and previously used loading budget remain the caller's responsibility.
#' @param non_private Whether to compute and return the ordinary PCA reference.
#'   The default is `TRUE`. With `FALSE`, ordinary PCA is not computed and the
#'   `nonprivate` result is omitted.
#'
#' @details
#' Let \eqn{Y = X V} denote the shared projected score matrix after preprocessing.
#' Each method independently permutes its rows and forms
#' \eqn{m = \lfloor n/2 \rfloor} disjoint pairs, using
#' \deqn{W_{jl} = (Y_{b_j,l} - Y_{a_j,l})^2/2.}
#' One randomly selected row is unused when \eqn{n} is odd. The methods estimate
#' means of the columns of \eqn{W}, without an additional sample-size rescaling.
#' For fixed scores, the unclipped pair mean has the usual sample variance as
#' its expectation over random pairing; an individual pairing need not equal
#' that variance.
#'
#' The supported methods differ in their treatment of \eqn{W}:
#' \itemize{
#'   \item `"clipped"` clips entries at `C_clip` and releases the vector of
#'   column means with Gaussian noise calibrated jointly across components.
#'   \item `"pmwm"` estimates two cutoffs per component using the one-sided
#'   unbounded quantile routine of
#'   \insertCite{durfee2023unbounded;textual}{dppca}, with public lower bound
#'   zero. The cutoffs are truncated to the public bounds from
#'   [pmwm_control()]. Each of the \eqn{2k} quantiles receives
#'   \eqn{\epsilon/(4k)}; the vector of winsorized means receives
#'   \eqn{\epsilon/2} and all of the method's \eqn{\delta}.
#'   \item `"huber"` first estimates a scale for each component using
#'   `m2_frac * epsilon / k`. This scale step uses only epsilon. The vector
#'   noisy-gradient procedure uses the remaining epsilon and all delta,
#'   with Gaussian calibration across components and iterations. At least
#'   eight observations are required; see [huber_control()].
#' }
#'
#' The data are preprocessed once. Directions are estimated once with
#' [dp_loading()] unless `V_dp` is supplied, and the first `k` directions
#' produce the score matrix shared by the requested methods. Loading is not
#' re-estimated within a scree method.
#'
#' The paired-score sensitivity calculations condition on the directions
#' and use fixed-size replacement adjacency. Common centering cancels in
#' the pair differences. Variable scaling must be fixed or public, or have
#' its privacy cost accounted for separately; sample standard deviations
#' computed from the input are not covered by these sensitivity calculations.
#' Supplying `V_dp` does not by itself establish its privacy guarantee.
#'
#' Multiple methods receive the same full scree allocation for comparison.
#' Releasing their results together requires composition across all methods,
#' counting loading once. The total recorded in `privacy` covers loading plus
#' one scree method, not a joint release of several methods. Repeated calls
#' likewise require separate accounting for newly released estimates.
#'
#' The optional `nonprivate` result uses ordinary PCA eigenvalues, rather than
#' variances along the DP directions. It is a comparison reference and is not
#' a private release. PVE is normalized over the `k` returned components and
#' sums to one when their scree values have a positive finite sum. Monotone
#' adjustment accesses only the estimated scree values.
#'
#' @return A named list containing one entry per requested method (`clipped`,
#'   `pmwm`, and/or `huber`). Each entry has `scree` and `pve` vectors named
#'   `PC1`, ..., `PCk`. If `non_private = TRUE`, the first entry, `nonprivate`,
#'   contains the corresponding ordinary PCA `scree` and `pve` reference.
#'   PVE is normalized over the `k` returned components in every entry.
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
    mono = TRUE,
    V_dp = NULL,
    non_private = TRUE
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
  validate_logical_value(non_private, "non_private")

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
  # V_used <- dp_loading(
  #   X = X_proc,
  #   eps = b$loading[["eps"]],
  #   delta = b$loading[["delta"]],
  #   center = FALSE,
  #   standardize = FALSE,
  #   cpp.option = cpp.option
  # )$private[, seq_len(k), drop = FALSE]
  V_used <- .get_dp_directions(
    X_proc,
    b$loading,
    cpp.option,
    V_dp
  )[, seq_len(k), drop = FALSE]
  Y <- X_proc %*% V_used

  pc_names <- paste0("PC", seq_len(k))
  out <- list()

  if (non_private) {
    scree_np <- .nonprivate_pca(
      X_proc
    )$eigenvalues[seq_len(k)]

    out$nonprivate <- list(
      scree = stats::setNames(
        scree_np,
        pc_names
      ),
      pve = stats::setNames(
        scree_to_pve(scree_np),
        pc_names
      )
    )
  }

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
#' Computes and plots one or more scree estimates using base R graphics.
#' The ordinary PCA curve is included only when `non_private = TRUE`.
#'
#' @inheritParams dp_scree
#' @param type Quantity to plot: `"pve"` for proportions of variance explained
#'   or `"scree"` for scree values. The default is `"pve"`.
#' @param plot_control Optional list from [scree_plot_control()]. If `NULL`,
#'   default plotting settings are used.
#'
#' @details
#' Calls [dp_scree()] once and plots the requested method results, together
#' with the optional `nonprivate` reference. All methods share one estimated
#' or supplied direction matrix and one projected score matrix. Each method
#' receives the full scree allocation; the composition and preprocessing
#' conditions described in [dp_scree()] also apply here.
#'
#' Use [scree_plot_control()] to change titles, axes, colors, line types,
#' point symbols, legend position, and text sizes. The usual labels are
#' `"Non-private"`, `"Clipped"`, `"PMWM"`, and `"Huber"` for included curves.
#' A plot with the non-private reference is intended for comparison and is
#' not itself a differentially private release.
#'
#' @return Invisibly returns the list produced by [dp_scree()]. The
#'   `nonprivate` component is present only when `non_private = TRUE`.
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
    plot_control = NULL,
    V_dp = NULL,
    non_private = TRUE
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
    mono = mono,
    V_dp = V_dp,
    non_private = non_private
  )

  series_names <- c(
    if (non_private) "nonprivate",
    method
  )

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
