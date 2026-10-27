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
#' @param X A numeric matrix or data frame with observations in rows and variables
#'   in columns. Non-numeric data-frame columns are excluded with a warning;
#'   numeric values must be finite. At least two observations and one numeric
#'   feature are required. If directions must be estimated, the pairwise-separation
#'   condition of [dp_loading()] also applies.
#' @param k Number of leading components to estimate; an integer from `1` to
#'   the number of retained numeric features.
#' @param privacy A [privacy_control()] object allocating budget to `loading`
#'   and `scree` (equal shares by default). Each requested scree method receives
#'   the full scree allocation.
#' @param method One or more scree estimators: `"clipped"`, `"pmwm"`, or `"huber"`
#'   (`"clipped"` by default).
#' @param control Tuning parameters from [clipped_control()], [pmwm_control()],
#'   or [huber_control()]. For multiple methods, supply a named list of
#'   method-specific controls; required tuning values must be supplied.
#' @param cpp.option Whether to use the compiled spherical-Kendall implementation
#'   (`TRUE` by default) when estimating private directions; ignored if
#'   `V_dp` is supplied.
#' @param mono Whether to apply nonincreasing isotonic regression (`TRUE` by
#'   default). Negative scree estimates are truncated to zero in either case.
#' @param V_dp Optional `p` by `p` private PC direction matrix (`NULL` by default),
#'   where `p` is the number of numeric PCA features. Columns must be
#'   orthonormal and rows must follow the feature order in `X`; matching
#'   row names are required when both sets of names exist.
#' @param non_private Whether to compute and store ordinary PCA reference results
#'   (`FALSE` by default). The reference contains ordinary PCA scree and PVE.
#'
#' @details
#' Let `Y = X V` be scores along the private directions. For each method,
#' `m = floor(n / 2)` disjoint pairs of rows are formed and the component-wise
#' values are
#' \deqn{W_{jl} = (Y_{b_j,l} - Y_{a_j,l})^2/2.}
#' An odd final observation is unused. The methods privately estimate the
#' column means of `W`, corresponding to variance along private directions
#' rather than necessarily to sample covariance eigenvalues.
#'
#' The estimators differ as follows:
#' \itemize{
#'   \item `"clipped"`: Clip the paired values at `C_clip`, then perturb their
#'   mean vector with Gaussian noise.
#'   \item `"pmwm"`: Estimate private lower and upper quantiles, winsorize
#'   the paired values, and perturb their mean vector with Gaussian noise.
#'   Quantile calls share half of the method's epsilon; the mean receives
#'   the remaining half and all delta.
#'   \item `"huber"`: Estimate private component scales and apply noisy Huber
#'   gradient updates. The scale step uses `m2_frac` of epsilon; the updates
#'   use the remainder and all delta. At least eight observations are required.
#' }
#'
#' Directions are estimated once unless `V_dp` is supplied. The input is not
#' automatically centered or standardized; sensitivity arguments use fixed-size
#' replacement adjacency and condition on the directions. Data-dependent
#' scaling requires separate privacy analysis. Sample size, schema, component
#' selection, and tuning parameters are assumed public.
#'
#' Each requested scree method uses the full scree allocation. Releasing
#' several methods together, or releasing results from repeated calls, requires
#' composition beyond the stored budget for a single method. The optional
#' ordinary PCA reference is not private.
#'
#' Scree values are truncated at zero, optionally made nonincreasing, and
#' normalized over the returned `k` components for PVE. A zero total yields
#' zero PVE. These steps are post-processing and do not reorder the directions.
#'
#' @return A `dp_scree` list containing one entry per requested method
#'   (`clipped`, `pmwm`, and/or `huber`). Each entry has `scree` and `pve` vectors named
#'   `PC1`, ..., `PCk`. The `nonprivate` entry holds the optional ordinary PCA
#'   reference (or `NULL` when omitted), and `metadata` stores input dimensions,
#'   methods and the recorded privacy budget.
#'   PVE is normalized over the `k` returned components in every entry.
#'   Use `summary()` for scree, PVE and cumulative PVE tables, and `plot()`
#'   to draw stored estimates without a new randomized computation.
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
#' privacy <- privacy_control(eps = 3, delta = 1e-4)
#'
#' # Estimate clipped private scree values and an ordinary PCA comparison.
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
    cpp.option = TRUE,
    mono = TRUE,
    V_dp = NULL,
    non_private = FALSE
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
  b <- privacy$budgets
  X_proc <- prep_matrix_for_pca(X)

  validate_scree_inputs(
    X = X_proc,
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

  # Use the same constructors as dp_pca() to validate each method's controls.
  # Validate before estimating directions so bad controls consume no DP draws.
  controls <- stats::setNames(lapply(method, function(m) {
    control_m <- if (length(method) == 1L) {
      control
    } else if (!is.null(control) && m %in% names(control)) {
      control[[m]]
    } else {
      NULL
    }
    control_m <- .merge_scree_control(m, control_m)
    control_constructor <- switch(
      m,
      clipped = clipped_control,
      pmwm = pmwm_control,
      huber = huber_control
    )
    do.call(control_constructor, control_m)
  }), method)

  # Validate input once, then share one private loading and score matrix.
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
    control_m <- controls[[m]]

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

  nonprivate <- out$nonprivate
  out$nonprivate <- NULL
  if (is.null(nonprivate)) {
    out["nonprivate"] <- list(NULL)
  } else {
    out <- c(list(nonprivate = nonprivate), out)
  }
  out$metadata <- list(
    n = nrow(X_proc), p = ncol(X_proc), k = k,
    methods = method, privacy = privacy
  )
  structure(out, class = "dp_scree")
}

#' Print private scree results
#'
#' @param x A `dp_scree` object.
#' @param ... Further arguments, reserved for compatibility with the generic.
#' @return Invisibly returns `x`.
#' @rdname dp_scree
#' @export
print.dp_scree <- function(x, ...) {
  methods <- setdiff(names(x), c("nonprivate", "metadata"))
  k <- length(x[[methods[1L]]]$scree)
  cat("Differentially Private Scree Estimates\n")
  cat("Data:", x$metadata$n, "observations x",
      x$metadata$p, "variables\n")
  cat("Components:", x$metadata$p, "PC directions;", k, "scree estimates\n")
  cat("Methods: scree =", paste(methods, collapse = ", "), "\n")
  budget <- x$metadata$privacy$total
  cat("Total privacy budget",
      if (length(methods) > 1L) " (per method)" else "",
      ": eps = ", format(budget[["eps"]], digits = 6),
      "; delta = ", format(budget[["delta"]], digits = 6), "\n", sep = "")
  cat("Non-private reference:",
      if (is.null(x$nonprivate)) "not stored" else "included", "\n")
  invisible(x)
}

#' Summarize stored private scree estimates
#'
#' @param object A `dp_scree` object.
#' @param ... Further arguments, reserved for compatibility with the generic.
#' @return A `summary.dp_scree` object containing one table per stored method.
#'   Each table gives scree values, PVE and cumulative PVE. PVE is normalized
#'   over the fitted components, as in [dp_scree()]. Printing the summary
#'   shows component importance by method and the privacy allocations.
#' @rdname dp_scree
#' @export
summary.dp_scree <- function(object, ...) {
  results <- object[setdiff(names(object), c("metadata", "nonprivate"))]
  if (!is.null(object$nonprivate)) {
    results$nonprivate <- object$nonprivate
  }
  tables <- lapply(results, function(result) {
    data.frame(
      scree = result$scree,
      pve = result$pve,
      cumulative_pve = cumsum(result$pve)
    )
  })
  structure(tables, class = "summary.dp_scree",
            n = object$metadata$n, p = object$metadata$p,
            privacy = object$metadata$privacy)
}

#' @param x A `summary.dp_scree` object.
#' @param digits Number of significant digits used for printing.
#' @rdname dp_scree
#' @export
print.summary.dp_scree <- function(x, digits = 4, ...) {
  methods <- setdiff(names(x), c("nonprivate", "metadata"))
  k <- nrow(x[[methods[1L]]])
  cat("DP Scree Estimates Summary\n")
  cat("Data:", attr(x, "n"), "observations x",
      attr(x, "p"), "variables\n")
  cat("Components:", attr(x, "p"), "PC directions;", k, "scree estimates\n")
  cat("Methods: scree =", paste(methods, collapse = ", "), "\n")

  for (method in methods) {
    result <- x[[method]]
    importance <- rbind(
      Eigenvalue = result$scree,
      PVE = result$pve,
      `Cumulative PVE` = result$cumulative_pve
    )
    colnames(importance) <- rownames(result)
    cat("\nDP component importance",
        if (length(methods) > 1L) paste0(" (", method, ")") else "",
        ":\n", sep = "")
    print(importance, digits = digits, ...)
  }

  cat("\nProportions are normalized over", k, "estimated components.\n")

  privacy <- attr(x, "privacy")
  budget <- rbind(
    Loading = privacy$budgets$loading,
    Scree = privacy$budgets$scree,
    Total = privacy$total
  )
  cat("\nAllocated privacy budget",
      if (length(methods) > 1L) " (per method)" else "",
      ":\n", sep = "")
  print(budget, digits = digits, ...)
  if (length(methods) > 1L) {
    cat("\nAll scree methods use the same privately estimated PC directions.\n")
  }
  cat("\nNon-private reference:",
      if (is.null(x$nonprivate)) "not stored" else "included", "\n")
  invisible(x)
}

#' Plot stored private scree estimates
#'
#' @param x A `dp_scree` object.
#' @param type One or both quantities to plot: `"pve"` for proportions of
#'   variance explained and `"scree"` for scree values. The default draws both
#'   side by side. Supply one value to draw only that quantity.
#' @param plot_control Optional list from [scree_plot_control()]. If `NULL`,
#'   default plotting settings are used.
#' @param ... Further arguments, reserved for compatibility with the generic.
#' @return Invisibly returns `x` after drawing the requested plot(s).
#' @details Uses stored estimates only, so repeated plotting consumes no new
#'   privacy budget. The ordinary reference is non-private comparison output
#'   and is displayed automatically when present in the fitted result.
#' @export
plot.dp_scree <- function(
    x,
    type = c("pve", "scree"),
    plot_control = NULL,
    ...
) {
  if (length(list(...))) {
    stop("Use `type` or `plot_control` to customize the plot.",
         call. = FALSE)
  }
  type <- unique(match.arg(type, c("pve", "scree"), several.ok = TRUE))
  non_private <- !is.null(x$nonprivate)
  if (length(type) > 1L) {
    old_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par), add = TRUE)
    graphics::par(mfrow = c(1L, length(type)))
  }

  for (kind in type) {
    control <- .merge_scree_plot_control(plot_control, kind)

    series_names <- setdiff(names(x), "metadata")
    if (!non_private) {
      series_names <- setdiff(series_names, "nonprivate")
    }

    y_list <- lapply(series_names, function(nm) {
      if (kind == "scree") {
        x[[nm]]$scree
      } else {
        x[[nm]]$pve
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

      if (kind == "pve") {
        ylim[1] <- max(0, ylim[1])
      }
    }

    idx <- seq_along(y_list[[1L]])

    nm1 <- series_names[1L]
    graphics::plot(
      idx,
      y_list[[nm1]],
      type = "b",
      col = unname(control$col[nm1]),
      lty = unname(control$lty[nm1]),
      lwd = control$lwd,
      pch = unname(control$pch[nm1]),
      cex = control$point_cex,
      xlab = control$xlab,
      ylab = control$ylab,
      main = control$title,
      ylim = ylim,
      xaxt = "n",
      cex.main = control$cex_main,
      cex.lab = control$cex_lab,
      cex.axis = control$cex_axis
    )

    graphics::axis(
      1,
      at = idx,
      labels = idx,
      cex.axis = control$cex_axis
    )

    if (length(series_names) > 1L) {
      for (nm in series_names[-1L]) {
        graphics::lines(
          idx,
          y_list[[nm]],
          type = "b",
          col = unname(control$col[nm]),
          lty = unname(control$lty[nm]),
          lwd = control$lwd,
          pch = unname(control$pch[nm]),
          cex = control$point_cex
        )
      }
    }

    graphics::legend(
      control$legend_position,
      legend = unname(control$legend_labels[series_names]),
      col = unname(control$col[series_names]),
      lty = unname(control$lty[series_names]),
      lwd = control$lwd,
      pch = unname(control$pch[series_names]),
      pt.cex = control$point_cex,
      bty = control$legend_bty,
      cex = control$cex_legend
    )
  }

  invisible(x)
}

# Scree estimation controls -------------------------------------------------

#' Control options for clipped scree estimation
#'
#' Creates a control list for `method = "clipped"` in [dp_scree()]
#' or the scree configuration of [dp_pca()].
#'
#' @param C_clip Required positive public clipping threshold for squared pair
#'   differences. Larger values reduce clipping bias but increase Gaussian
#'   noise; data-dependent selection requires separate privacy accounting.
#'
#' @details
#' The method clips the paired quantities at `C_clip` and adds Gaussian
#' noise to their vector of column means. For \eqn{m=\lfloor n/2\rfloor}
#' pairs and \eqn{k} components, the calibration uses vector sensitivity
#' \eqn{\sqrt{k}\,C_{clip}/m} under the conditions in [dp_scree()]. Larger
#' thresholds reduce clipping bias but increase the noise scale.
#'
#' @return A control list for `method = "clipped"`.
#' @seealso
#' [dp_scree()] for using clipped scree estimation.
#' [plot.dp_scree()] for plotting private scree estimates.
#'
#' @references
#' \insertRef{dwork2014algorithmic}{dppca}
#'
#' @examples
#' clipped_control(C_clip = 3)
#'
#' @export
clipped_control <- function(C_clip) {
  if (missing(C_clip)) {
    stop("`C_clip` must be supplied.", call. = FALSE)
  }
  if (!is.numeric(C_clip) || length(C_clip) != 1L || !is.finite(C_clip) || C_clip <= 0) {
    stop("`C_clip` must be a positive number.", call. = FALSE)
  }

  list(C_clip = C_clip)
}

#' Control options for private modified winsorized scree estimation
#'
#' Creates a control list for `method = "pmwm"` in [dp_scree()]
#' or the scree configuration of [dp_pca()].
#'
#' @param a,b Required finite public bounds for private quantile cutoffs,
#'   with `a < b`, on the scale of squared pair differences.
#' @param trim_const Required positive public trimming constant.
#' @param eta Required public tuning value in `[0, 0.5)`.
#' @param beta Geometric quantile-search base greater than `1`
#'   (`1.001` by default).
#' @param split_mode Whether to split score pairs between quantile and mean
#'   estimation (`FALSE` by default).
#'
#' @details
#' PMWM estimates lower and upper private quantiles of squared pair
#' differences, truncates the cutoffs to `[a, b]`, winsorizes the data, and
#' perturbs the mean vector. The cutoffs, not the original data, are bounded.
#' If the estimated upper cutoff is smaller than the lower, the cutoffs are
#' reset to `a` and `b`.
#'
#' With `k` components, each of `2k` private quantiles receives
#' `epsilon / (4k)`; Gaussian mean estimation receives the remaining half
#' of epsilon and all delta. If `n_q` is the number of paired rows used for
#' quantile estimation, the trimming proportion is
#' \deqn{p=\min\{\max(\mathrm{trim\_const}/n_q,\eta),0.49\}.}
#' With `split_mode = TRUE`, the pair rows are split between quantile and
#' mean estimation (at least four observations are required); otherwise all
#' rows are reused under the same privacy allocation. Tunings must be public
#' or selected through an accounted private procedure.
#'
#' @return A control list for `method = "pmwm"`.
#' @seealso
#' [dp_scree()] for computing differentially private scree estimates using these
#' control options.
#' [plot.dp_scree()] for plotting scree estimates.
#'
#' @references
#' \insertRef{ramsay2025pmw}{dppca}
#'
#' \insertRef{durfee2023unbounded}{dppca}
#'
#' @examples
#' pmwm_control(a = 0, b = 20, trim_const = 10, eta = 0.01)
#'
#' @export
pmwm_control <- function(
    a, b, trim_const, eta,
    beta = 1.001,
    split_mode = FALSE
) {
  if (missing(a) || missing(b) || missing(trim_const) || missing(eta)) {
    stop(
      "`a`, `b`, `trim_const`, and `eta` must be supplied.",
      call. = FALSE
    )
  }
  if (!is.numeric(a) || length(a) != 1L || !is.finite(a)) {
    stop("`a` must be a finite number.", call. = FALSE)
  }
  if (!is.numeric(b) || length(b) != 1L || !is.finite(b)) {
    stop("`b` must be a finite number.", call. = FALSE)
  }
  if (a >= b) {
    stop("`a` must be smaller than `b`.", call. = FALSE)
  }
  if (
    !is.numeric(trim_const) || length(trim_const) != 1L ||
    !is.finite(trim_const) || trim_const <= 0
  ) {
    stop("`trim_const` must be a positive number.", call. = FALSE)
  }
  if (!is.numeric(eta) || length(eta) != 1L ||
      !is.finite(eta) || eta < 0 || eta >= 0.5) {
    stop(
      "`eta` must be a number in `[0, 0.5)`.",
      call. = FALSE
    )
  }
  if (!is.numeric(beta) || length(beta) != 1L || !is.finite(beta) || beta <= 1) {
    stop("`beta` must be a number greater than 1.", call. = FALSE)
  }
  validate_logical_value(split_mode, "split_mode")

  list(
    beta = beta,
    a = a,
    b = b,
    trim_const = trim_const,
    eta = eta,
    split_mode = split_mode
  )
}

#' Control options for Huber scree estimation
#'
#' Creates a control list for `method = "huber"` in [dp_scree()]
#' or the scree configuration of [dp_pca()].
#'
#' @param k_min_m2,k_max_m2 Required public integer bounds for the dyadic
#'   scale grid, with `k_min_m2 < k_max_m2`.
#' @param m2_frac Required fraction in `(0, 1)` of epsilon used for the
#'   private scale step.
#' @param mu0 Public initial mean estimate (`0` by default).
#' @param eta0 Positive public gradient step size (`1` by default).
#' @param T Number of gradient iterations; `NULL` (default) uses
#'   `ceiling(log(m))`, where `m = floor(n / 2)`.
#' @param M Number of scale-estimation blocks; `NULL` (default) uses
#'   `floor(sqrt(m / 2))`, capped by available input pairs.
#'
#' @details
#' The Huber method privately estimates a robust scale for each component
#' using a noisy dyadic histogram, then performs projected noisy gradient
#' updates to estimate the mean of squared pair differences. The scale proxy
#' sets the Huber threshold; it is not a raw second moment. At least eight
#' observations are required.
#'
#' The private scale step uses `m2_frac * epsilon / k` for each component
#' and no delta. Joint gradient updates use the remaining epsilon and all
#' delta, calibrated across components and iterations. Each update is projected
#' to nonnegative values; optional monotone adjustment follows estimation.
#' Tuning parameters must be public or have separate privacy accounting.
#'
#' @return A control list for `method = "huber"`.
#' @seealso
#' [dp_scree()] for computing differentially private scree estimates using these
#' control options.
#' [plot.dp_scree()] for plotting scree estimates.
#'
#' @references
#' \insertRef{yu2024gaussian}{dppca}
#'
#' @examples
#' huber_control(k_min_m2 = -10, k_max_m2 = 10, m2_frac = 1 / 4)
#'
#' @export
huber_control <- function(
    k_min_m2, k_max_m2, m2_frac,
    mu0 = 0,
    eta0 = 1,
    T = NULL,
    M = NULL
) {
  if (missing(k_min_m2) || missing(k_max_m2) || missing(m2_frac)) {
    stop(
      "`k_min_m2`, `k_max_m2`, and `m2_frac` must be supplied.",
      call. = FALSE
    )
  }
  if (
    !is.numeric(k_min_m2) || length(k_min_m2) != 1L ||
    !is.finite(k_min_m2) || k_min_m2 != as.integer(k_min_m2)
  ) {
    stop("`k_min_m2` must be a finite integer.", call. = FALSE)
  }
  if (
    !is.numeric(k_max_m2) || length(k_max_m2) != 1L ||
    !is.finite(k_max_m2) || k_max_m2 != as.integer(k_max_m2)
  ) {
    stop("`k_max_m2` must be a finite integer.", call. = FALSE)
  }
  if (k_min_m2 >= k_max_m2) {
    stop("`k_min_m2` must be smaller than `k_max_m2`.", call. = FALSE)
  }
  if (
    !is.numeric(m2_frac) || length(m2_frac) != 1L ||
    !is.finite(m2_frac) || m2_frac <= 0 || m2_frac >= 1
  ) {
    stop("`m2_frac` must be a number in `(0, 1)`.", call. = FALSE)
  }
  if (!is.numeric(mu0) || length(mu0) != 1L || !is.finite(mu0)) {
    stop("`mu0` must be a finite number.", call. = FALSE)
  }
  if (!is.numeric(eta0) || length(eta0) != 1L || !is.finite(eta0) || eta0 <= 0) {
    stop("`eta0` must be a positive number.", call. = FALSE)
  }
  if (!is.null(T)) {
    validate_positive_integer(T, "T")
  }
  if (!is.null(M)) {
    validate_positive_integer(M, "M")
  }

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

#' Default control options for scree estimation
#'
#' @param method Scree estimation method.
#'
#' @return A method-specific control list containing only non-data-dependent
#'   defaults.
#'
#' @noRd
.default_scree_control <- function(method) {
  switch(
    method,
    clipped = list(),
    pmwm = list(beta = 1.001, split_mode = FALSE),
    huber = list(mu0 = 0, eta0 = 1, T = NULL, M = NULL)
  )
}

#' Merge user-supplied scree controls with defaults
#'
#' @param method Scree estimation method.
#' @param control Optional user-supplied control list.
#'
#' @return A complete method-specific control list.
#'
#' @noRd
.merge_scree_control <- function(method, control) {
  default <- .default_scree_control(method)

  if (is.null(control)) {
    control <- default
  } else {
    if (!is.list(control)) {
      stop(
        "`control` must be a control list created by `clipped_control()`, ",
        "`pmwm_control()`, or `huber_control()`.",
        call. = FALSE
      )
    }
    control <- utils::modifyList(default, control, keep.null = TRUE)
  }

  if (method == "clipped") {
    .validate_required_control(control, "C_clip", "clipped_control(C_clip = ...)")
  }

  if (method == "pmwm") {
    .validate_required_control(
      control,
      c("a", "b", "trim_const", "eta"),
      "pmwm_control(a = ..., b = ..., trim_const = ..., eta = ...)"
    )
  }

  if (method == "huber") {
    .validate_required_control(
      control,
      c("k_min_m2", "k_max_m2", "m2_frac"),
      "huber_control(k_min_m2 = ..., k_max_m2 = ..., m2_frac = ...)"
    )
  }

  control
}

#' Validate that required control entries are present
#'
#' @param control A control list.
#' @param required Required names.
#' @param constructor Constructor call to suggest in the error message.
#'
#' @return Invisibly returns `TRUE`.
#'
#' @noRd
.validate_required_control <- function(control, required, constructor) {
  missing_names <- required[
    !required %in% names(control) |
      vapply(control[required], is.null, logical(1))
  ]

  if (length(missing_names) > 0L) {
    stop(
      "Missing required control parameter(s): `",
      paste(missing_names, collapse = "`, `"),
      "`. Provide `control = ", constructor, "`.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}


