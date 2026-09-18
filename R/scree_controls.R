# ============================================================
# scree_controls.R
# control constructors for DP scree estimators
# ============================================================


#' Control options for clipped scree estimation
#'
#' Creates a control list for `method = "clipped"` in [dp_scree()],
#' [dp_scree_plot()], or the scree configuration of [dp_pca()].
#'
#' @param C_clip Positive public clipping threshold for squared pair differences
#'   \eqn{(Y_{b_j,l} - Y_{a_j,l})^2/2}. It has no default and should be chosen
#'   for this scale without an unaccounted use of private data.
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
#' [dp_scree_plot()] for plotting private scree estimates.
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
#' Creates a control list for `method = "pmwm"` in [dp_scree()],
#' [dp_scree_plot()], or the scree configuration of [dp_pca()].
#'
#' @param a,b Finite public post-processing bounds for the private cutoffs,
#'   with `a < b`. They have no defaults and refer to the scale of the squared
#'   pair differences used in [dp_scree()].
#' @param trim_const Positive public constant controlling the baseline trimming
#'   proportion. It has no default.
#' @param eta Public lower bound on the practical trimming proportion, in
#'   `[0, 0.5)`. It has no default.
#' @param beta Geometric search-grid base greater than `1` for private quantile
#'   estimation. The default is `1.001`.
#' @param split_mode Whether to split the pair rows into quantile and mean
#'   subsets. The default is `FALSE`, which reuses all pairs in both steps.
#'
#' @details
#' The PMWM construction follows
#' \insertCite{ramsay2025pmw;textual}{dppca}. It estimates lower and upper
#' quantiles of the squared pair differences using the one-sided unbounded
#' routine of \insertCite{durfee2023unbounded;textual}{dppca}, with public
#' lower bound zero. Both private cutoffs are truncated to `[a, b]`; the
#' paired quantities are then winsorized and their mean vector is perturbed
#' with Gaussian noise.
#'
#' For \eqn{k} components, each of the \eqn{2k} pure-DP quantile calls receives
#' \eqn{\epsilon/(4k)}. The vector mean step receives \eqn{\epsilon/2} and
#' all of the method's delta budget. Gaussian noise is calibrated jointly
#' across components, using the private cutoff widths. See [dp_scree()] for
#' preprocessing and release-accounting conditions.
#'
#' Let \eqn{n_q} be the number of pair rows used for quantile estimation.
#' The trimming proportion is
#' \deqn{p=\min\{\max(\mathrm{trim\_const}/n_q,\eta),0.49\}.}
#' With `split_mode = TRUE`, the \eqn{\lfloor n/2\rfloor} pair rows are
#' randomly split between quantile and mean estimation, requiring at least
#' four observations. With `FALSE`, all pair rows are reused in both steps.
#'
#' Smaller `beta` gives a finer grid but may require more search steps. The
#' parameters without defaults should be chosen as public tuning inputs;
#' choosing them from private data requires separate privacy accounting.
#'
#' @return A control list for `method = "pmwm"`.
#' @seealso
#' [dp_scree()] for computing differentially private scree estimates using these
#' control options.
#' [dp_scree_plot()] for plotting scree estimates.
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
  .validate_control_logical_value(split_mode, "split_mode")

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
#' Creates a control list for `method = "huber"` in [dp_scree()],
#' [dp_scree_plot()], or the scree configuration of [dp_pca()].
#'
#' @param k_min_m2,k_max_m2 Public integer bounds for the dyadic scale histogram,
#'   with `k_min_m2 < k_max_m2`. The candidate scales are \eqn{2^r} for
#'   \eqn{k_{min}\le r\le k_{max}}. Both bounds must be supplied.
#' @param m2_frac Fraction in `(0, 1)` of the method's epsilon budget assigned
#'   to the pure-DP scale step. It must be supplied.
#' @param mu0 Finite public initial value for noisy gradient descent.
#'   The default is `0`.
#' @param eta0 Positive public fixed step size. The default is `1`.
#' @param T Optional positive integer number of iterations. If `NULL`, uses
#'   \eqn{\lceil\log m\rceil}, where \eqn{m=\lfloor n/2\rfloor} is the
#'   number of disjoint score pairs.
#' @param M Optional positive integer number of scale-estimation blocks. If
#'   `NULL`, uses \eqn{\lfloor\sqrt{m/2}\rfloor}. The scale helper caps this
#'   value at the number of available pairs of its input values.
#'
#' @details
#' The method applies a Huber-type robust mean procedure to squared pair
#' score differences, following the approach of
#' \insertCite{yu2024gaussian;textual}{dppca}. It requires at least four score
#' pairs, or eight observations. For each component, a private scale proxy
#' is computed by pairing the squared score differences again, forming
#' block medians, and selecting a dyadic scale with a noisy histogram.
#' The scale proxy determines a Huber clipping threshold.
#'
#' The scale step uses `m2_frac * epsilon / k` for each of the `k` components
#' and consumes no delta. Joint noisy gradient descent uses
#' `(1 - m2_frac) * epsilon` and all of the method's delta. Its Gaussian
#' calibration accounts for all components and iterations. Every iterate is
#' projected to the nonnegative orthant; optional monotone adjustment is
#' applied after the last iteration.
#'
#' `k_min_m2`, `k_max_m2`, and `M` control scale estimation, while `mu0`,
#' `eta0`, and `T` control gradient descent. Tuning choices should be public
#' or have their privacy costs accounted for separately. See [dp_scree()]
#' for the preprocessing and release-accounting conditions.
#'
#' @return A control list for `method = "huber"`.
#' @seealso
#' [dp_scree()] for computing differentially private scree estimates using these
#' control options.
#' [dp_scree_plot()] for plotting scree estimates.
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
    .validate_control_positive_integer(T, "T")
  }
  if (!is.null(M)) {
    .validate_control_positive_integer(M, "M")
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

#' Extract a control list for a plot method
#'
#' @param control Optional control list or named list of control lists.
#' @param method Scree estimation method.
#'
#' @return A method-specific control list or `NULL`.
#'
#' @noRd
.extract_plot_control <- function(control, method) {
  if (is.null(control)) {
    return(NULL)
  }

  if (!is.list(control)) {
    stop(
      "`control` must be a control list, or a named list of control lists ",
      "for `dp_scree_plot()`.",
      call. = FALSE
    )
  }

  if (!is.null(control[[method]]) && is.list(control[[method]])) {
    return(control[[method]])
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

#' Validate a logical value
#'
#' @param x Object to validate.
#' @param arg Argument name.
#'
#' @return Invisibly returns `TRUE`.
#'
#' @noRd
.validate_control_logical_value <- function(x, arg) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop("`", arg, "` must be `TRUE` or `FALSE`.", call. = FALSE)
  }

  invisible(TRUE)
}

#' Validate a positive integer value
#'
#' @param x Object to validate.
#' @param arg Argument name.
#'
#' @return Invisibly returns `TRUE`.
#'
#' @noRd
.validate_control_positive_integer <- function(x, arg) {
  if (
    !is.numeric(x) || length(x) != 1L || !is.finite(x) ||
    x <= 0 || x != as.integer(x)
  ) {
    stop("`", arg, "` must be a positive integer.", call. = FALSE)
  }

  invisible(TRUE)
}
