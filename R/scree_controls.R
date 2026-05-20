# ============================================================
# scree_controls.R
# control constructors for DP scree estimators
# ============================================================


#' Control options for clipped scree estimation
#'
#' Creates a control list for the clipped-mean scree estimator used by
#' [dp_scree()] and [dp_scree_plot()] when `method = "clipped"`.
#'
#' @param C_clip Positive clipping threshold for squared principal component
#'   scores. This value has no default because it should be chosen according to
#'   the scale of the data.
#'
#' @details
#' The clipped method estimates each scree value by clipping the squared scores
#' at `C_clip` and then applying a sensitivity-calibrated Gaussian mechanism
#' \insertCite{dwork2014algorithmic}{dppca}. Larger values of `C_clip` reduce
#' clipping bias but increase the sensitivity, and therefore the scale of
#' the privacy noise.
#'
#' @return A list of control options for `method = "clipped"`.
#'
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
#' Creates a control list for the private modified winsorized mean scree
#' estimator used by [dp_scree()] and [dp_scree_plot()] when `method = "pmwm"`.
#'
#' @param a,b Finite lower and upper search bounds supplied to the private
#'   quantile routine. The private lower and upper clipping cutoffs are searched
#'   within this range. These values have no defaults because they should be
#'   chosen on the scale of squared principal component scores.
#' @param trim_const Positive number controlling the baseline clipping level in
#'   the practical clipping proportion. This value has no default.
#' @param eta Nonnegative number controlling the expected contamination level in
#'   the practical clipping proportion. This value has no default.
#' @param beta Positive number greater than `1` defining the log-binning base
#'   used by the private quantile routine. The default is `1.001`.
#' @param split_mode A logical value indicating whether to split the sample into
#'   quantile-estimation and mean-estimation subsets. The default is `TRUE`.
#'
#' @details
#' The PMWM method privately estimates lower and upper tail cutoffs, winsorizes
#' the squared scores to those cutoffs, and then releases a noisy winsorized
#' mean. It is based on the private modified winsorized mean of
#' \insertCite{ramsay2025pmw;textual}{dppca}.
#'
#' The implementation used here is an R adaptation of the publicly available
#' Python implementation accompanying \insertCite{ramsay2025pmw;textual}{dppca}.
#' The adaptation is used for scree estimation by applying the PMWM estimator to
#' squared principal component scores.
#'
#' The PMWM scree estimator uses additional control parameters for private
#' quantile estimation and winsorization. The parameter `beta` determines the
#' spacing of the geometric search grid used by the private quantile estimator
#' and must satisfy \eqn{\beta > 1}. Smaller values of `beta` give a finer grid
#' but may increase computation.
#'
#' The bounds `a` and `b` define the lower and upper search range supplied to the
#' private quantile routine. The private lower and upper winsorization cutoffs
#' are searched within this range. These bounds should be chosen on the scale of
#' the squared principal component scores.
#'
#' The parameters `trim_const` and `eta` determine the practical clipping
#' proportion used by the modified winsorized mean. If \eqn{n_q} denotes the
#' number of observations used for private quantile estimation, the clipping
#' proportion is
#' \deqn{
#'   p = \min\left\{
#'     \max\left(\frac{\mathrm{trim\_const}}{n_q}, \eta\right),
#'     0.49
#'   \right\}.
#' }
#' Here, `trim_const / n_q` controls the baseline clipping level, while `eta`
#' gives a lower bound reflecting the expected contamination level.
#'
#' If `split_mode = TRUE`, the sample is split into two parts: one part is used
#' for private quantile estimation and the other part is used for the winsorized
#' mean step. If `split_mode = FALSE`, all observations are used in both steps.
#'
#' The parameters `a`, `b`, `trim_const`, and `eta` are intentionally not given
#' defaults. They are data- and robustness-dependent choices and should be set
#' deliberately by the user.
#'
#' @return A list of control options for `method = "pmwm"`.
#'
#' @seealso
#' [dp_scree()] for computing differentially private scree estimates using these
#' control options.
#' [dp_scree_plot()] for plotting scree estimates.
#'
#' @references
#' \insertRef{ramsay2025pmw}{dppca}
#'
#' @examples
#' pmwm_control(a = 0, b = 20, trim_const = 10, eta = 0.01)
#' pmwm_control(
#'   a = 0,
#'   b = 20,
#'   trim_const = 10,
#'   eta = 0.01,
#'   beta = 1.001,
#'   split_mode = FALSE
#' )
#'
#' @export
pmwm_control <- function(
    a, b, trim_const, eta,
    beta = 1.001,
    split_mode = TRUE
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
  if (!is.numeric(eta) || length(eta) != 1L || !is.finite(eta) || eta < 0) {
    stop("`eta` must be a nonnegative number.", call. = FALSE)
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
#' Creates a control list for the Huber-type private scree estimator used by
#' [dp_scree()] and [dp_scree_plot()] when `method = "huber"`.
#'
#' @param k_min_m2,k_max_m2 Integers defining the lower and upper dyadic bin
#'   indices used in the private second-moment scale step. The histogram searches
#'   over scale levels \eqn{2^k} for \eqn{k_{\min} \le k \le k_{\max}}. These
#'   values have no defaults because they should be chosen according to the scale
#'   of the data.
#' @param m2_frac Number in `(0, 1)` defining the fraction of the Huber scree
#'   privacy parameters allocated to the private second-moment scale step. This
#'   value has no default.
#' @param mu0 Numeric initial value for Huber noisy gradient descent. The default
#'   is `0`.
#' @param eta0 Positive number defining the fixed step size for Huber noisy
#'   gradient descent. The default is `1`.
#' @param T Optional positive integer defining the number of noisy gradient
#'   descent iterations. If `NULL`, the implementation uses
#'   \eqn{\lceil \log n \rceil}, where \eqn{n} is the number of observations.
#' @param M Optional positive integer defining the number of blocks used in the
#'   private second-moment scale step. If `NULL`, the implementation uses
#'   \eqn{\lfloor \sqrt{n} / 2 \rfloor}, where \eqn{n} is the number of
#'   observations.
#'
#' @details
#' The Huber method estimates the mean of squared principal component scores by
#' noisy gradient descent on the Huber loss. It follows the Huber-type private
#' robust mean approach of \insertCite{yu2024gaussian;textual}{dppca}.
#'
#' The method first privately estimates a scale proxy for the squared scores,
#' denoted by \eqn{m_2}. This scale proxy is then used to choose the Huber
#' robustification level for noisy gradient descent. The parameters
#' `k_min_m2`, `k_max_m2`, and `m2_frac` control this private scale-proxy step,
#' while `mu0`, `eta0`, `T`, and `M` control the subsequent noisy gradient
#' descent routine.
#'
#' The dyadic indices `k_min_m2` and `k_max_m2` define the search range for the
#' private histogram used to estimate \eqn{m_2}. The histogram searches over
#' candidate scale levels \eqn{2^k} satisfying
#' \eqn{k_{\min} \le k \le k_{\max}}. Because this range depends on the scale of
#' the squared scores, these arguments are intentionally not given defaults.
#'
#' The argument `m2_frac` determines how the Huber scree privacy parameters are
#' split between the private scale-proxy step and the noisy gradient descent
#' step. If \eqn{(\epsilon_{\mathrm{scree}}, \delta_{\mathrm{scree}})} denotes
#' the privacy parameters available for Huber scree estimation, then
#' \eqn{m2_frac \cdot (\epsilon_{\mathrm{scree}},
#' \delta_{\mathrm{scree}})} is used to privately estimate \eqn{m_2}, while
#' \eqn{(1 - m2_frac) \cdot (\epsilon_{\mathrm{scree}},
#' \delta_{\mathrm{scree}})} is used for Huber noisy gradient descent.
#'
#' The remaining parameters have default values. The default `mu0 = 0` is the
#' initial value for noisy gradient descent, and the default `eta0 = 1` is the
#' fixed step size. If `T = NULL`, the number of noisy gradient descent
#' iterations is chosen as \eqn{\lceil \log n \rceil}. If `M = NULL`, the number
#' of blocks used in the private estimator of \eqn{m_2} is chosen as
#' \eqn{\lfloor \sqrt{n} / 2 \rfloor}.
#'
#' @return A list of control options for `method = "huber"`.
#'
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
#' huber_control(
#'   k_min_m2 = -10,
#'   k_max_m2 = 10,
#'   m2_frac = 1 / 4,
#'   T = 50,
#'   M = 20
#' )
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
    pmwm = list(beta = 1.001, split_mode = TRUE),
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
