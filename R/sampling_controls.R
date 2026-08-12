# ============================================================
# sampling_controls.R
# Control options for synthetic score sampling
# ============================================================

#' Control options for synthetic score sampling
#'
#' Creates a control list for generating synthetic score points from a released
#' private score histogram. Sampling is used only as post-processing of the
#' private histogram; it does not alter the histogram privacy mechanism.
#'
#' @param method Character vector specifying the sampling construction.
#'   `"center"` selects a histogram bin according to its private probability
#'   and samples around that bin center using Gaussian jitter.
#'   `"uniform"` selects a histogram bin according to its private probability,
#'   samples uniformly inside the selected bin, and then adds Gaussian jitter.
#'   The default is `"uniform"`. Both methods can be requested with
#'   `c("center", "uniform")`.
#' @param sample_size Optional positive integer giving the synthetic sample
#'   size for each requested sampling construction. For [dp_score_plot()], this
#'   is the number of pooled synthetic score points. For
#'   [dp_score_plot_group()], this is the total number of synthetic points
#'   across all groups and is allocated in proportion to the observed group
#'   sizes. If `NULL`, the total number of input observations is used. This
#'   default therefore assumes that the sample size and group sizes are not
#'   sensitive.
#' @param bandwidth_scale Gaussian smoothing scale. If `NULL`, method-specific
#'   defaults are used: `0.5` for `"center"` and `0.25` for `"uniform"`.
#'   A single nonnegative number applies the same scale to all requested
#'   sampling methods. A named numeric vector can specify method-specific
#'   values, for example `c(center = 0.5, uniform = 0.2)`.
#'
#' @details
#' Let \eqn{\widetilde p_k} denote the released private probability of
#' histogram bin \eqn{B_k}, and let \eqn{c_k} be its center. Both sampling
#' constructions first draw a bin index
#' \deqn{J \sim \mathrm{Categorical}(\widetilde p_1,\ldots,\widetilde p_m).}
#'
#' For `method = "center"`, a synthetic score is generated as
#' \deqn{S^* = c_J + Z.}
#'
#' For `method = "uniform"`, a point is first drawn uniformly from the selected
#' bin and then jittered:
#' \deqn{U \mid J=k \sim \mathrm{Uniform}(B_k), \qquad S^* = U + Z.}
#'
#' In both cases,
#' \deqn{
#' Z \sim N_2(0,H), \qquad
#' H = c^2 \mathrm{diag}(\Delta_x^2,\Delta_y^2),
#' }
#' where \eqn{\Delta_x} and \eqn{\Delta_y} are the histogram bin widths and
#' \eqn{c} is `bandwidth_scale`.
#'
#' The sampling step uses only the released private histogram and fixed
#' histogram geometry. Consequently, when its tuning parameters are specified
#' without consulting the original score data, it is a post-processing
#' operation and requires no additional privacy budget.
#'
#' @return A named list with components `method`, `sample_size`, and
#'   `bandwidth_scale`.
#'
#' @examples
#' sampling_control()
#'
#' sampling_control(
#'   method = c("center", "uniform"),
#'   sample_size = 1000,
#'   bandwidth_scale = c(center = 0.5, uniform = 0.2)
#' )
#'
#' @export
sampling_control <- function(
    method = "uniform",
    sample_size = NULL,
    bandwidth_scale = NULL
) {
  method <- unique(
    match.arg(
      method,
      choices = c("center", "uniform"),
      several.ok = TRUE
    )
  )

  .validate_sampling_size(sample_size)
  scales <- .resolve_sampling_bandwidth(
    method = method,
    bandwidth_scale = bandwidth_scale
  )

  out <- list(
    method = method,
    sample_size = sample_size,
    bandwidth_scale = scales
  )

  class(out) <- c("dppca_sampling_control", "list")
  out
}


#' Default sampling control
#'
#' @return A sampling-control list using the package defaults.
#' @noRd
.default_sampling_control <- function() {
  sampling_control(
    method = "uniform",
    sample_size = NULL,
    bandwidth_scale = NULL
  )
}


#' Resolve and validate a sampling control
#'
#' @param control Optional sampling-control list.
#' @param n Number of input observations.
#'
#' @return A validated list with a resolved positive integer `sample_size`.
#' @noRd
.resolve_sampling_control <- function(control, n) {
  if (is.null(control)) {
    control <- .default_sampling_control()
  }

  if (!is.list(control)) {
    stop(
      "`sampling_control` must be a list created by `sampling_control()`.",
      call. = FALSE
    )
  }

  required <- c("method", "sample_size", "bandwidth_scale")
  missing <- setdiff(required, names(control))
  if (length(missing) > 0L) {
    stop(
      "`sampling_control` is missing component(s): ",
      paste(missing, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  method <- unique(
    match.arg(
      control$method,
      choices = c("center", "uniform"),
      several.ok = TRUE
    )
  )

  .validate_sampling_size(control$sample_size)

  if (
    !is.numeric(n) ||
    length(n) != 1L ||
    !is.finite(n) ||
    n < 1 ||
    n != as.integer(n)
  ) {
    stop("`n` must be a positive integer.", call. = FALSE)
  }

  sample_size <- if (is.null(control$sample_size)) {
    as.integer(n)
  } else {
    as.integer(control$sample_size)
  }

  scales <- .resolve_sampling_bandwidth(
    method = method,
    bandwidth_scale = control$bandwidth_scale
  )

  list(
    method = method,
    sample_size = sample_size,
    bandwidth_scale = scales
  )
}


#' Resolve method-specific bandwidth scales
#'
#' @param method Requested sampling methods.
#' @param bandwidth_scale User-supplied scale specification.
#'
#' @return Named numeric vector with one scale per requested method.
#' @noRd
.resolve_sampling_bandwidth <- function(method, bandwidth_scale) {
  defaults <- c(
    center = 0.5,
    uniform = 0.25
  )

  method <- unique(
    match.arg(
      method,
      choices = names(defaults),
      several.ok = TRUE
    )
  )

  if (is.null(bandwidth_scale)) {
    return(defaults[method])
  }

  if (
    !is.numeric(bandwidth_scale) ||
    length(bandwidth_scale) < 1L ||
    anyNA(bandwidth_scale) ||
    any(!is.finite(bandwidth_scale)) ||
    any(bandwidth_scale < 0)
  ) {
    stop(
      "`bandwidth_scale` must be `NULL` or contain finite nonnegative numbers.",
      call. = FALSE
    )
  }

  if (length(bandwidth_scale) == 1L) {
    return(
      stats::setNames(
        rep(as.numeric(bandwidth_scale), length(method)),
        method
      )
    )
  }

  nm <- names(bandwidth_scale)
  if (is.null(nm) || any(nm == "")) {
    stop(
      "When `bandwidth_scale` has length greater than one, it must be a named ",
      "vector with names `center` and/or `uniform`.",
      call. = FALSE
    )
  }

  if (anyDuplicated(nm)) {
    stop("`bandwidth_scale` must not contain duplicated names.", call. = FALSE)
  }

  invalid <- setdiff(nm, names(defaults))
  if (length(invalid) > 0L) {
    stop(
      "Unsupported `bandwidth_scale` name(s): ",
      paste(invalid, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  missing_method <- setdiff(method, nm)
  if (length(missing_method) > 0L) {
    stop(
      "`bandwidth_scale` must provide a value for each requested sampling ",
      "method. Missing: ",
      paste(missing_method, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  as.numeric(bandwidth_scale[method]) |>
    stats::setNames(method)
}


#' Validate synthetic sample size
#'
#' @param sample_size Sample-size specification.
#'
#' @return Invisibly returns `TRUE`.
#' @noRd
.validate_sampling_size <- function(sample_size) {
  if (is.null(sample_size)) {
    return(invisible(TRUE))
  }

  if (
    !is.numeric(sample_size) ||
    length(sample_size) != 1L ||
    !is.finite(sample_size) ||
    sample_size < 1 ||
    sample_size != as.integer(sample_size)
  ) {
    stop(
      "`sample_size` must be `NULL` or a positive integer.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}
