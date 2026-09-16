# Shared privacy budgets for PCA, scree and score estimation.

#' Allocate privacy budgets for shared PCA estimation
#'
#' @param eps,delta Total privacy parameters. `eps` must be positive and
#'   finite, and `delta` must lie in `(0, 1)`. Leave both `NULL` when
#'   supplying absolute component budgets instead.
#' @param split Named positive proportions that sum to one. The same
#'   proportions apply to `eps` and `delta`. Include `loading` and either
#'   `scree`, `score`, or both. The default allocates one third to each of
#'   `loading`, `scree` and `score`, as required by [dp_pca()]. For
#'   [dp_scree()] use, for example, `c(loading = 0.5, scree = 0.5)`; for
#'   [dp_score()] use `c(loading = 0.5, score = 0.5)`.
#' @param loading,scree,score Optional absolute component budgets, each a
#'   named `c(eps = ..., delta = ...)` vector. Supply `loading` with either
#'   `scree`, `score`, or both, without total budgets or a `split` argument.
#'   Each component must have positive finite `eps` and `delta` in `(0, 1)`;
#'   the component deltas must sum to less than one. Absolute budgets may
#'   allocate different proportions of `eps` and `delta`.
#'
#' @return A `dppca_privacy_control` list with `total`, a named
#'   `c(eps, delta)` vector, and `components`, a list of the supplied
#'   component budget vectors in `loading`, `scree`, `score` order.
#'
#' @details
#' Use the component set required by the estimation function: `loading`
#'   and `scree` for scree functions, `loading` and `score` for score
#'   functions, or all three for [dp_pca()]. An estimation function rejects
#'   a different component set rather than discarding or redistributing a
#'   budget. Loading is estimated once and shared within that function call.
#'
#' When comparing several scree or score methods, each method receives the
#' corresponding full component allocation. The allocation is not divided
#' across methods, so releasing several methods together requires composition.
#' The returned totals sum each named component once. This constructor records
#' allocations; it does not establish a privacy guarantee for preprocessing
#' or the complete estimation procedure.
#'
#' @examples
#' # Integrated PCA: equal loading, scree and score allocations.
#' privacy_control(eps = 6, delta = 1e-4)
#'
#' # Scree estimation: equal loading and scree allocations.
#' privacy_control(
#'   eps = 6,
#'   delta = 1e-4,
#'   split = c(loading = 0.5, scree = 0.5)
#' )
#'
#' # Score estimation: a different allocation for its two components.
#' privacy_control(
#'   eps = 6,
#'   delta = 1e-4,
#'   split = c(loading = 0.4, score = 0.6)
#' )
#'
#' # Absolute budgets, without total budgets or split.
#' privacy_control(
#'   loading = c(eps = 2, delta = 2e-5),
#'   scree = c(eps = 4, delta = 8e-5)
#' )
#'
#' @export
privacy_control <- function(
    eps = NULL,
    delta = NULL,
    split = c(loading = 1 / 3, scree = 1 / 3, score = 1 / 3),
    loading = NULL,
    scree = NULL,
    score = NULL
) {
  allowed_components <- c("loading", "scree", "score")

  check_components <- function(
      component_names,
      argument
  ) {
    if (
      !is.character(component_names) ||
      !(length(component_names) %in% c(2L, 3L)) ||
      anyNA(component_names) || anyDuplicated(component_names) ||
      any(!component_names %in% allowed_components) ||
      !"loading" %in% component_names
    ) {
      stop(
        "`", argument, "` must name loading and scree, loading and score, ",
        "or all three components exactly once.",
        call. = FALSE
      )
    }

    allowed_components[allowed_components %in% component_names]
  }

  check_number <- function(
      value,
      argument,
      upper = Inf
  ) {
    if (
      !is.numeric(value) || is.complex(value) || length(value) != 1L ||
      !is.finite(value) || value <= 0 || value >= upper
    ) {
      requirement <- if (is.finite(upper)) {
        " must be a finite number in (0, 1)."
      } else {
        " must be a positive finite number."
      }
      stop("`", argument, "`", requirement, call. = FALSE)
    }

    as.numeric(value)
  }

  check_budget <- function(
      budget,
      component
  ) {
    if (
      !is.numeric(budget) || is.complex(budget) ||
      !is.null(dim(budget)) || length(budget) != 2L ||
      is.null(names(budget)) || anyNA(names(budget)) ||
      anyDuplicated(names(budget)) ||
      !setequal(names(budget), c("eps", "delta"))
    ) {
      stop(
        "`", component, "` must be a named c(eps = ..., delta = ...) vector.",
        call. = FALSE
      )
    }

    c(
      eps = check_number(budget[["eps"]], paste0(component, "$eps")),
      delta = check_number(budget[["delta"]], paste0(component, "$delta"), 1)
    )
  }

  direct_budgets <- list(loading = loading, scree = scree, score = score)
  supplied <- !vapply(direct_budgets, is.null, logical(1))

  if (any(supplied)) {
    if (!is.null(eps) || !is.null(delta) || !missing(split)) {
      stop(
        "Use either eps and delta with split, or absolute component budgets; ",
        "do not combine the two modes.",
        call. = FALSE
      )
    }

    components <- check_components(names(direct_budgets)[supplied], "budgets")
    budgets <- direct_budgets[components]
  } else {
    eps <- check_number(eps, "eps")
    delta <- check_number(delta, "delta", 1)

    if (
      !is.numeric(split) || is.complex(split) || !is.null(dim(split)) ||
      any(!is.finite(split)) || any(split <= 0) ||
      !is.finite(sum(split)) || abs(sum(split) - 1) > 1e-8
    ) {
      stop("`split` must contain positive finite proportions summing to one.", call. = FALSE)
    }
    components <- check_components(names(split), "split")
    split <- split[components]

    budgets <- lapply(split, function(share) {
      c(eps = eps * share, delta = delta * share)
    })
  }

  # Recheck computed shares too, including zero values from underflow.
  for (component in components) {
    budgets[[component]] <- check_budget(budgets[[component]], component)
  }

  total_eps <- sum(vapply(budgets, function(budget) budget[["eps"]], numeric(1)))
  total_delta <- sum(vapply(budgets, function(budget) budget[["delta"]], numeric(1)))
  total_eps <- check_number(total_eps, "total eps")

  if (!is.finite(total_delta) || total_delta >= 1) {
    stop("The sum of component deltas must be less than one.", call. = FALSE)
  }

  structure(
    list(
      total = c(eps = total_eps, delta = total_delta),
      components = budgets
    ),
    class = "dppca_privacy_control"
  )
}

#' Validate the privacy object required by an estimation function
#'
#' @param privacy Output of [privacy_control()].
#' @param components Exact component names required by the caller.
#'
#' @return A revalidated `dppca_privacy_control` object with totals recomputed
#'   from all its component budgets. No components are dropped or reallocated.
#' @noRd
.validate_privacy_control <- function(
    privacy,
    components
) {
  allowed_components <- c("loading", "scree", "score")
  if (
    !is.character(components) || !(length(components) %in% c(2L, 3L)) ||
    anyNA(components) || anyDuplicated(components) ||
    any(!components %in% allowed_components) || !"loading" %in% components
  ) {
    stop(
      "`components` must request loading and scree, loading and score, ",
      "or all three components exactly once.",
      call. = FALSE
    )
  }
  components <- allowed_components[allowed_components %in% components]

  share <- if (length(components) == 3L) "1 / 3" else "0.5"
  split_example <- paste(paste0(components, " = ", share), collapse = ", ")
  example <- paste0(
    "privacy_control(eps = ..., delta = ..., split = c(", split_example, "))"
  )
  guidance <- paste0(
    "This function requires exactly ", paste(components, collapse = ", "),
    " budgets. Supply ", example, "."
  )

  if (!is.list(privacy) || !inherits(privacy, "dppca_privacy_control")) {
    stop(guidance, call. = FALSE)
  }

  budgets <- privacy[["components"]]
  if (
    !is.list(budgets) || is.null(names(budgets)) || anyNA(names(budgets)) ||
    anyDuplicated(names(budgets)) || !setequal(names(budgets), components)
  ) {
    stop(guidance, call. = FALSE)
  }

  tryCatch(
    do.call(privacy_control, budgets),
    error = function(error) {
      stop(conditionMessage(error), " ", guidance, call. = FALSE)
    }
  )
}
