# Shared privacy budgets for PCA, scree and score estimation.

#' Allocate privacy budgets for shared PCA estimation
#'
#' @param eps,delta Total privacy parameters. `eps` must be positive and
#'   finite, and `delta` must lie in `(0, 1)`. Leave both `NULL` when
#'   supplying absolute component budgets instead.
#' @param split Named positive proportions that sum to one. The same
#'   proportions apply to `eps` and `delta`. Include `loading` alone, or
#'   `loading` with `scree`, `score`, or both. When `split` is omitted,
#'   the fitting function automatically divides the total budget equally
#'   among its required components: `loading` for [dp_loading()],
#'   `loading` and `scree` for [dp_scree()], `loading` and `score` for
#'   [dp_score()], or all three for [dp_pca()]. An explicitly supplied
#'   `split` is used exactly as specified, without automatic redistribution.
#' @param loading,scree,score Optional absolute component budgets, each a
#'   named `c(eps = ..., delta = ...)` vector. Supply `loading` alone or
#'   with `scree`, `score`, or both, without total budgets or a `split` argument.
#'   Each component must have positive finite `eps` and `delta` in `(0, 1)`;
#'   the component deltas must sum to less than one. Absolute budgets may
#'   allocate different proportions of `eps` and `delta`.
#'
#' @return A `dppca_privacy_control` list with `total`, a named
#'   `c(eps, delta)` vector, and `budgets`, a list of the supplied
#'   component budget vectors in `loading`, `scree`, `score` order.
#'
#' @details
#' The estimation function resolves the required component budgets:
#' `loading` for [dp_loading()], `loading` plus `scree` for [dp_scree()],
#' `loading` plus `score` for [dp_score()], and all three for [dp_pca()].
#' With total `eps` and `delta` and no explicit `split`, the constructor
#' initially shows equal thirds, but the fitting function resolves the actual
#' equal shares. An explicit split or absolute budget must match the required
#' components.
#'
#' With supplied `V_dp`, the loading allocation records the prior direction
#' release; it is not spent again or redistributed. Score epsilon is divided
#' 35% to the private frame and 65% to the histogram; all score delta goes to
#' the histogram. When multiple scree or score estimators are requested, each
#' receives its respective full allocation. Releasing them together requires
#' privacy composition; the recorded totals are allocations, not cumulative
#' privacy expenditure. Preprocessing and supplied directions require their
#' own privacy justification.
#'
#' @examples
#' # With no split, the fitting function resolves equal shares:
#' # loading only (100%), scree/score (50%-50%), or PCA (one-third each).
#' privacy_control(eps = 6, delta = 1e-4)
#'
#' # Alternatively, specify the loading-only absolute budget.
#' privacy_control(loading = c(eps = 6, delta = 1e-4))
#'
#' # Explicit splits customize allocation for scree or score.
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
      !(length(component_names) %in% 1:3) ||
      anyNA(component_names) || anyDuplicated(component_names) ||
      any(!component_names %in% allowed_components) ||
      !"loading" %in% component_names
    ) {
      stop(
        "`", argument, "` must name loading alone, loading and scree, ",
        "loading and score, or all three components exactly once.",
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
  implicit_split <- !any(supplied) && missing(split)

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

  result <- structure(
    list(
      total = c(eps = total_eps, delta = total_delta),
      budgets = budgets
    ),
    class = "dppca_privacy_control"
  )

  # Distinguish an implicit three-way default from an explicit allocation.
  if (implicit_split) attr(result, "implicit_split") <- TRUE
  result
}

#' Validate the privacy object required by an estimation function
#'
#' @param privacy Output of [privacy_control()].
#' @param components Exact component names required by the caller.
#'
#' @return A revalidated `dppca_privacy_control` object. An implicit
#'   default split is resolved equally across the requested components;
#'   explicit component allocations are validated without redistribution.
#' @noRd
.validate_privacy_control <- function(
    privacy,
    components
) {
  allowed_components <- c("loading", "scree", "score")
  if (
    !is.character(components) || !(length(components) %in% 1:3) ||
    anyNA(components) || anyDuplicated(components) ||
    any(!components %in% allowed_components) || !"loading" %in% components
  ) {
    stop(
      "`components` must request loading alone, loading and scree, ",
      "loading and score, or all three components exactly once.",
      call. = FALSE
    )
  }
  components <- allowed_components[allowed_components %in% components]

  share <- if (length(components) == 3L) "1 / 3" else "0.5"
  split_example <- paste(paste0(components, " = ", share), collapse = ", ")
  example <- if (length(components) == 1L) {
    "privacy_control(loading = c(eps = ..., delta = ...))"
  } else {
    paste0(
      "privacy_control(eps = ..., delta = ..., split = c(", split_example, "))"
    )
  }
  guidance <- paste0(
    "This function requires exactly ", paste(components, collapse = ", "),
    " budget component(s). Use privacy_control(eps = ..., delta = ...) ",
    "for automatic equal allocation, or supply ", example, "."
  )

  if (!is.list(privacy) || !inherits(privacy, "dppca_privacy_control")) {
    stop(guidance, call. = FALSE)
  }

  budgets <- privacy[["budgets"]]

  # Only an omitted split is context-dependent. The constructor retains
  # its three-way default plus this marker until a fitting function knows
  # which components it needs. Validate the original amounts first.
  if (identical(attr(privacy, "implicit_split"), TRUE)) {
    if (
      !is.list(budgets) || is.null(names(budgets)) || anyNA(names(budgets)) ||
      anyDuplicated(names(budgets)) || !setequal(names(budgets), allowed_components)
    ) {
      stop(guidance, call. = FALSE)
    }
    checked <- tryCatch(
      do.call(privacy_control, budgets),
      error = function(error) {
        stop(conditionMessage(error), " ", guidance, call. = FALSE)
      }
    )
    equal_split <- stats::setNames(
      rep(1 / length(components), length(components)), components
    )
    return(privacy_control(
      eps = checked$total[["eps"]],
      delta = checked$total[["delta"]],
      split = equal_split
    ))
  }

  # Explicit splits and absolute component budgets must match exactly.
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
