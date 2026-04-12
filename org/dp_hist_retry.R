#' @keywords internal
.dp_hist_retry <- function(
    X,
    center,
    scale.,
    dp_pca_flag,
    cpp.option,
    axes,
    eps_total,
    delta_total,
    eps_ratio,
    delta_ratio,
    inflate,
    q_frame,
    m_x,
    m_y,
    bin_method,
    mechanism,
    sampling,
    verbose = TRUE
) {
  attempt <- 1L

  repeat {
    res <- tryCatch(
      dp_hist(
        X            = X,
        center       = center,
        scale.       = scale.,
        dp_pca_flag  = dp_pca_flag,
        cpp.option   = cpp.option,
        axes         = axes,
        eps_total    = eps_total,
        delta_total  = delta_total,
        eps_ratio    = eps_ratio,
        delta_ratio  = delta_ratio,
        inflate      = inflate,
        q_frame      = q_frame,
        m_x          = m_x,
        m_y          = m_y,
        bin_method   = bin_method,
        mechanism    = mechanism,
        sampling     = sampling
      ),
      error = function(e) e
    )
    if (!inherits(res, "error")) {
      if (verbose && attempt > 1L) {
        message(
          "[dppca] dp_hist succeeded after ",
          attempt, " attempt", if (attempt > 1L) "s." else "."
        )
      }
      return(res)
    }

    msg <- conditionMessage(res)

    is_breaks_error <-
      grepl("breaks.*고유", msg) ||
      grepl("breaks are not unique", msg, ignore.case = TRUE)

    if (!is_breaks_error) {
      stop(res)
    }

    if (verbose) {
      message(
        "[dppca] Non-unique breaks in dp_hist (attempt ",
        attempt, "), retrying with new DP noise..."
      )
    }

    attempt <- attempt + 1L
  }
}

