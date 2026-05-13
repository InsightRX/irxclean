# =============================================================================
# Model stability assessment across iterative removal runs
# =============================================================================

#' Check model stability across removal iterations
#'
#' Evaluates whether model parameters reached stable estimates as observations
#' are iteratively removed.  Rather than flagging every iteration where a
#' parameter drifted from baseline (which conflates expected convergence to a
#' new value with genuine instability), this function classifies each
#' parameter as:
#'
#' \describe{
#'   \item{stable}{Low variability across iterations.}
#'   \item{drifted}{Monotonic shift to a new value (few reversals); the
#'     parameter moved but settled.  This is expected behaviour when true
#'     outliers are removed.}
#'   \item{unstable}{RSE > 50\% with at least three trend reversals, or
#'     CV of the final third of iterations > 10\%.  Warrants investigation.}
#' }
#'
#' @param results A results list from \code{remove_erroneous_obs()} or
#'   \code{load_removal_results()}.
#' @param threshold_pct Numeric. RSE \% threshold above which a parameter
#'   is considered to have high variability.  Also used (as
#'   \code{threshold_pct / 4}) to assess late-iteration stability.
#'   Default 20 (i.e. 20\%).
#' @param late_window_fraction Numeric. Fraction of iterations (from the end)
#'   used to assess late-stage CV stability.  Default \code{1/3} (final
#'   third, minimum 2 iterations).
#' @param verbose Logical. Print a stability summary to the console?
#'   Default TRUE.
#'
#' @return An object of class \code{irxclean_stability} (a named list):
#' \describe{
#'   \item{param_summary}{Data frame with one row per parameter containing
#'     \code{PARAMETER}, \code{MEAN}, \code{FIRST}, \code{LAST},
#'     \code{RSE_PCT} (SD / |mean| x 100), \code{N_REVERSALS} (trend
#'     direction changes), \code{CV_LATE_PCT} (CV of the final third of
#'     iterations), \code{DRIFT_PCT} (|first - last| / |first| x 100), and
#'     \code{STATUS} (stable / drifted / unstable).}
#'   \item{parameter_drift}{Data frame of per-iteration values with
#'     \code{PCT_CHANGE} from baseline and \code{STATUS}, used for the
#'     trajectory plot.}
#'   \item{unstable_params}{Character vector of parameters classified as
#'     unstable.}
#'   \item{flagged_iterations}{Integer vector of iterations associated with
#'     unstable parameters.}
#'   \item{n_flagged_params}{Number of unstable parameters.}
#'   \item{rmse_trend}{Character: \code{"decreasing"}, \code{"stable"}, or
#'     \code{"increasing"} based on a linear trend of nRMSE over iterations.}
#'   \item{threshold_pct}{The threshold used.}
#' }
#'
#' @seealso [plot.irxclean_stability()], [stability_flag_table()],
#'   [remove_erroneous_obs()]
#'
#' @examples
#' \dontrun{
#' results   <- load_removal_results("err_rem1")
#' stability <- check_model_stability(results)
#' print(stability)
#' plot(stability)
#' stability_flag_table(stability)
#' }
#'
#' @export
check_model_stability <- function(
  results,
  threshold_pct        = 20,
  late_window_fraction = 1/3,
  verbose              = TRUE
) {
  .validate_results(results)

  par_df  <- results$par
  rmse_df <- results$rmse

  param_cols <- grep("^(THETA|OMEGA|SIGMA)", names(par_df), value = TRUE)
  n_iter     <- max(par_df$ITERATION, na.rm = TRUE)
  labels     <- results$param_labels

  # Remove fixed parameters (constant across all iterations)
  is_fixed <- vapply(par_df[, param_cols, drop = FALSE], function(x) {
    v <- x[!is.na(x)]
    length(v) == 0L || all(v == v[1L])
  }, logical(1L))
  param_cols <- param_cols[!is_fixed]

  # For OMEGAs: restrict to diagonal elements, then drop SAME via labels
  theta_cols <- param_cols[grepl("^THETA", param_cols)]
  sigma_cols <- param_cols[grepl("^SIGMA", param_cols)]
  omega_cols <- param_cols[grepl("^OMEGA", param_cols)]

  diag_omega <- omega_cols[
    grepl("^OMEGA\\.\\d+\\.\\d+$", omega_cols) &
      vapply(strsplit(omega_cols, "\\."), function(x) {
        length(x) == 3L && x[2L] == x[3L]
      }, logical(1L))
  ]
  if (length(diag_omega) == 0L) diag_omega <- omega_cols

  if (!is.null(labels) && length(labels$omegas) > 0) {
    labeled_omega <- intersect(diag_omega, names(labels$omegas))
    if (length(labeled_omega) > 0) diag_omega <- labeled_omega
  }
  param_cols <- c(theta_cols, diag_omega, sigma_cols)

  # Window for late-stage stability (final fraction of iterations, min 2)
  n_late <- max(2L, floor(n_iter * late_window_fraction))

  # -- Per-parameter summary ---------------------------------------------------
  summary_rows <- lapply(param_cols, function(param) {
    vals <- par_df[[param]][order(par_df$ITERATION)]
    vals <- vals[!is.na(vals)]
    if (length(vals) < 2L) return(NULL)

    mean_val <- mean(vals)
    rse_pct  <- if (abs(mean_val) > 0) {
      sd(vals) / abs(mean_val) * 100
    } else {
      NA_real_
    }

    # Trend reversals: sign changes in consecutive differences
    diffs <- diff(vals)
    n_rev <- if (length(diffs) >= 2L) {
      nz <- diffs[diffs != 0]
      if (length(nz) >= 2L) {
        sum(diff(sign(nz)) != 0L)
      } else {
        0L
      }
    } else {
      0L
    }

    # Late-stage CV (last n_late iterations)
    late_vals <- tail(vals, n_late)
    cv_late <- if (length(late_vals) >= 2L && abs(mean(late_vals)) > 0) {
      sd(late_vals) / abs(mean(late_vals)) * 100
    } else {
      0
    }

    # Overall drift first -> last
    drift_pct <- if (length(vals) >= 2L && abs(vals[1L]) > 0) {
      abs(vals[length(vals)] - vals[1L]) / abs(vals[1L]) * 100
    } else {
      NA_real_
    }

    # Classification:
    # unstable - high RSE with multiple reversals, OR still drifting late
    # drifted  - large shift but monotonic (settled to new value)
    # stable   - low variability throughout
    status <- if (
      (!is.na(rse_pct) && rse_pct > 50 && n_rev >= 3L) ||
        cv_late > 10
    ) {
      "unstable"
    } else if (!is.na(drift_pct) && drift_pct > threshold_pct / 2) {
      "drifted"
    } else {
      "stable"
    }

    data.frame(
      PARAMETER   = param,
      MEAN        = mean_val,
      FIRST       = vals[1L],
      LAST        = vals[length(vals)],
      RSE_PCT     = rse_pct,
      N_REVERSALS = n_rev,
      CV_LATE_PCT = cv_late,
      DRIFT_PCT   = drift_pct,
      STATUS      = status,
      stringsAsFactors = FALSE
    )
  })
  summary_rows  <- summary_rows[!vapply(summary_rows, is.null, logical(1L))]
  param_summary <- if (length(summary_rows) > 0) {
    do.call(rbind, summary_rows)
  } else {
    NULL
  }

  unstable_params <- if (!is.null(param_summary)) {
    param_summary$PARAMETER[param_summary$STATUS == "unstable"]
  } else {
    character(0L)
  }

  # Per-iteration trajectory (for plot; STATUS applied at parameter level)
  status_map <- if (!is.null(param_summary)) {
    setNames(param_summary$STATUS, param_summary$PARAMETER)
  } else {
    character(0L)
  }

  drift_rows <- lapply(param_cols, function(param) {
    vals     <- par_df[[param]]
    baseline <- vals[par_df$ITERATION == 1L]
    if (length(baseline) == 0L || is.na(baseline) || baseline == 0) {
      return(NULL)
    }
    status <- if (param %in% names(status_map)) {
      status_map[[param]]
    } else {
      "stable"
    }
    rows <- lapply(seq_len(n_iter), function(i) {
      v <- vals[par_df$ITERATION == i]
      if (length(v) == 0L || is.na(v)) return(NULL)
      data.frame(
        ITERATION  = i,
        PARAMETER  = param,
        VALUE      = v,
        BASELINE   = baseline,
        PCT_CHANGE = 100 * (v - baseline) / abs(baseline),
        STATUS     = status,
        stringsAsFactors = FALSE
      )
    })
    do.call(rbind, rows[!vapply(rows, is.null, logical(1L))])
  })
  drift_rows <- drift_rows[!vapply(drift_rows, is.null, logical(1L))]
  drift_df   <- if (length(drift_rows) > 0) {
    do.call(rbind, drift_rows)
  } else {
    NULL
  }

  flagged_iters <- if (!is.null(drift_df) && length(unstable_params) > 0) {
    sort(unique(drift_df$ITERATION[drift_df$STATUS == "unstable"]))
  } else {
    integer(0L)
  }

  n_unstable <- length(unstable_params)

  # -- Apply human-readable parameter labels ------------------------------------
  if (!is.null(labels)) {
    all_lbl <- c(
      if (length(labels$thetas) > 0) labels$thetas else character(0L),
      if (length(labels$omegas) > 0) labels$omegas else character(0L)
    )
    if (length(all_lbl) > 0) {
      .lbl <- function(x) ifelse(x %in% names(all_lbl), all_lbl[x], x)
      if (!is.null(param_summary)) {
        param_summary$PARAMETER <- .lbl(param_summary$PARAMETER)
      }
      if (!is.null(drift_df)) {
        drift_df$PARAMETER <- .lbl(drift_df$PARAMETER)
      }
      unstable_params <- .lbl(unstable_params)
      if (length(status_map) > 0) {
        status_map <- setNames(status_map, .lbl(names(status_map)))
      }
    }
  }

  # -- nRMSE trend -------------------------------------------------------------
  rmse_trend <- "stable"
  if (!is.null(rmse_df) && nrow(rmse_df) >= 3L) {
    rmse_dedup <- rmse_df[!duplicated(rmse_df$ITERATION), ]
    fit        <- stats::lm(rmse ~ ITERATION, data = rmse_dedup)
    slope      <- stats::coef(fit)[["ITERATION"]]
    rmse_trend <- if (slope < -1e-6) {
      "decreasing"
    } else if (slope > 1e-6) {
      "increasing"
    } else {
      "stable"
    }
  }

  if (verbose) {
    message(sprintf(
      "\nirxclean :: Model stability check (threshold: %g%%)\n%s",
      threshold_pct, paste(rep("-", 50), collapse = "")
    ))
    if (!is.null(param_summary)) {
      message(sprintf(
        "  Parameters assessed : %d  (fixed excluded)",
        nrow(param_summary)
      ))
      tally <- table(param_summary$STATUS)
      for (s in c("stable", "drifted", "unstable")) {
        n <- if (s %in% names(tally)) tally[[s]] else 0L
        message(sprintf("    %-10s : %d", s, n))
      }
    }
    if (n_unstable == 0L) {
      message("  All parameters reached a stable estimate.")
    } else {
      message(sprintf(
        "  Unstable: %s", paste(unstable_params, collapse = ", ")
      ))
    }
    message(sprintf("  nRMSE trend: %s", rmse_trend))
    if (rmse_trend == "increasing") {
      message(paste0(
        "  WARNING: nRMSE is increasing - ",
        "removals may be worsening model fit."
      ))
    }
  }

  structure(
    list(
      param_summary      = param_summary,
      parameter_drift    = drift_df,
      unstable_params    = unstable_params,
      flagged_iterations = flagged_iters,
      n_flagged_params   = n_unstable,
      rmse_trend         = rmse_trend,
      threshold_pct      = threshold_pct
    ),
    class = "irxclean_stability"
  )
}


#' Print method for irxclean_stability objects
#' @param x An `irxclean_stability` object.
#' @param ... Ignored.
#' @export
print.irxclean_stability <- function(x, ...) {
  cat(sprintf(
    "irxclean model stability (threshold: %g%%)\n", x$threshold_pct
  ))
  if (!is.null(x$param_summary) && nrow(x$param_summary) > 0) {
    tally <- table(x$param_summary$STATUS)
    for (s in c("stable", "drifted", "unstable")) {
      n <- if (s %in% names(tally)) tally[[s]] else 0L
      cat(sprintf("  %-10s : %d parameter(s)\n", s, n))
    }
    cat("\nParameter summary:\n")
    show_cols <- c(
      "PARAMETER", "RSE_PCT", "N_REVERSALS", "CV_LATE_PCT",
      "DRIFT_PCT", "STATUS"
    )
    show_cols <- intersect(show_cols, names(x$param_summary))
    print(
      x$param_summary[, show_cols, drop = FALSE],
      row.names = FALSE, digits = 3
    )
  } else {
    cat("  No free parameters to assess.\n")
  }
  cat(sprintf("\nnRMSE trend: %s\n", x$rmse_trend))
  invisible(x)
}


#' Plot method for irxclean_stability objects
#'
#' Produces two panels:
#' \enumerate{
#'   \item Parameter trajectories (\% change from baseline) coloured by
#'     stability status (stable = teal, drifted = grey, unstable = red).
#'   \item RSE bar chart per parameter, sorted by RSE and coloured by status.
#' }
#'
#' @param x An `irxclean_stability` object from [check_model_stability()].
#' @param ... Ignored.
#'
#' @return A list of \code{ggplot} objects: \code{$parameter_drift} and
#'   \code{$param_rse}.
#'
#' @export
plot.irxclean_stability <- function(x, ...) {
  status_colours <- c(
    "stable"   = "#003B4F",
    "drifted"  = "#5E5E5E",
    "unstable" = "#AD0000"
  )
  plots <- list()

  # -- Trajectory plot ---------------------------------------------------------
  if (!is.null(x$parameter_drift) && nrow(x$parameter_drift) > 0) {
    pd <- x$parameter_drift
    pd$STATUS <- factor(pd$STATUS, levels = c("stable", "drifted", "unstable"))

    # Invisible phantom points anchor each panel to at least -25% / +25%
    # (equivalent to 0.75x / 1.25x baseline), keeping stable panels readable
    # while allowing dramatic panels to expand freely beyond that range.
    phantom_df <- data.frame(
      PARAMETER  = rep(unique(pd$PARAMETER), each = 5L),
      PCT_CHANGE = rep(c(-25, -10, 0, 10, 25), times = length(unique(pd$PARAMETER))),
      ITERATION  = 1L,
      STATUS     = factor("stable", levels = c("stable", "drifted", "unstable"))
    )

    plots$parameter_drift <- ggplot2::ggplot(
      pd,
      ggplot2::aes(
        x      = .data$ITERATION,
        y      = .data$PCT_CHANGE,
        colour = .data$STATUS,
        group  = .data$PARAMETER
      )
    ) +
      ggplot2::geom_blank(data = phantom_df) +
      ggplot2::geom_hline(
        yintercept = 0, colour = "#5E5E5E", linewidth = 0.4
      ) +
      ggplot2::geom_line(alpha = 0.8) +
      ggplot2::geom_point() +
      ggplot2::scale_colour_manual(values = status_colours, drop = FALSE) +
      ggplot2::facet_wrap(~PARAMETER, scales = "free_y") +
      ggplot2::labs(
        title    = "Parameter trajectories across removal iterations",
        subtitle = paste0(
          "% change from iteration-1 baseline  |  ",
          "unstable = bouncing/unsettled  |  drifted = settled to new value"
        ),
        x      = "Iteration",
        y      = "% change from baseline",
        colour = "Status"
      ) +
      .clean_theme()
  }

  # -- RSE bar chart -----------------------------------------------------------
  if (!is.null(x$param_summary) && nrow(x$param_summary) > 0) {
    ps <- x$param_summary
    ps$STATUS    <- factor(
      ps$STATUS, levels = c("stable", "drifted", "unstable")
    )
    ps$PARAMETER <- factor(
      ps$PARAMETER, levels = ps$PARAMETER[order(ps$RSE_PCT)]
    )
    plots$param_rse <- ggplot2::ggplot(
      ps,
      ggplot2::aes(
        x    = .data$PARAMETER,
        y    = .data$RSE_PCT,
        fill = .data$STATUS
      )
    ) +
      ggplot2::geom_col() +
      ggplot2::geom_vline(
        xintercept = x$threshold_pct,
        linetype = 2, colour = "#AD0000"
      ) +
      ggplot2::coord_flip() +
      ggplot2::scale_fill_manual(values = status_colours, drop = FALSE) +
      ggplot2::labs(
        title    = "Parameter RSE across removal iterations",
        subtitle = sprintf(
          "RSE = SD / |mean| x 100%%  |  threshold = %g%%", x$threshold_pct
        ),
        x    = NULL,
        y    = "RSE (%)",
        fill = "Status"
      ) +
      .clean_theme()
  }

  plots
}


#' Summarise unstable parameters in a tidy table
#'
#' Returns the rows of the parameter summary for parameters classified as
#' \code{"unstable"} (i.e. those that did not reach a stable estimate).
#'
#' @param stability An `irxclean_stability` object.
#'
#' @return A data frame of unstable parameters with stability metrics.
#'
#' @export
stability_flag_table <- function(stability) {
  if (!inherits(stability, "irxclean_stability")) {
    stop("Input must be an `irxclean_stability` object.")
  }
  if (is.null(stability$param_summary)) return(data.frame())
  stability$param_summary[stability$param_summary$STATUS == "unstable", ]
}
