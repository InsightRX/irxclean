# =============================================================================
# Post-removal demographic comparison plots
# =============================================================================

#' Compare demographics of subjects with removed observations vs. whole dataset
#'
#' After iterative removal, this function checks whether subjects whose
#' observations were removed are demographically representative of the full
#' dataset, and whether removals cluster at a particular position on the PK
#' curve (by time after dose).  Systematic differences may indicate that the
#' removal algorithm is biased against a specific patient population or a
#' specific part of the concentration-time profile.
#'
#' Statistical tests applied:
#' \itemize{
#'   \item Continuous covariates: two-sided Wilcoxon rank-sum test.
#'   \item Categorical covariates: chi-squared test of proportions.
#'   \item TAD bins: per-bin one-sided binomial test vs. the overall removal
#'     rate, Bonferroni-corrected.
#' }
#'
#' @param data A data frame containing the original NONMEM dataset (one row
#'   per record).
#' @param results A results list from [remove_erroneous_obs()] or
#'   [load_removal_results()].
#' @param continuous_cols Character vector of continuous covariate column names
#'   (e.g. \code{c("AGE", "WT")}).
#' @param categorical_cols Character vector of categorical covariate column
#'   names (e.g. \code{"SEX"}).
#' @param id_col Character. Subject identifier column.  Default \code{"ID"}.
#' @param time_col Character. Time column used to compute TAD when no explicit
#'   TAD column is found.  Default \code{"TIME"}.
#' @param evid_col Character. EVID column used to distinguish observations
#'   (EVID == 0) from dose records (EVID == 1).  Default \code{"EVID"}.  If
#'   not present, all rows are treated as observations.
#' @param tad_col Character or \code{NULL}.  Name of a pre-computed time-after-
#'   dose column.  \code{NULL} (default) triggers auto-detection: the function
#'   looks for columns named \code{"TAD"}, \code{"TAFD"}, or \code{"TSFD"}; if
#'   none are found it computes TAD from the dose records in \code{data}; if
#'   dosing records cannot be identified it falls back to \code{time_col}.
#' @param tad_bin_width Numeric. Width of TAD bins for the curve-position
#'   analysis.  \code{NULL} (default) auto-selects a round value from the
#'   TAD range.
#' @param alpha Numeric. Significance level.  Default \code{0.05}.
#'
#' @return An object of class \code{irxclean_demographics} (a named list):
#' \describe{
#'   \item{plots}{Named list of \code{ggplot} objects, one per covariate.}
#'   \item{panels}{List of \code{patchwork} grids, each containing up to 9
#'     covariate plots.  Significant covariates have red bold titles.}
#'   \item{tad_plot}{A \code{ggplot} showing the proportion of removed
#'     observations per TAD bin, with Bonferroni-corrected per-bin binomial
#'     tests (elevated bins in red).}
#'   \item{stats}{Data frame: \code{covariate}, \code{test}, \code{p_value},
#'     \code{significant}.}
#'   \item{removed_ids}{IDs with removed observations.}
#'   \item{n_removed_subjects}{Number of subjects with any removal.}
#' }
#'
#' @seealso [plot_demographic_comparison()], [remove_erroneous_obs()]
#'
#' @examples
#' \dontrun{
#' original_data <- read.csv("my_data.csv")
#' results       <- load_removal_results("err_rem1")
#' demo <- plot_demographic_comparison(
#'   data             = original_data,
#'   results          = results,
#'   continuous_cols  = c("AGE", "WT", "HT"),
#'   categorical_cols = "SEX"
#' )
#' demo$panels[[1]]
#' demo$tad_plot
#' demo$stats
#' }
#'
#' @export
plot_demographic_comparison <- function(
  data,
  results,
  continuous_cols  = NULL,
  categorical_cols = NULL,
  id_col           = "ID",
  time_col         = "TIME",
  evid_col         = "EVID",
  tad_col          = NULL,
  tad_bin_width    = NULL,
  alpha            = 0.05
) {
  .validate_results(results)
  if (is.null(continuous_cols) && is.null(categorical_cols)) {
    stop("Provide at least one of `continuous_cols` or `categorical_cols`.")
  }

  removed_ids <- unique(results$rem[[id_col]])

  # One row per subject
  subj_data <- data[!duplicated(data[[id_col]]), , drop = FALSE]
  subj_data$GROUP <- ifelse(
    subj_data[[id_col]] %in% removed_ids,
    "Removed", "Retained"
  )
  subj_data$GROUP <- factor(subj_data$GROUP, levels = c("Retained", "Removed"))

  all_cols <- c(continuous_cols, categorical_cols)
  bad_cols <- setdiff(all_cols, names(subj_data))
  if (length(bad_cols) > 0) {
    warning(paste("Columns not found:", paste(bad_cols, collapse = ", ")))
    continuous_cols  <- intersect(continuous_cols, names(subj_data))
    categorical_cols <- intersect(categorical_cols, names(subj_data))
  }

  grp_colours <- c("Retained" = "#003B4F", "Removed" = "#AD0000")

  plots <- list()
  stats <- list()

  # -- Continuous covariates ---------------------------------------------------
  for (cov in continuous_cols) {
    plot_data <- subj_data[, c("GROUP", cov), drop = FALSE]
    names(plot_data)[2] <- "VALUE"
    plot_data <- plot_data[!is.na(plot_data$VALUE), ]

    removed_vals  <- plot_data$VALUE[plot_data$GROUP == "Removed"]
    retained_vals <- plot_data$VALUE[plot_data$GROUP == "Retained"]

    test_res <- tryCatch(
      stats::wilcox.test(removed_vals, retained_vals, exact = FALSE),
      error = function(e) list(p.value = NA_real_)
    )
    p_val <- test_res$p.value
    sig   <- !is.na(p_val) && p_val < alpha

    p <- ggplot2::ggplot(
      plot_data,
      ggplot2::aes(x = .data$GROUP, y = .data$VALUE, fill = .data$GROUP)
    ) +
      ggplot2::geom_violin(alpha = 0.55, trim = TRUE, colour = NA) +
      ggplot2::geom_jitter(
        ggplot2::aes(colour = .data$GROUP),
        alpha = 0.3, width = 0.15, size = 1.5, show.legend = FALSE
      ) +
      ggplot2::scale_fill_manual(values = grp_colours, guide = "none") +
      ggplot2::scale_colour_manual(values = grp_colours, guide = "none") +
      ggplot2::labs(
        title    = cov,
        subtitle = sprintf(
          "Wilcoxon p = %s%s", .fmt_p(p_val), if (sig) " *" else ""
        ),
        x = NULL,
        y = cov
      ) +
      .clean_theme()

    if (sig) {
      p <- p + ggplot2::theme(
        plot.title = ggplot2::element_text(colour = "#AD0000", face = "bold")
      )
    }

    plots[[cov]] <- p
    stats[[cov]] <- data.frame(
      covariate   = cov,
      test        = "Wilcoxon rank-sum",
      p_value     = p_val,
      significant = sig,
      stringsAsFactors = FALSE
    )
  }

  # -- Categorical covariates --------------------------------------------------
  for (cov in categorical_cols) {
    plot_data <- subj_data[, c("GROUP", cov), drop = FALSE]
    names(plot_data)[2] <- "CATEGORY"
    plot_data$CATEGORY <- as.character(plot_data$CATEGORY)
    plot_data <- plot_data[!is.na(plot_data$CATEGORY), ]

    contingency <- table(plot_data$GROUP, plot_data$CATEGORY)
    test_res <- tryCatch(
      stats::chisq.test(contingency, correct = FALSE),
      error = function(e) list(p.value = NA_real_)
    )
    p_val <- test_res$p.value
    sig   <- !is.na(p_val) && p_val < alpha

    prop_df  <- as.data.frame(prop.table(contingency, margin = 1))
    names(prop_df) <- c("GROUP", "CATEGORY", "PROPORTION")
    count_df <- as.data.frame(contingency)
    names(count_df) <- c("GROUP", "CATEGORY", "COUNT")
    bar_df <- merge(prop_df, count_df, by = c("GROUP", "CATEGORY"))

    p <- ggplot2::ggplot(
      bar_df,
      ggplot2::aes(
        x    = .data$CATEGORY,
        y    = .data$PROPORTION,
        fill = .data$GROUP
      )
    ) +
      ggplot2::geom_col(position = ggplot2::position_dodge(0.8), alpha = 0.8) +
      ggplot2::scale_fill_manual(values = grp_colours) +
      ggplot2::labs(
        title    = cov,
        subtitle = sprintf(
          "Chi-squared p = %s%s", .fmt_p(p_val), if (sig) " *" else ""
        ),
        x    = cov,
        y    = "Proportion",
        fill = "Group"
      ) +
      .clean_theme()

    if (sig) {
      p <- p + ggplot2::theme(
        plot.title = ggplot2::element_text(colour = "#AD0000", face = "bold")
      )
    }

    plots[[cov]] <- p
    stats[[cov]] <- data.frame(
      covariate   = cov,
      test        = "Chi-squared",
      p_value     = p_val,
      significant = sig,
      stringsAsFactors = FALSE
    )
  }

  # -- Panel grids (≤9 plots per panel) ----------------------------------------
  n_plots <- length(plots)
  panels  <- if (n_plots > 0) {
    lapply(seq_len(ceiling(n_plots / 9L)), function(pi) {
      idx  <- seq(
        from = (pi - 1L) * 9L + 1L,
        to   = min(pi * 9L, n_plots)
      )
      n_in <- length(idx)
      patchwork::wrap_plots(
        plots[idx],
        ncol = min(3L, n_in)
      )
    })
  } else {
    list()
  }

  # -- TAD analysis ------------------------------------------------------------
  tad_plot <- .plot_tad_removals(
    data          = data,
    results       = results,
    id_col        = id_col,
    time_col      = time_col,
    evid_col      = evid_col,
    tad_col       = tad_col,
    tad_bin_width = tad_bin_width,
    alpha         = alpha
  )

  stats_df <- if (length(stats) > 0) do.call(rbind, stats) else data.frame()

  structure(
    list(
      plots              = plots,
      panels             = panels,
      tad_plot           = tad_plot,
      stats              = stats_df,
      removed_ids        = removed_ids,
      n_removed_subjects = length(removed_ids)
    ),
    class = "irxclean_demographics"
  )
}


#' Print method for irxclean_demographics objects
#' @param x An `irxclean_demographics` object.
#' @param ... Ignored.
#' @export
print.irxclean_demographics <- function(x, ...) {
  cat(sprintf(
    "irxclean demographic comparison (%d subjects with removed obs)\n",
    x$n_removed_subjects
  ))
  cat("\nStatistical tests:\n")
  print(x$stats, row.names = FALSE)
  invisible(x)
}


# Internal: resolve TAD values for each observation.
# Priority: explicit tad_col > auto-detected column > computed from doses >
#   fall back to time_col (with warning).
# Returns list(obs, tad_col, label) or NULL.
.resolve_tad <- function(data, results, id_col, time_col, evid_col, tad_col) {
  # Helper: filter to observations only
  .obs_only <- function(df) {
    if (evid_col %in% names(df)) df[df[[evid_col]] == 0, ] else df
  }

  # Helper: compute TAD per subject given dose times in the dataset
  .compute_tad_col <- function(df) {
    if (!all(c(id_col, time_col, evid_col) %in% names(df))) return(NULL)
    # Dose indicator: EVID == 1, or if no EVID==1 rows, AMT > 0
    dose_flag <- df[[evid_col]] == 1
    if (!any(dose_flag, na.rm = TRUE) && "AMT" %in% names(df)) {
      dose_flag <- df[["AMT"]] > 0 & !is.na(df[["AMT"]])
    }
    if (!any(dose_flag, na.rm = TRUE)) return(NULL)

    ids      <- unique(df[[id_col]])
    tad_list <- lapply(ids, function(id) {
      subj <- df[df[[id_col]] == id, ]
      subj <- subj[order(subj[[time_col]]), ]
      dose_times <- subj[[time_col]][dose_flag[df[[id_col]] == id]]
      obs_rows   <- subj[subj[[evid_col]] == 0, ]
      if (nrow(obs_rows) == 0 || length(dose_times) == 0) return(NULL)
      obs_rows$TAD_COMPUTED <- vapply(obs_rows[[time_col]], function(t) {
        prior <- dose_times[dose_times <= t]
        if (length(prior) == 0L) NA_real_ else t - max(prior)
      }, numeric(1L))
      obs_rows
    })
    tad_list <- tad_list[!vapply(tad_list, is.null, logical(1L))]
    if (length(tad_list) == 0L) return(NULL)
    obs <- do.call(rbind, tad_list)
    obs[!is.na(obs$TAD_COMPUTED), ]
  }

  # 1. Explicit tad_col
  if (!is.null(tad_col) && tad_col %in% names(data)) {
    obs <- .obs_only(data)
    obs <- obs[!is.na(obs[[tad_col]]), ]
    return(list(obs = obs, tad_col = tad_col, label = "Time after dose"))
  }

  # 2. Auto-detect common TAD column names
  tad_candidates <- c("TAD", "TAFD", "TSFD", "tad", "tafd")
  found_tad <- intersect(tad_candidates, names(data))
  if (length(found_tad) > 0) {
    tc  <- found_tad[1L]
    obs <- .obs_only(data)
    obs <- obs[!is.na(obs[[tc]]), ]
    return(list(obs = obs, tad_col = tc, label = "Time after dose"))
  }

  # 3. Compute TAD from dosing records
  obs_computed <- .compute_tad_col(data)
  if (!is.null(obs_computed) && nrow(obs_computed) > 0) {
    message(
      "TAD column not found - computed from dose records (EVID == 1 / AMT > 0)."
    )
    return(list(
      obs     = obs_computed,
      tad_col = "TAD_COMPUTED",
      label   = "Time after dose"
    ))
  }

  # 4. Fall back to time_col (e.g. single-dose where TIME == TAD)
  obs <- .obs_only(data)
  if (!time_col %in% names(obs) || nrow(obs) == 0) return(NULL)
  warning(
    "No TAD column found and no dose records identified. ",
    "Using ", time_col, " as a proxy for time after dose."
  )
  list(obs = obs, tad_col = time_col, label = "Time")
}


# Internal: TAD bin plot — removed observations by curve position
.plot_tad_removals <- function(
  data,
  results,
  id_col,
  time_col,
  evid_col,
  tad_col,
  tad_bin_width,
  alpha
) {
  if (!time_col %in% names(results$rem)) return(NULL)

  resolved <- .resolve_tad(
    data    = data,
    results = results,
    id_col  = id_col,
    time_col = time_col,
    evid_col = evid_col,
    tad_col  = tad_col
  )
  if (is.null(resolved)) return(NULL)

  obs     <- resolved$obs
  tc      <- resolved$tad_col
  x_label <- resolved$label

  if (!tc %in% names(obs) || nrow(obs) == 0) return(NULL)

  # Removed observation TAD values
  # Match removed obs to full data by TIME (best available link)
  rem_times <- results$rem[[time_col]]

  # Compute TAD for removed observations
  rem_ids <- if (id_col %in% names(results$rem)) {
    results$rem[[id_col]]
  } else {
    rep(NA_integer_, length(rem_times))
  }

  rem_tad <- if (tc == time_col) {
    rem_times
  } else if (tc == "TAD_COMPUTED") {
    # Re-compute TAD for removed observations using dose records
    vapply(seq_along(rem_ids), function(k) {
      if (is.na(rem_ids[k])) return(NA_real_)
      subj <- data[data[[id_col]] == rem_ids[k], ]
      if (nrow(subj) == 0 || !evid_col %in% names(subj)) return(NA_real_)
      subj     <- subj[order(subj[[time_col]]), ]
      dose_flg <- subj[[evid_col]] == 1
      if (!any(dose_flg) && "AMT" %in% names(subj)) {
        dose_flg <- subj[["AMT"]] > 0 & !is.na(subj[["AMT"]])
      }
      dose_t <- subj[[time_col]][dose_flg]
      if (length(dose_t) == 0) return(NA_real_)
      prior  <- dose_t[dose_t <= rem_times[k]]
      if (length(prior) == 0) NA_real_ else rem_times[k] - max(prior)
    }, numeric(1L))
  } else {
    # Match TAD by ID + TIME from the resolved obs data frame
    vapply(seq_along(rem_ids), function(k) {
      if (is.na(rem_ids[k])) return(NA_real_)
      rows <- obs[
        obs[[id_col]] == rem_ids[k] &
          abs(obs[[time_col]] - rem_times[k]) < 1e-6, ,
        drop = FALSE
      ]
      if (nrow(rows) == 0) NA_real_ else rows[[tc]][1L]
    }, numeric(1L))
  }

  # Drop NAs
  rem_tad <- rem_tad[!is.na(rem_tad)]
  if (length(rem_tad) == 0) return(NULL)

  obs_tad <- obs[[tc]]
  obs_tad <- obs_tad[!is.na(obs_tad)]

  max_tad <- max(obs_tad, na.rm = TRUE)

  # Auto bin width
  if (is.null(tad_bin_width)) {
    raw_w    <- max_tad / 20
    mag      <- 10^floor(log10(max(raw_w, 0.01)))
    nice     <- c(1, 2, 5, 10) * mag
    tad_bin_width <- nice[which.min(abs(nice - raw_w))]
    tad_bin_width <- max(tad_bin_width, 0.25)
  }

  # Bin assignment
  obs_bins <- floor(obs_tad / tad_bin_width) * tad_bin_width
  rem_bins <- floor(rem_tad / tad_bin_width) * tad_bin_width

  bin_total   <- as.data.frame(table(TIME_BIN = obs_bins))
  bin_total$TIME_BIN <- as.numeric(as.character(bin_total$TIME_BIN))
  names(bin_total)[2] <- "N_TOTAL"

  bin_removed <- as.data.frame(table(TIME_BIN = rem_bins))
  bin_removed$TIME_BIN <- as.numeric(as.character(bin_removed$TIME_BIN))
  names(bin_removed)[2] <- "N_REMOVED"

  bin_df <- merge(bin_total, bin_removed, by = "TIME_BIN", all.x = TRUE)
  bin_df$N_REMOVED[is.na(bin_df$N_REMOVED)] <- 0L
  bin_df$PROP_REMOVED <- bin_df$N_REMOVED / bin_df$N_TOTAL

  # Overall removal rate
  p_overall <- sum(bin_df$N_REMOVED) / sum(bin_df$N_TOTAL)

  # Per-bin Bonferroni-corrected one-sided binomial test
  n_bins  <- nrow(bin_df)
  alpha_c <- alpha / max(n_bins, 1L)

  bin_df$ELEVATED <- vapply(seq_len(n_bins), function(i) {
    n_rem <- bin_df$N_REMOVED[i]
    n_tot <- bin_df$N_TOTAL[i]
    if (n_tot < 3L || n_rem == 0L) return(FALSE)
    bt <- tryCatch(
      stats::binom.test(n_rem, n_tot, p = p_overall, alternative = "greater"),
      error = function(e) list(p.value = NA_real_)
    )
    isTRUE(!is.na(bt$p.value) && bt$p.value < alpha_c)
  }, logical(1L))

  bin_df$STATUS <- ifelse(bin_df$ELEVATED, "Elevated", "Normal")

  ggplot2::ggplot(
    bin_df,
    ggplot2::aes(
      x    = .data$TIME_BIN,
      y    = .data$PROP_REMOVED,
      fill = .data$STATUS
    )
  ) +
    ggplot2::geom_col(
      alpha = 0.85,
      width = tad_bin_width * 0.88
    ) +
    ggplot2::geom_hline(
      yintercept = p_overall,
      linetype   = 2,
      colour     = "#5E5E5E",
      linewidth  = 0.6
    ) +
    ggplot2::annotate(
      "text",
      x      = max(bin_df$TIME_BIN) * 0.98,
      y      = p_overall,
      label  = sprintf("Overall: %.1f%%", p_overall * 100),
      vjust  = -0.4,
      hjust  = 1,
      size   = 3,
      colour = "#5E5E5E"
    ) +
    ggplot2::scale_fill_manual(
      values = c("Elevated" = "#AD0000", "Normal" = "#003B4F"),
      guide  = "none"
    ) +
    ggplot2::labs(
      title    = "Removed observations by position on the PK curve",
      subtitle = sprintf(
        paste0(
          "%s bins (width = %g)  |  ",
          "red = elevated vs. overall rate (%.1f%%, Bonferroni-corrected)"
        ),
        x_label, tad_bin_width, p_overall * 100
      ),
      x = sprintf("%s (bin width = %g)", x_label, tad_bin_width),
      y = "Proportion of observations removed"
    ) +
    .clean_theme()
}


# Internal helper: format p-values
#' @keywords internal
.fmt_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) return("< 0.001")
  sprintf("%.3f", p)
}
