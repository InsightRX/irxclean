# =============================================================================
# Gross (model-independent) data quality assessment
# =============================================================================

#' Assess data quality before model-dependent cleaning
#'
#' Performs a series of model-independent checks on a NONMEM-format dataset:
#'
#' \enumerate{
#'   \item **Covariate outliers** - flags subjects whose covariate values fall
#'     outside \eqn{Q_1 - k \cdot IQR} or \eqn{Q_3 + k \cdot IQR} (where
#'     \eqn{k} = `iqr_multiplier`).
#'   \item **Concentration TAD-bin outliers** - observations whose DV falls
#'     outside IQR-based bounds within their time-after-dose (TAD) bin.
#'     TAD is read from a TAD column (auto-detected or specified via
#'     `tad_col`); if absent it is computed from dose records (EVID == 1).
#'   \item **Concentration elevations without a preceding dose** - flags
#'     consecutive observation pairs within a subject where the concentration
#'     rises by more than `conc_increase_threshold`-fold with no intervening
#'     dose event.
#'   \item **Missing data summary** - tabulates `NA` rates by column.
#' }
#'
#' @param data A data frame in NONMEM format.
#' @param columns Named list of column labels.  Must include at minimum:
#'   `ID`, `TIME`, `DV`, `EVID`, `MDV`.
#' @param covariate_cols Character vector of continuous covariate column names
#'   to screen for outliers.  `NULL` skips this check.
#' @param categorical_cols Character vector of categorical covariate column
#'   names.  Currently used only for completeness in the missing-data summary.
#' @param iqr_multiplier Numeric. IQR multiplier for outlier bounds, applied to
#'   both covariate outlier detection and concentration bin outlier detection.
#'   Default `3`.
#' @param tad_col Character or \code{NULL}.  Name of a pre-computed
#'   time-after-dose column.  \code{NULL} (default) triggers auto-detection:
#'   looks for columns named \code{"TAD"}, \code{"TAFD"}, or \code{"TSFD"};
#'   if none found, computes TAD from dose records (EVID == 1 / AMT > 0);
#'   if no dose records exist, falls back to \code{TIME}.
#' @param tad_bin_width Numeric. Width of TAD bins used when checking for
#'   concentration outliers within each TAD window.  \code{NULL} (default)
#'   auto-selects a round value from the TAD range.  Bins with fewer than
#'   4 observations are skipped.
#' @param conc_increase_threshold Numeric. Minimum fold-increase in
#'   concentration (without intervening dose) to flag as suspicious.
#'   Default `1.5`.
#' @param conc_floor Numeric. Minimum value for the preceding concentration
#'   when evaluating fold-increases.  Pairs where the earlier DV is at or
#'   below this value are skipped (ratio is undefined / noise-dominated near
#'   the assay LOQ).  Set to your assay LLOQ.  Default `0` (no floor).
#' @param verbose Logical. Print a summary of findings?  Default `TRUE`.
#'
#' @return An object of class \code{irxclean_quality} (a named list) with
#'   elements:
#' \describe{
#'   \item{covariate_outliers}{Data frame of flagged covariate values, or
#'     \code{NULL} if none detected / not requested.}
#'   \item{concentration_bin_outliers}{Data frame of observations whose DV is
#'     extreme relative to other observations in the same TAD bin, or
#'     \code{NULL} if none detected.  Columns: ID, TIME, TAD, TAD_BIN, DV,
#'     LOWER_BOUND, UPPER_BOUND, N_IN_BIN.}
#'   \item{conc_scatter_data}{Data frame of all observation rows with columns
#'     \code{ID}, \code{TIME}, \code{TAD}, \code{DV}, and \code{OUTLIER}
#'     (logical).  Used for the DV vs TAD scatter plot.}
#'   \item{concentration_flags}{Data frame of suspicious concentration
#'     elevations without an intervening dose, or \code{NULL} if none detected.}
#'   \item{missing_summary}{Data frame with columns \code{column},
#'     \code{n_missing}, and \code{pct_missing}.}
#'   \item{n_subjects}{Number of unique subjects.}
#'   \item{n_observations}{Number of observation rows (EVID==0, MDV==0).}
#' }
#'
#' @seealso [plot.irxclean_quality()], [remove_erroneous_obs()]
#'
#' @examples
#' \dontrun{
#' dat <- read.csv("my_data.csv")
#' qc  <- assess_data_quality(
#'   dat,
#'   covariate_cols   = c("AGE", "WT", "HT"),
#'   categorical_cols = "SEX"
#' )
#' print(qc)
#' plot(qc)
#' }
#'
#' @export
assess_data_quality <- function(
    data,
    columns = list(
      ID   = "ID",
      TIME = "TIME",
      DV   = "DV",
      EVID = "EVID",
      MDV  = "MDV"
    ),
    covariate_cols          = NULL,
    categorical_cols        = NULL,
    iqr_multiplier          = 3,
    tad_col                 = NULL,
    tad_bin_width           = NULL,
    conc_increase_threshold = 1.5,
    conc_floor              = 0,
    verbose                 = TRUE
) {
  # -- Validate columns ---------------------------------------------------------
  required <- c("ID", "TIME", "DV", "EVID", "MDV")
  missing_cols <- setdiff(required, names(columns))
  if (length(missing_cols) > 0) {
    stop(paste("columns list missing keys:", paste(missing_cols, collapse = ", ")))
  }
  for (key in required) {
    if (!columns[[key]] %in% names(data)) {
      stop(sprintf("Column '%s' ('%s') not found in data.", key, columns[[key]]))
    }
  }

  id_col   <- columns$ID
  evid_col <- columns$EVID
  mdv_col  <- columns$MDV

  n_subj <- length(unique(data[[id_col]]))
  n_obs  <- sum(data[[evid_col]] == 0 & data[[mdv_col]] == 0, na.rm = TRUE)

  # -- Missing data -------------------------------------------------------------
  miss_summary <- .summarise_missing(data)

  # -- Covariate outliers -------------------------------------------------------
  cov_outliers <- NULL
  if (!is.null(covariate_cols)) {
    bad_cov <- setdiff(covariate_cols, names(data))
    if (length(bad_cov) > 0) {
      warning(paste("Covariate columns not found in data:", paste(bad_cov, collapse = ", ")))
      covariate_cols <- intersect(covariate_cols, names(data))
    }
    if (length(covariate_cols) > 0) {
      cov_outliers <- .flag_covariate_outliers(data, id_col, covariate_cols, iqr_multiplier)
    }
  }

  # -- Resolve TAD for observations ---------------------------------------------
  obs_all <- data[data[[evid_col]] == 0 & data[[mdv_col]] == 0, ]
  tad_resolved <- .resolve_gross_tad(
    obs       = obs_all,
    all_data  = data,
    columns   = columns,
    tad_col   = tad_col
  )
  tad_vals  <- tad_resolved$tad
  tad_label <- tad_resolved$label

  # Auto bin width
  if (is.null(tad_bin_width)) {
    max_tad <- if (length(tad_vals) > 0) max(tad_vals, na.rm = TRUE) else NA_real_
    if (is.finite(max_tad) && max_tad > 0) {
      raw_w     <- max_tad / 20
      mag       <- 10^floor(log10(max(raw_w, 0.01)))
      nice      <- c(1, 2, 5, 10) * mag
      tad_bin_width <- nice[which.min(abs(nice - raw_w))]
    } else {
      tad_bin_width <- 1
    }
    tad_bin_width <- max(tad_bin_width, 0.25)
  }

  # -- Concentration TAD-bin outliers ------------------------------------------
  bin_result <- .flag_conc_bin_outliers(
    obs       = obs_all,
    tad_vals  = tad_vals,
    columns   = columns,
    iqr_mult  = iqr_multiplier,
    bin_width = tad_bin_width
  )
  conc_bin_outliers  <- bin_result$flagged
  conc_scatter_data  <- bin_result$scatter

  # -- Concentration elevations -------------------------------------------------
  conc_flags <- .flag_concentration_increases(
    data, columns, conc_increase_threshold, conc_floor
  )

  # -- Verbose summary ----------------------------------------------------------
  if (verbose) {
    message(sprintf(
      "\nirxclean :: Gross data quality assessment\n%s",
      paste(rep("-", 50), collapse = "")
    ))
    message(sprintf("  Subjects    : %d", n_subj))
    message(sprintf("  Observations: %d", n_obs))
    if (!is.null(cov_outliers) && nrow(cov_outliers) > 0) {
      message(sprintf("  Covariate outliers flagged      : %d subjects", length(unique(cov_outliers[[id_col]]))))
    } else {
      message("  Covariate outliers flagged      : none")
    }
    if (!is.null(conc_bin_outliers) && nrow(conc_bin_outliers) > 0) {
      message(sprintf(
        "  Concentration TAD-bin outliers  : %d observations", nrow(conc_bin_outliers)
      ))
    } else {
      message(sprintf(
        "  Concentration TAD-bin outliers  : none (%s bin width = %g)",
        tad_label, tad_bin_width
      ))
    }
    if (!is.null(conc_flags) && nrow(conc_flags) > 0) {
      message(sprintf("  Concentration elevation flags   : %d events", nrow(conc_flags)))
    } else {
      message("  Concentration elevation flags   : none")
    }
    n_miss_cols <- sum(miss_summary$n_missing > 0)
    message(sprintf("  Columns with missing values     : %d", n_miss_cols))
  }

  structure(
    list(
      covariate_outliers         = cov_outliers,
      concentration_bin_outliers = conc_bin_outliers,
      conc_scatter_data          = conc_scatter_data,
      concentration_flags        = conc_flags,
      missing_summary            = miss_summary,
      n_subjects                 = n_subj,
      n_observations             = n_obs,
      columns                    = columns,
      covariate_cols             = covariate_cols,
      categorical_cols           = categorical_cols,
      iqr_multiplier             = iqr_multiplier,
      tad_bin_width              = tad_bin_width,
      tad_label                  = tad_label
    ),
    class = "irxclean_quality"
  )
}


#' Print method for irxclean_quality objects
#' @param x An `irxclean_quality` object.
#' @param ... Ignored.
#' @export
print.irxclean_quality <- function(x, ...) {
  cat("irxclean gross data quality assessment\n")
  cat(sprintf("  Subjects                  : %d\n", x$n_subjects))
  cat(sprintf("  Observations              : %d\n", x$n_observations))
  n_cov <- if (!is.null(x$covariate_outliers)) nrow(x$covariate_outliers) else 0L
  cat(sprintf("  Covariate outlier rows    : %d\n", n_cov))
  n_bin <- if (!is.null(x$concentration_bin_outliers)) nrow(x$concentration_bin_outliers) else 0L
  cat(sprintf("  Concentration bin outliers: %d\n", n_bin))
  n_cf  <- if (!is.null(x$concentration_flags)) nrow(x$concentration_flags) else 0L
  cat(sprintf("  Concentration elev. flags : %d\n", n_cf))
  n_miss <- sum(x$missing_summary$n_missing > 0)
  cat(sprintf("  Columns with NAs          : %d\n", n_miss))
  invisible(x)
}


#' Plot method for irxclean_quality objects
#'
#' Generates a panel of up to three diagnostic plots:
#' covariate distributions with outlier bounds, flagged concentration profiles,
#' and missing-data rates.
#'
#' @param x An `irxclean_quality` object from [assess_data_quality()].
#' @param ... Ignored.
#'
#' @return A `ggplot` object (or a list of `ggplot` objects).
#'
#' @export
plot.irxclean_quality <- function(x, ...) {
  plots <- list()

  # Missing data bar chart
  miss <- x$missing_summary[x$missing_summary$n_missing > 0, ]
  if (nrow(miss) > 0) {
    plots$missing <- ggplot2::ggplot(
      miss,
      ggplot2::aes(
        x    = stats::reorder(.data$column, .data$pct_missing),
        y    = .data$pct_missing
      )
    ) +
      ggplot2::geom_col(fill = "#003B4F") +
      ggplot2::coord_flip() +
      ggplot2::labs(
        title = "Missing data by column",
        x     = NULL,
        y     = "% missing"
      ) +
      .clean_theme()
  }

  # DV vs TAD scatter: all obs (low alpha) + outliers (high alpha, red)
  if (!is.null(x$conc_scatter_data) && nrow(x$conc_scatter_data) > 0) {
    id_col  <- x$columns$ID
    scat    <- x$conc_scatter_data
    scat$LABEL <- ifelse(scat$OUTLIER, "Outlier", "Normal")
    scat$LABEL <- factor(scat$LABEL, levels = c("Normal", "Outlier"))

    # Separate for layering (outliers drawn on top)
    scat_normal  <- scat[!scat$OUTLIER, ]
    scat_outlier <- scat[scat$OUTLIER, ]

    p_scat <- ggplot2::ggplot(
      scat,
      ggplot2::aes(x = .data$TAD, y = .data$DV)
    ) +
      ggplot2::geom_point(
        data   = scat_normal,
        colour = "#003B4F",
        alpha  = 0.08,
        size   = 1.2
      ) +
      ggplot2::geom_point(
        data   = scat_outlier,
        colour = "#AD0000",
        alpha  = 0.9,
        size   = 2
      ) +
      ggplot2::geom_text(
        data   = scat_outlier,
        ggplot2::aes(label = .data[[id_col]]),
        colour = "#AD0000",
        size   = 2.5,
        vjust  = -0.7
      ) +
      ggplot2::labs(
        title    = "Concentration outliers by TAD bin",
        subtitle = sprintf(
          "IQR multiplier = %.1f  |  %s bin width = %g  |  red = flagged outlier",
          x$iqr_multiplier,
          x$tad_label,
          x$tad_bin_width
        ),
        x = x$tad_label,
        y = "DV"
      ) +
      .clean_theme()

    plots$concentration_bin_outliers <- p_scat
  }

  # Covariate outlier dot plots
  if (!is.null(x$covariate_outliers) && nrow(x$covariate_outliers) > 0) {
    id_col <- x$columns$ID
    cov_data <- x$covariate_outliers
    plots$covariate_outliers <- ggplot2::ggplot(
      cov_data,
      ggplot2::aes(
        x = .data$COVARIATE,
        y = .data$VALUE
      )
    ) +
      ggplot2::geom_point(colour = "#AD0000", size = 2) +
      ggplot2::geom_errorbar(
        ggplot2::aes(ymin = .data$LOWER_BOUND, ymax = .data$UPPER_BOUND),
        colour = "#5E5E5E", width = 0.2
      ) +
      ggplot2::labs(
        title    = "Covariate outliers",
        subtitle = sprintf("IQR multiplier = %.1f", x$covariate_outliers$IQR_MULT[1]),
        x        = "Covariate",
        y        = "Value"
      ) +
      ggplot2::facet_wrap(~COVARIATE, scales = "free") +
      .clean_theme()
  }

  # Concentration flags
  if (!is.null(x$concentration_flags) && nrow(x$concentration_flags) > 0) {
    cf <- x$concentration_flags
    plots$conc_flags <- ggplot2::ggplot(
      cf,
      ggplot2::aes(x = .data$RATIO)
    ) +
      ggplot2::geom_histogram(
        bins   = max(5L, min(20L, ceiling(log2(nrow(cf)) + 1L))),
        fill   = "#00769E",
        colour = "white"
      ) +
      ggplot2::labs(
        title    = "Concentration increases without intervening dose",
        subtitle = "Distribution of concentration fold-increases",
        x        = "Fold-increase (C_next / C_prev)",
        y        = "Count"
      ) +
      .clean_theme()
  }

  if (length(plots) == 0) {
    message("No issues detected; no plots generated.")
    return(invisible(NULL))
  }
  if (length(plots) == 1) return(plots[[1]])
  plots
}


# =============================================================================
# Internal helpers
# =============================================================================

#' @keywords internal
.summarise_missing <- function(data) {
  n_total <- nrow(data)
  missing_counts <- vapply(data, function(col) sum(is.na(col)), integer(1L))
  data.frame(
    column    = names(missing_counts),
    n_missing = as.integer(missing_counts),
    pct_missing = round(100 * missing_counts / n_total, 2),
    stringsAsFactors = FALSE
  )
}


#' @keywords internal
.flag_covariate_outliers <- function(data, id_col, covariate_cols, iqr_mult) {
  # Reduce to one row per subject (covariates should be time-invariant)
  subj_data <- data[!duplicated(data[[id_col]]), c(id_col, covariate_cols), drop = FALSE]

  flagged <- list()
  for (col in covariate_cols) {
    vals <- subj_data[[col]]
    vals_clean <- vals[!is.na(vals)]
    if (length(vals_clean) < 4L) next  # too few subjects to compute IQR

    q1  <- stats::quantile(vals_clean, 0.25)
    q3  <- stats::quantile(vals_clean, 0.75)
    iqr <- q3 - q1
    lo  <- q1 - iqr_mult * iqr
    hi  <- q3 + iqr_mult * iqr
    # Cap lower bound at 0 for non-negative quantities (age, weight, conc, etc.)
    if (all(vals_clean >= 0, na.rm = TRUE)) lo <- max(lo, 0)

    idx <- !is.na(vals) & (vals < lo | vals > hi)
    if (!any(idx)) next

    out_rows <- subj_data[idx, , drop = FALSE]
    flagged[[col]] <- data.frame(
      ID          = out_rows[[id_col]],
      COVARIATE   = col,
      VALUE       = out_rows[[col]],
      LOWER_BOUND = lo,
      UPPER_BOUND = hi,
      IQR_MULT    = iqr_mult,
      stringsAsFactors = FALSE
    )
    names(flagged[[col]])[1L] <- id_col
  }

  if (length(flagged) == 0L) return(NULL)
  do.call(rbind, flagged)
}


#' @keywords internal
.resolve_gross_tad <- function(obs, all_data, columns, tad_col) {
  id_col   <- columns$ID
  time_col <- columns$TIME
  evid_col <- columns$EVID

  if (nrow(obs) == 0L) {
    return(list(tad = numeric(0L), label = "Time after dose"))
  }

  # 1. Explicit tad_col
  if (!is.null(tad_col) && tad_col %in% names(obs)) {
    return(list(tad = obs[[tad_col]], label = "Time after dose"))
  }

  # 2. Auto-detect common TAD column names
  found <- intersect(c("TAD", "TAFD", "TSFD", "tad", "tafd"), names(obs))
  if (length(found) > 0) {
    return(list(tad = obs[[found[1L]]], label = "Time after dose"))
  }

  # 3. Compute from dose records (EVID == 1 or AMT > 0)
  if (evid_col %in% names(all_data)) {
    dose_flag <- all_data[[evid_col]] == 1
    if (!any(dose_flag, na.rm = TRUE) && "AMT" %in% names(all_data)) {
      dose_flag <- all_data[["AMT"]] > 0 & !is.na(all_data[["AMT"]])
    }
    if (any(dose_flag, na.rm = TRUE)) {
      message(
        "TAD column not found - computing from dose records (EVID==1 / AMT>0)."
      )
      ids      <- unique(obs[[id_col]])
      tad_list <- lapply(ids, function(id) {
        subj_all  <- all_data[all_data[[id_col]] == id, ]
        subj_all  <- subj_all[order(subj_all[[time_col]]), ]
        dose_times <- subj_all[[time_col]][dose_flag[all_data[[id_col]] == id]]
        obs_rows  <- obs[obs[[id_col]] == id, ]
        if (nrow(obs_rows) == 0L || length(dose_times) == 0L) {
          return(data.frame(
            row_idx = which(obs[[id_col]] == id),
            TAD     = NA_real_
          ))
        }
        tad_vals <- vapply(obs_rows[[time_col]], function(t) {
          prior <- dose_times[dose_times <= t]
          if (length(prior) == 0L) NA_real_ else t - max(prior)
        }, numeric(1L))
        data.frame(row_idx = which(obs[[id_col]] == id), TAD = tad_vals)
      })
      tad_df <- do.call(rbind, tad_list)
      if (is.null(tad_df) || nrow(tad_df) == 0L) {
        return(list(tad = rep(NA_real_, nrow(obs)), label = "Time after dose"))
      }
      tad_df <- tad_df[order(tad_df$row_idx), ]
      return(list(tad = tad_df$TAD, label = "Time after dose"))
    }
  }

  # 4. Fall back to TIME
  warning(
    "No TAD column found and no dose records identified. ",
    "Using TIME as proxy for time after dose."
  )
  list(tad = obs[[time_col]], label = "Time")
}


#' @keywords internal
.flag_conc_bin_outliers <- function(
  obs, tad_vals, columns, iqr_mult, bin_width
) {
  id_col   <- columns$ID
  time_col <- columns$TIME
  dv_col   <- columns$DV

  if (nrow(obs) == 0L) {
    scatter <- data.frame(
      ID = character(), TIME = numeric(), TAD = numeric(),
      DV = numeric(), OUTLIER = logical()
    )
    names(scatter)[1L] <- id_col
    return(list(flagged = NULL, scatter = scatter))
  }

  obs$TAD     <- tad_vals
  obs$TAD_BIN <- floor(obs$TAD / bin_width) * bin_width
  obs$OUTLIER <- FALSE

  flagged <- list()
  for (bin in sort(unique(obs$TAD_BIN))) {
    if (is.na(bin)) next
    bin_rows   <- obs[!is.na(obs$TAD_BIN) & obs$TAD_BIN == bin, ]
    vals       <- bin_rows[[dv_col]]
    vals_clean <- vals[!is.na(vals)]
    if (length(vals_clean) < 4L) next

    q1  <- stats::quantile(vals_clean, 0.25)
    q3  <- stats::quantile(vals_clean, 0.75)
    iqr <- q3 - q1
    lo  <- q1 - iqr_mult * iqr
    hi  <- q3 + iqr_mult * iqr

    idx_flag <- !is.na(vals) & (vals < lo | vals > hi)
    if (!any(idx_flag)) next

    # Mark in obs
    bin_row_idx <- which(!is.na(obs$TAD_BIN) & obs$TAD_BIN == bin)
    obs$OUTLIER[bin_row_idx[idx_flag]] <- TRUE

    out_rows <- bin_rows[idx_flag, , drop = FALSE]
    flagged[[length(flagged) + 1L]] <- data.frame(
      ID          = out_rows[[id_col]],
      TIME        = out_rows[[time_col]],
      TAD         = out_rows$TAD,
      TAD_BIN     = bin,
      DV          = out_rows[[dv_col]],
      LOWER_BOUND = lo,
      UPPER_BOUND = hi,
      N_IN_BIN    = length(vals_clean),
      stringsAsFactors = FALSE
    )
  }

  flagged_df <- if (length(flagged) > 0) {
    out <- do.call(rbind, flagged)
    names(out)[1L] <- id_col
    out
  } else {
    NULL
  }

  # Scatter data: all observations with TAD and OUTLIER flag
  scatter <- data.frame(
    ID      = obs[[id_col]],
    TIME    = obs[[time_col]],
    TAD     = obs$TAD,
    DV      = obs[[dv_col]],
    OUTLIER = obs$OUTLIER,
    stringsAsFactors = FALSE
  )
  names(scatter)[1L] <- id_col

  list(flagged = flagged_df, scatter = scatter)
}


#' @keywords internal
.flag_concentration_increases <- function(data, columns, threshold,
                                          conc_floor = 0) {
  id_col   <- columns$ID
  time_col <- columns$TIME
  dv_col   <- columns$DV
  evid_col <- columns$EVID
  mdv_col  <- columns$MDV

  flagged <- list()

  for (id in unique(data[[id_col]])) {
    subj <- data[data[[id_col]] == id, ]
    subj <- subj[order(subj[[time_col]]), ]

    # Observation rows only
    obs <- subj[subj[[evid_col]] == 0 & subj[[mdv_col]] == 0, ]
    if (nrow(obs) < 2L) next

    for (j in seq(2L, nrow(obs))) {
      t1 <- obs[[time_col]][j - 1L]
      t2 <- obs[[time_col]][j]
      c1 <- obs[[dv_col]][j - 1L]
      c2 <- obs[[dv_col]][j]

      # Skip if either value is missing, non-positive, or at/below the LOQ
      # floor (ratio is noise-dominated when c1 is near zero)
      if (is.na(c1) || is.na(c2) || c1 <= conc_floor) next

      # Any dose between t1 and t2?
      between_dose <- any(
        subj[[evid_col]] == 1 &
          subj[[time_col]] >  t1 &
          subj[[time_col]] <= t2
      )
      if (!between_dose && (c2 / c1) > threshold) {
        flagged[[length(flagged) + 1L]] <- data.frame(
          ID       = id,
          TIME_PREV = t1,
          DV_PREV  = c1,
          TIME     = t2,
          DV       = c2,
          RATIO    = c2 / c1,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (length(flagged) == 0L) return(NULL)
  out <- do.call(rbind, flagged)
  names(out)[1L] <- id_col
  out
}
