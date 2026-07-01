# =============================================================================
# Pre-modelling exclusion criteria
# =============================================================================

#' Apply structured exclusion criteria to a NONMEM-format dataset
#'
#' Applies a battery of model-independent exclusion checks and returns a
#' flagged dataset together with a per-criterion summary.  This function is
#' designed for the \emph{initial} data-cleaning step run once after data
#' wrangling, before a stable base model is available.
#'
#' Each criterion adds a logical \code{excl_*} column to the data.  Rows where
#' \emph{any} exclusion column is \code{TRUE} are removed to produce
#' \code{data_clean}.
#'
#' @param data A data frame in NONMEM format.
#' @param columns Named list mapping logical column roles to actual column
#'   names.  Required keys: \code{ID}, \code{TIME}, \code{DV}, \code{EVID},
#'   \code{MDV}.  Optional: \code{AMT}, \code{RATE} (needed for dose-outlier
#'   and infusion checks).
#' @param check_no_doses Logical. Flag all rows for subjects who have no dose
#'   records (EVID == 1)?  Default \code{TRUE}.
#' @param check_no_tdm Logical. Flag all rows for subjects who have no TDM
#'   observations (EVID == 0, MDV == 0)?  Default \code{TRUE}.
#' @param check_during_infusion Logical. Flag observation rows collected
#'   during an ongoing infusion?  Requires \code{AMT} and \code{RATE} columns.
#'   Default \code{TRUE}.
#' @param check_conc_increase Logical. Flag observations where concentration
#'   rises by more than \code{conc_increase_threshold}-fold without an
#'   intervening dose?  Default \code{TRUE}.
#' @param conc_increase_threshold Numeric. Fold-increase threshold for
#'   flagging concentration elevations.  Default \code{2.0}.
#' @param conc_floor Numeric. Minimum preceding concentration for evaluating
#'   fold-increases (pairs where DV_prev <= \code{conc_floor} are skipped).
#'   Default \code{0}.
#' @param check_dose_outlier Logical. Flag dose records whose AMT is an
#'   IQR-based outlier across all dose records?  Requires \code{AMT} column.
#'   Default \code{TRUE}.
#' @param dose_iqr_multiplier Numeric. IQR multiplier for dose-outlier
#'   detection.  Default \code{5}.
#' @param custom_criteria Named list of row-level logical vectors (one per
#'   criterion, length equal to \code{nrow(data)}).  Each vector is added as
#'   an \code{excl_<name>} column.  \code{NULL} skips custom criteria.
#' @param verbose Logical. Print a summary to the console?  Default \code{TRUE}.
#'
#' @return An object of class \code{irxclean_exclusions} (a named list):
#' \describe{
#'   \item{data_flagged}{The original data frame with \code{excl_*} logical
#'     columns appended.}
#'   \item{data_clean}{Rows from \code{data} where no \code{excl_*} column
#'     is \code{TRUE}.}
#'   \item{exclusion_summary}{Data frame with columns \code{criterion},
#'     \code{n_records}, \code{n_subjects}, and \code{pct_records} (percentage
#'     of total observation rows affected).}
#'   \item{criteria_applied}{Character vector of \code{excl_*} column names
#'     added.}
#' }
#'
#' @seealso [assess_data_quality()], [remove_erroneous_obs()]
#'
#' @examples
#' \dontrun{
#' dat <- read.csv("my_data.csv")
#' excl <- apply_exclusion_criteria(
#'   dat,
#'   check_dose_outlier = TRUE,
#'   dose_iqr_multiplier = 3
#' )
#' print(excl)
#' head(excl$data_clean)
#' }
#'
#' @export
apply_exclusion_criteria <- function(
  data,
  columns = list(
    ID   = "ID",
    TIME = "TIME",
    DV   = "DV",
    EVID = "EVID",
    MDV  = "MDV",
    AMT  = "AMT",
    RATE = "RATE"
  ),
  check_no_doses          = TRUE,
  check_no_tdm            = TRUE,
  check_during_infusion   = TRUE,
  check_conc_increase     = TRUE,
  conc_increase_threshold = 2.0,
  conc_floor              = 0,
  check_dose_outlier      = TRUE,
  dose_iqr_multiplier     = 5,
  custom_criteria         = NULL,
  verbose                 = TRUE
) {
  # -- Validate required columns -----------------------------------------------
  required <- c("ID", "TIME", "DV", "EVID", "MDV")
  missing_keys <- setdiff(required, names(columns))
  if (length(missing_keys) > 0) {
    stop(paste("columns list missing keys:", paste(missing_keys, collapse = ", ")))
  }
  for (key in required) {
    if (!columns[[key]] %in% names(data)) {
      stop(sprintf("Column '%s' ('%s') not found in data.", key, columns[[key]]))
    }
  }

  id_col   <- columns$ID
  time_col <- columns$TIME
  evid_col <- columns$EVID
  mdv_col  <- columns$MDV
  amt_col  <- if (!is.null(columns$AMT) && columns$AMT %in% names(data))
    columns$AMT else NULL
  rate_col <- if (!is.null(columns$RATE) && columns$RATE %in% names(data))
    columns$RATE else NULL

  flagged          <- data
  criteria_applied <- character(0L)

  is_obs     <- data[[evid_col]] == 0 & data[[mdv_col]] == 0
  n_total_obs <- sum(is_obs, na.rm = TRUE)
  all_ids    <- unique(data[[id_col]])

  # -- 1. Subjects with no dose records ----------------------------------------
  if (check_no_doses) {
    subj_has_dose <- vapply(all_ids, function(id) {
      any(data[[evid_col]][data[[id_col]] == id] == 1, na.rm = TRUE)
    }, logical(1L))
    no_dose_ids <- all_ids[!subj_has_dose]
    flagged$excl_no_doses <- data[[id_col]] %in% no_dose_ids
    criteria_applied <- c(criteria_applied, "excl_no_doses")
  }

  # -- 2. Subjects with no TDM observations ------------------------------------
  if (check_no_tdm) {
    subj_has_obs <- vapply(all_ids, function(id) {
      mask <- data[[id_col]] == id
      any(data[[evid_col]][mask] == 0 & data[[mdv_col]][mask] == 0, na.rm = TRUE)
    }, logical(1L))
    no_tdm_ids <- all_ids[!subj_has_obs]
    flagged$excl_no_tdm <- data[[id_col]] %in% no_tdm_ids
    criteria_applied <- c(criteria_applied, "excl_no_tdm")
  }

  # -- 3. Observations during ongoing infusion ---------------------------------
  if (check_during_infusion) {
    flagged$excl_during_infusion <- FALSE
    if (!is.null(amt_col) && !is.null(rate_col)) {
      for (id in all_ids) {
        subj_idx  <- which(data[[id_col]] == id)
        ord       <- order(data[[time_col]][subj_idx])
        subj_idx  <- subj_idx[ord]
        subj      <- data[subj_idx, ]

        for (j in seq_len(nrow(subj))) {
          if (subj[[evid_col]][j] != 0L || subj[[mdv_col]][j] != 0L) next
          t_obs <- subj[[time_col]][j]

          # Most recent prior infusion dose (RATE > 0, at or before obs)
          prior <- which(
            subj[[evid_col]] == 1L &
              !is.na(subj[[rate_col]]) & subj[[rate_col]] > 0 &
              subj[[time_col]] <= t_obs
          )
          if (length(prior) == 0L) next

          last_d <- prior[length(prior)]
          t_dose <- subj[[time_col]][last_d]
          amt    <- subj[[amt_col]][last_d]
          rate   <- subj[[rate_col]][last_d]

          if (!is.na(amt) && !is.na(rate) && rate > 0 &&
              t_obs < t_dose + (amt / rate)) {
            flagged$excl_during_infusion[subj_idx[j]] <- TRUE
          }
        }
      }
    } else {
      if (verbose) {
        message(
          "  check_during_infusion skipped: AMT and/or RATE column not found."
        )
      }
    }
    criteria_applied <- c(criteria_applied, "excl_during_infusion")
  }

  # -- 4. Concentration increases without intervening dose ---------------------
  if (check_conc_increase) {
    flagged$excl_conc_increase <- FALSE
    cf <- .flag_concentration_increases(
      data, columns, conc_increase_threshold, conc_floor
    )
    if (!is.null(cf) && nrow(cf) > 0) {
      for (k in seq_len(nrow(cf))) {
        idx <- which(
          data[[id_col]] == cf[[id_col]][k] &
            data[[time_col]] == cf$TIME[k] &
            data[[evid_col]] == 0L & data[[mdv_col]] == 0L
        )
        if (length(idx) > 0L) flagged$excl_conc_increase[idx] <- TRUE
      }
    }
    criteria_applied <- c(criteria_applied, "excl_conc_increase")
  }

  # -- 5. Dose AMT IQR outlier -------------------------------------------------
  if (check_dose_outlier) {
    flagged$excl_dose_outlier <- FALSE
    if (!is.null(amt_col)) {
      dose_rows <- which(data[[evid_col]] == 1L)
      if (length(dose_rows) >= 4L) {
        amts       <- data[[amt_col]][dose_rows]
        amts_clean <- amts[!is.na(amts)]
        if (length(amts_clean) >= 4L) {
          q1  <- stats::quantile(amts_clean, 0.25)
          q3  <- stats::quantile(amts_clean, 0.75)
          iqr <- q3 - q1
          lo  <- max(q1 - dose_iqr_multiplier * iqr, 0)
          hi  <- q3 + dose_iqr_multiplier * iqr
          out_idx <- dose_rows[!is.na(amts) & (amts < lo | amts > hi)]
          flagged$excl_dose_outlier[out_idx] <- TRUE
        }
      }
    } else {
      if (verbose) {
        message("  check_dose_outlier skipped: AMT column not found.")
      }
    }
    criteria_applied <- c(criteria_applied, "excl_dose_outlier")
  }

  # -- 6. Custom criteria -------------------------------------------------------
  if (!is.null(custom_criteria)) {
    if (!is.list(custom_criteria) || is.null(names(custom_criteria))) {
      stop("`custom_criteria` must be a named list.")
    }
    for (crit_name in names(custom_criteria)) {
      col_name <- paste0("excl_", crit_name)
      vec <- custom_criteria[[crit_name]]
      if (length(vec) != nrow(data)) {
        stop(sprintf(
          "Custom criterion '%s' has length %d; expected %d (nrow(data)).",
          crit_name, length(vec), nrow(data)
        ))
      }
      flagged[[col_name]] <- as.logical(vec)
      criteria_applied <- c(criteria_applied, col_name)
    }
  }

  # -- Build exclusion summary --------------------------------------------------
  excl_cols <- criteria_applied

  summary_rows <- lapply(excl_cols, function(col) {
    flag_vec <- flagged[[col]]
    n_rec    <- sum(flag_vec, na.rm = TRUE)
    n_subj   <- length(unique(data[[id_col]][!is.na(flag_vec) & flag_vec]))
    pct_rec  <- if (n_total_obs > 0) {
      round(100 * sum(flag_vec & is_obs, na.rm = TRUE) / n_total_obs, 2)
    } else {
      NA_real_
    }
    data.frame(
      criterion   = col,
      n_records   = n_rec,
      n_subjects  = n_subj,
      pct_records = pct_rec,
      stringsAsFactors = FALSE
    )
  })
  excl_summary <- if (length(summary_rows) > 0) {
    do.call(rbind, summary_rows)
  } else {
    data.frame(
      criterion = character(), n_records = integer(),
      n_subjects = integer(), pct_records = numeric(),
      stringsAsFactors = FALSE
    )
  }

  # -- Clean data: rows where ANY excl_* is TRUE --------------------------------
  any_excl <- if (length(excl_cols) > 0) {
    Reduce(`|`, lapply(excl_cols, function(col) {
      v <- flagged[[col]]
      !is.na(v) & v
    }))
  } else {
    rep(FALSE, nrow(data))
  }
  data_clean <- data[!any_excl, , drop = FALSE]

  # -- Verbose output -----------------------------------------------------------
  if (verbose) {
    message(sprintf(
      "\nirxclean :: Exclusion criteria applied\n%s",
      paste(rep("-", 50), collapse = "")
    ))
    message(sprintf("  Total rows      : %d", nrow(data)))
    message(sprintf("  Subjects        : %d", length(all_ids)))
    message(sprintf("  Observation rows: %d", n_total_obs))
    message(sprintf("  Rows excluded   : %d", sum(any_excl)))
    message(sprintf("  Rows retained   : %d", nrow(data_clean)))
    if (nrow(excl_summary) > 0) {
      for (i in seq_len(nrow(excl_summary))) {
        message(sprintf(
          "    %-30s: %d records, %d subjects (%.1f%% of obs)",
          excl_summary$criterion[i],
          excl_summary$n_records[i],
          excl_summary$n_subjects[i],
          excl_summary$pct_records[i]
        ))
      }
    }
  }

  structure(
    list(
      data_flagged      = flagged,
      data_clean        = data_clean,
      exclusion_summary = excl_summary,
      criteria_applied  = criteria_applied
    ),
    class = "irxclean_exclusions"
  )
}


#' Print method for irxclean_exclusions objects
#' @param x An \code{irxclean_exclusions} object.
#' @param ... Ignored.
#' @export
print.irxclean_exclusions <- function(x, ...) {
  cat("irxclean exclusion criteria summary\n")
  cat(sprintf("  Criteria applied : %d\n", length(x$criteria_applied)))
  cat(sprintf("  Total rows       : %d\n", nrow(x$data_flagged)))
  cat(sprintf("  Rows retained    : %d\n", nrow(x$data_clean)))
  cat(sprintf("  Rows excluded    : %d\n",
              nrow(x$data_flagged) - nrow(x$data_clean)))
  if (nrow(x$exclusion_summary) > 0) {
    cat("\nExclusion summary:\n")
    print(x$exclusion_summary, row.names = FALSE)
  }
  invisible(x)
}
