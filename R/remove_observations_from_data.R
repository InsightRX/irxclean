#' Apply user-reviewed observation exclusions to produce a cleaned dataset
#'
#' After running \code{remove_erroneous_obs()} and reviewing the flagged
#' observations, use this function to apply exclusions and produce a cleaned
#' dataset suitable for final analysis.  The top `n` flagged observations
#' (by iteration order) are excluded.
#'
#' @section User review:
#' Flagged observations should be investigated before calling this function.
#' A flagged observation is not necessarily erroneous; it may represent a
#' genuine patient population subgroup or an unusual but valid measurement.
#' Final exclusion decisions rest with the analyst.
#'
#' Matching is performed on rounded (`3` decimal places) `ID`, `TIME`, and
#' `DV` values.  A warning is raised if the number of rows removed does not
#' equal `n`.
#'
#' @param data A data frame containing the original NONMEM dataset.
#' @param results A results list from [remove_erroneous_obs()] or
#'   [load_removal_results()].
#' @param n Integer. Number of observations to remove (taken from the first
#'   `n` iterations in `results$rem`).  Default: all removed observations.
#' @param columns Named list specifying column labels for matching.
#'   Default: `list(ID = "ID", TIME = "TIME", DV = "DV")`.
#'
#' @return The input `data` frame with the identified observations removed.
#'
#' @examples
#' \dontrun{
#' original_data <- read.csv("my_data.csv")
#' results       <- load_removal_results("err_rem1")
#' cleaned       <- remove_observations_from_data(original_data, results, n = 5)
#' }
#'
#' @export
remove_observations_from_data <- function(
    data,
    results,
    n       = NULL,
    columns = list(ID = "ID", TIME = "TIME", DV = "DV")
) {
  .validate_results(results)

  rm <- results$rem
  if (!is.null(n)) {
    n <- as.integer(n)
    if (is.na(n) || n < 1L) stop("`n` must be a positive integer.")
    rm <- dplyr::slice(rm, seq_len(min(n, nrow(rm))))
  }

  id_col   <- columns$ID
  time_col <- columns$TIME
  dv_col   <- columns$DV

  new_data <- data |>
    dplyr::filter(
      !(
        round(.data[[id_col]],   3) %in% round(rm[[id_col]],   3) &
        round(.data[[time_col]], 3) %in% round(rm[[time_col]], 3) &
        round(.data[[dv_col]],   3) %in% round(rm[[dv_col]],   3)
      )
    )

  n_requested <- nrow(rm)
  n_removed   <- nrow(data) - nrow(new_data)

  if (n_removed != n_requested) {
    warning(sprintf(
      "Removed %d rows but requested %d.  Check that column names match.",
      n_removed, n_requested
    ))
  }

  new_data
}
