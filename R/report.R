#' Generate an automated HTML cleaning report
#'
#' Renders an R Markdown report summarising the iterative outlier-removal
#' workflow.  The report is structured into two sections:
#' \enumerate{
#'   \item \strong{Key Results} — always visible: observations removed,
#'     pOFV change, error-parameter stability (if detected), and nRMSE trend.
#'   \item \strong{Detailed Results} — collapsible: full THETA, OMEGA, and
#'     SIGMA stability facet plots plus the model stability assessment.
#' }
#' An optional demographic comparison section is rendered when \code{data}
#' and covariate columns are supplied.
#'
#' @param results A results list from [remove_erroneous_obs()] or
#'   [load_removal_results()].  **Required.**
#' @param data A data frame containing the original NONMEM dataset.  Required
#'   for demographic comparison; `NULL` skips that section.
#' @param prop_error_param Character or \code{NULL}.  Column name of the
#'   proportional residual error parameter to highlight in the Key Results
#'   section.  Accepts the column name as stored in \code{results$par} (e.g.
#'   \code{"THETA7"} or \code{"SIGMA.1.1"}), or NONMEM index notation such as
#'   \code{"THETA(7)"} or \code{"SIGMA(1)"} which is normalised automatically.
#'   \code{NULL} triggers auto-detection from \code{results$param_labels};
#'   the detected column is printed when \code{quiet = FALSE}.
#' @param add_error_param Character or \code{NULL}.  Column name of the
#'   additive residual error parameter.  Accepts the same formats as
#'   \code{prop_error_param}.  \code{NULL} triggers auto-detection.
#' @param continuous_cov_cols Character vector of continuous covariate column
#'   names for demographic comparison.  `NULL` skips continuous covariates.
#' @param categorical_cov_cols Character vector of categorical covariate column
#'   names.  `NULL` skips categorical covariates.
#' @param id_col Character. Subject ID column name.  Default `"ID"`.
#' @param time_col Character. Time column for TAD computation.  Default
#'   \code{"TIME"}.
#' @param evid_col Character. EVID column for identifying dose records.
#'   Default \code{"EVID"}.
#' @param tad_col Character or \code{NULL}.  Pre-computed TAD column name;
#'   \code{NULL} triggers auto-detection / computation.  Default \code{NULL}.
#' @param tad_bin_width Numeric or \code{NULL}.  TAD bin width for the
#'   curve-position plot.  \code{NULL} auto-selects.  Default \code{NULL}.
#' @param stability_threshold_pct Numeric. Parameter drift threshold (\%) passed
#'   to [check_model_stability()].  Default `20`.
#' @param title Character. Report title.
#'   Default `"irxclean Data Cleaning Report"`.
#' @param output_file Character. Output file name.
#'   Default `"irxclean_report.html"`.
#' @param output_dir Character. Directory for the output file.
#'   Default: current working directory.
#' @param open_report Logical. Open the report in the default browser after
#'   rendering?  Default `FALSE`.
#' @param quiet Logical. Suppress rmarkdown rendering messages?
#'   Default `TRUE`.
#'
#' @return Invisibly returns the path to the rendered HTML file.
#'
#' @section Requirements:
#' Requires the `rmarkdown` and `knitr` packages (listed in `Suggests`).
#'
#' @examples
#' \dontrun{
#' results <- load_removal_results("err_rem1")
#' data    <- read.csv("my_data.csv")
#'
#' generate_report(
#'   results              = results,
#'   data                 = data,
#'   continuous_cov_cols  = c("AGE", "WT", "HT"),
#'   categorical_cov_cols = "SEX",
#'   title                = "Busulfan PK Cleaning Report"
#' )
#' }
#'
#' @export
generate_report <- function(
    results,
    data                    = NULL,
    prop_error_param        = NULL,
    add_error_param         = NULL,
    continuous_cov_cols     = NULL,
    categorical_cov_cols    = NULL,
    id_col                  = "ID",
    time_col                = "TIME",
    evid_col                = "EVID",
    tad_col                 = NULL,
    tad_bin_width           = NULL,
    stability_threshold_pct = 20,
    title                   = "irxclean Data Cleaning Report",
    output_file             = "irxclean_report.html",
    output_dir              = ".",
    open_report             = FALSE,
    quiet                   = TRUE
) {
  if (!requireNamespace("rmarkdown", quietly = TRUE)) {
    stop("Package 'rmarkdown' is required. Install with: install.packages('rmarkdown')")
  }

  .validate_results(results)

  # Normalize override formats: "THETA(7)" -> "THETA7", "SIGMA(2)" -> "SIGMA.2.2"
  prop_error_param <- .parse_param_name(prop_error_param)
  add_error_param  <- .parse_param_name(add_error_param)

  # Auto-detect error parameters for any that are still NULL
  if (is.null(prop_error_param) || is.null(add_error_param)) {
    detected <- .find_error_params(
      results$param_labels, results$par,
      verbose = !quiet
    )
    if (is.null(prop_error_param)) prop_error_param <- detected$prop_col
    if (is.null(add_error_param))  add_error_param  <- detected$add_col
  }

  template <- system.file(
    "rmarkdown", "irxclean_report.Rmd",
    package = "irxclean"
  )
  if (!nzchar(template)) {
    stop("Report template not found. Re-install the irxclean package.")
  }

  # Pass data to the template via a temporary environment stored in global
  # options to avoid file serialisation issues
  report_env <- new.env(parent = emptyenv())
  report_env$results                 <- results
  report_env$data                    <- data
  report_env$prop_error_param        <- prop_error_param
  report_env$add_error_param         <- add_error_param
  report_env$continuous_cov_cols     <- continuous_cov_cols
  report_env$categorical_cov_cols    <- categorical_cov_cols
  report_env$id_col                  <- id_col
  report_env$time_col                <- time_col
  report_env$evid_col                <- evid_col
  report_env$tad_col                 <- tad_col
  report_env$tad_bin_width           <- tad_bin_width
  report_env$stability_threshold_pct <- stability_threshold_pct
  report_env$report_title            <- title

  old_opt <- getOption("irxclean.report_env")
  options(irxclean.report_env = report_env)
  on.exit(options(irxclean.report_env = old_opt), add = TRUE)

  output_path <- file.path(normalizePath(output_dir), output_file)

  rmarkdown::render(
    input       = template,
    output_file = output_path,
    quiet       = quiet,
    envir       = new.env(parent = globalenv())
  )

  message(sprintf("Report written to: %s", output_path))

  if (open_report) {
    utils::browseURL(output_path)
  }

  invisible(output_path)
}
