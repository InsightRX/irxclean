#' Iteratively detect potentially erroneous pharmacokinetic observations
#'
#' Performs iterative NONMEM-based detection of potentially erroneous
#' observations in a pharmacokinetic dataset.  The ranked list of flagged
#' observations is intended for user review before any exclusions are applied
#' to the final analysis dataset (see \code{remove_observations_from_data()}).
#' At each iteration the algorithm:
#' \enumerate{
#'   \item Runs NONMEM (via PsN `execute`) on the current dataset.
#'   \item Identifies the observation with the highest absolute CWRES.
#'   \item Removes that observation from the dataset.
#'   \item Saves parameter estimates and diagnostics for later plotting.
#' }
#'
#' @param dat File path to the NONMEM-ready `.csv` dataset (with or without
#'   the `.csv` extension, relative or absolute), **or** an R data frame that
#'   is already loaded in memory.
#' @param mod File path to the NONMEM `.mod` control stream (with or without
#'   the `.mod` extension, relative or absolute), **or** a \code{nm_model}
#'   list returned by \code{nm_read_model()}.
#' @param run_id Character. Unique run identifier (e.g. `"err_rem1"`).
#' @param n Integer. Number of observations to remove.
#' @param columns Named list specifying column labels in the dataset.
#'   Default: `list(ID = "ID", TIME = "TIME", DV = "DV")`.
#' @param verbose Logical. Print progress messages?  Default `TRUE`.
#' @param save_results Logical. Save an `.RDS` results file to disk?
#'   Default `TRUE`.
#' @param save_temp_dir Logical. Copy the temporary working directory to
#'   the current working directory (`TRUE`) or delete it after completion
#'   (`FALSE`, default)?
#' @param stability_check Logical. Append model stability metrics to the
#'   returned results?  Default `TRUE`.  See [check_model_stability()].
#'
#' @return A list containing:
#' \describe{
#'   \item{par}{Data frame of population parameter estimates per iteration.}
#'   \item{rem}{Data frame of removed observations with their iteration.}
#'   \item{phi}{Data frame of individual objective function (OBJ) values
#'     per iteration.}
#'   \item{rmse}{Data frame of RMSE per iteration.}
#'   \item{stability}{(If \code{stability_check = TRUE}) Stability assessment
#'     from \code{check_model_stability()}.}
#'   \item{param_labels}{Named list with elements \code{thetas} and
#'     \code{omegas}: human-readable parameter labels parsed from the
#'     \code{.mod} file comments, used automatically in plots.}
#' }
#'
#' @section Specifying inputs -- three supported modes:
#'
#' \strong{Mode 1 -- bare file name (files in the working directory)}
#'
#' Pass the stem of the filename without the extension.  The files must
#' exist in the current working directory.  This was the original interface
#' and remains the simplest option when your working directory is already set
#' to the folder that contains the data and model.
#'
#' \preformatted{
#' results <- remove_erroneous_obs(
#'   dat    = "my_data",    # reads my_data.csv from getwd()
#'   mod    = "my_model",   # reads my_model.mod from getwd()
#'   run_id = "clean_run1",
#'   n      = 10
#' )
#' }
#'
#' \strong{Mode 2 -- file path (relative or absolute, extension optional)}
#'
#' Pass a relative or absolute path.  The \code{.csv} / \code{.mod}
#' extension is stripped and re-appended automatically, so you can include
#' or omit it.  Useful when files live in a different directory from your
#' R session.
#'
#' \preformatted{
#' results <- remove_erroneous_obs(
#'   dat    = "/projects/busulfan/data/bu_data.csv",
#'   mod    = "/projects/busulfan/models/bu_base",
#'   run_id = "clean_run1",
#'   n      = 10
#' )
#' }
#'
#' \strong{Mode 3 -- R objects already loaded in memory}
#'
#' Pass a \code{data.frame} for \code{dat} and/or a \code{nm_model} list
#' (from \code{nm_read_model()}) for \code{mod}.  This is useful when you
#' have already read and pre-processed your data in R (e.g. after applying
#' \code{apply_exclusion_criteria()}), or when you want to modify the
#' model programmatically before running.  NONMEM still runs from temporary
#' files on disk -- the objects are written there automatically.
#'
#' \preformatted{
#' pk_clean <- apply_exclusion_criteria(...)$data_clean
#' nm_mod   <- nm_read_model("/projects/busulfan/models/bu_base.mod")
#'
#' results <- remove_erroneous_obs(
#'   dat    = pk_clean,   # data.frame -- no file needed
#'   mod    = nm_mod,     # nm_model list -- no file copy needed
#'   run_id = "clean_run1",
#'   n      = 10
#' )
#' }
#'
#' The three modes can be mixed: you may supply an in-memory data frame for
#' \code{dat} while pointing \code{mod} at a file path, or vice versa.
#'
#' @section PsN requirement:
#' `execute` and `sumo` must be available on the system `PATH`.
#'
#' @examples
#' \dontrun{
#' # Mode 1: bare file names (files in working directory)
#' results <- remove_erroneous_obs(
#'   dat    = "my_data",
#'   mod    = "my_model",
#'   run_id = "err_rem1",
#'   n      = 5
#' )
#'
#' # Mode 2: full file paths (extension optional)
#' results <- remove_erroneous_obs(
#'   dat    = "/projects/pk/my_data.csv",
#'   mod    = "/projects/pk/my_model.mod",
#'   run_id = "err_rem1",
#'   n      = 5
#' )
#'
#' # Mode 3: R objects already in memory
#' pk_data  <- read.csv("/projects/pk/my_data.csv")
#' nm_model <- nm_read_model("/projects/pk/my_model.mod")
#' results  <- remove_erroneous_obs(
#'   dat    = pk_data,
#'   mod    = nm_model,
#'   run_id = "err_rem1",
#'   n      = 5
#' )
#'
#' # Mix and match: in-memory data + file path for model
#' pk_clean <- apply_exclusion_criteria(pk_data)$data_clean
#' results  <- remove_erroneous_obs(
#'   dat    = pk_clean,
#'   mod    = "/projects/pk/my_model",
#'   run_id = "err_rem1",
#'   n      = 5
#' )
#'
#' plot_removal_metrics(results, metric = "pOFV")
#' }
#'
#' @seealso [plot_removal_metrics()], [check_model_stability()],
#'   [remove_observations_from_data()], [generate_report()]
#'
#' @export
remove_erroneous_obs <- function(
    dat,
    mod,
    run_id,
    n,
    columns        = list(ID = "ID", TIME = "TIME", DV = "DV"),
    verbose        = TRUE,
    save_results   = TRUE,
    save_temp_dir  = FALSE,
    stability_check = TRUE
) {
  # -- Input validation --------------------------------------------------------
  if (!is.character(run_id)) stop("`run_id` must be a character string.")
  n <- as.integer(n)
  if (is.na(n) || n < 1L) stop("`n` must be a positive integer.")

  # Resolve dat: file path or data frame
  if (is.character(dat)) {
    dat_stem <- tools::file_path_sans_ext(dat)
    dat_file <- normalizePath(paste0(dat_stem, ".csv"), mustWork = FALSE)
    if (!file.exists(dat_file)) stop(sprintf("Dataset not found: %s", dat_file))
    dat_obj  <- NULL
  } else if (is.data.frame(dat)) {
    dat_obj  <- dat
    dat_file <- NULL
  } else {
    stop("`dat` must be a file path (character) or a data frame.")
  }

  # Resolve mod: file path or nm_model list
  if (is.character(mod)) {
    mod_stem <- tools::file_path_sans_ext(mod)
    mod_file <- normalizePath(paste0(mod_stem, ".mod"), mustWork = FALSE)
    if (!file.exists(mod_file)) stop(sprintf("Model file not found: %s", mod_file))
    mod_obj  <- NULL
  } else if (is.list(mod)) {
    mod_obj  <- mod
    mod_file <- NULL
  } else {
    stop("`mod` must be a file path (character) or an nm_model list from nm_read_model().")
  }

  original_wd <- getwd()

  # -- Temporary working directory ---------------------------------------------
  temp_dir <- tempfile(pattern = paste0("irxclean_", run_id, "_"))
  dir.create(temp_dir, recursive = TRUE)

  # Always restore wd and optionally clean up temp dir
  on.exit({
    setwd(original_wd)
    if (!save_temp_dir) {
      unlink(temp_dir, recursive = TRUE, force = TRUE)
      if (verbose) message("Temporary directory deleted.")
    } else {
      dest <- file.path(original_wd, basename(temp_dir))
      file.copy(from = temp_dir, to = original_wd, recursive = TRUE, overwrite = TRUE)
      if (verbose) message(sprintf("Temporary directory saved to: %s", dest))
    }
  }, add = TRUE)

  if (verbose) message(sprintf("Temporary directory: %s", temp_dir))
  setwd(temp_dir)

  # -- Prepare dataset and model ------------------------------------------------
  # Data: copy from disk or use in-memory data frame directly
  if (!is.null(dat_file)) {
    file.copy(from = dat_file, to = basename(dat_file))
    data <- utils::read.csv(basename(dat_file))
  } else {
    data <- dat_obj
  }

  # Model: copy from disk + parse, or use supplied nm_model list directly
  if (!is.null(mod_file)) {
    mod_base  <- basename(mod_file)
    file.copy(from = mod_file, to = mod_base)
    nm_model  <- nm_read_model(mod_base)
    labels_src <- mod_base
  } else {
    mod_base   <- NULL
    nm_model   <- mod_obj
    labels_src <- NULL  # will fall back to new_mod_file after it is written
  }

  total_dv  <- sum(data$EVID == 0, na.rm = TRUE)

  removed_observations <- data.frame(
    ID        = integer(),
    TIME      = numeric(),
    DV        = numeric(),
    CWRES     = numeric(),
    ITERATION = integer(),
    stringsAsFactors = FALSE
  )

  par_list  <- vector("list", n)
  phi_list  <- vector("list", n)
  rem_list  <- vector("list", n)
  rmse_list <- vector("list", n)

  # Modify model: update $DATA and $TABLE blocks
  if (verbose && length(nm_model$DATA) > 1) {
    message(sprintf("Found %d $DATA blocks - using the first.", length(nm_model$DATA)))
  }
  if (verbose && length(nm_model$TABLE) > 1) {
    message(sprintf("Found %d $TABLE blocks - using the first.", length(nm_model$TABLE)))
  }

  nm_model$DATA  <- sprintf("$DATA %s.csv IGNORE@", run_id)
  nm_model$TABLE <- sprintf(
    "$TABLE ID EVID MDV TIME DV PRED CWRES ONEHEADER NOPRINT FILE = %s.tab",
    run_id
  )

  new_mod_file <- paste0(run_id, ".mod")
  if (verbose) message("Writing modified model file for iterative runs.")
  mod_tmp_file <- paste0(new_mod_file, ".tmp")
  nm_write_model(nm_model, mod_tmp_file, overwrite = TRUE)
  file.rename(mod_tmp_file, new_mod_file)

  # -- Iterative removal loop ---------------------------------------------------
  for (i in seq_len(n)) {
    csv_path     <- paste0(run_id, ".csv")
    csv_tmp_path <- paste0(csv_path, ".tmp")
    writeLines(paste0("@", paste(names(data), collapse = ",")), csv_tmp_path)
    utils::write.table(data, csv_tmp_path, append = TRUE, sep = ",",
                       col.names = FALSE, row.names = FALSE, quote = FALSE)
    file.rename(csv_tmp_path, csv_path)

    if (verbose) {
      message(sprintf(
        "\n-- Iteration %d of %d (dataset rows: %d) ----------------------",
        i, n, nrow(data)
      ))
    }

    fit_dir      <- paste0("iteration_", run_id, "_", i - 1L)
    exit_execute <- system(
      command       = paste0("execute ", new_mod_file, " --dir=", fit_dir),
      ignore.stdout = !verbose,
      ignore.stderr = !verbose
    )
    if (exit_execute != 0L) {
      stop(sprintf(
        "PsN `execute` failed (exit code %d) at iteration %d. Check output in %s.",
        exit_execute, i, fit_dir
      ))
    }
    exit_sumo <- system(
      command       = paste0("sumo ", run_id, ".lst"),
      ignore.stdout = !verbose,
      ignore.stderr = !verbose
    )
    if (exit_sumo != 0L) {
      warning(sprintf(
        "PsN `sumo` returned non-zero exit code %d at iteration %d.",
        exit_sumo, i
      ))
    }

    # Parameter estimates
    pars   <- nm_read_pars(run_id)
    phifile <- paste0(run_id, ".phi")
    if (!file.exists(phifile)) stop(sprintf("PHI file not found: %s", phifile))

    OFV <- utils::read.table(phifile, skip = 1, header = TRUE)
    par <- data.frame(pars) |>
      dplyr::mutate(
        ITERATION = i,
        nOFV = sum(OFV$OBJ) / (total_dv - i)
      )
    par_list[[i]] <- par

    phi <- utils::read.table(phifile, skip = 1, header = TRUE) |>
      dplyr::select("ID", "OBJ") |>
      dplyr::mutate(ITERATION = i)
    phi_list[[i]] <- phi

    # Results table
    tabfile <- paste0(run_id, ".tab")
    if (!file.exists(tabfile)) stop(sprintf("Table file not found: %s", tabfile))

    result <- tryCatch(
      suppressWarnings(vpc::read_table_nm(tabfile)),
      error = function(e) {
        if (verbose) message("vpc::read_table_nm failed; falling back to read.table.")
        utils::read.table(tabfile, header = TRUE, skip = 1)
      }
    )

    # Observation with the highest |CWRES|
    result_rem <- result |>
      dplyr::mutate(.row_idx = dplyr::row_number()) |>
      dplyr::filter(.data$MDV == 0) |>
      dplyr::arrange(dplyr::desc(abs(.data$CWRES))) |>
      dplyr::slice(1L) |>
      dplyr::mutate(across(c("TIME", "DV", "ID"), as.numeric))

    if (verbose) {
      message(sprintf(
        "Removing: ID=%s, TIME=%.3f, DV=%.4f, CWRES=%.4f",
        result_rem$ID, result_rem$TIME, result_rem$DV, result_rem$CWRES
      ))
    }

    removed_observations <- rbind(
      removed_observations,
      data.frame(
        ID        = result_rem$ID,
        TIME      = result_rem$TIME,
        DV        = result_rem$DV,
        CWRES     = result_rem$CWRES,
        ITERATION = i,
        stringsAsFactors = FALSE
      )
    )

    # Remove row from dataset for next iteration
    current_nm_data <- utils::read.csv(paste0(run_id, ".csv"))
    row_idx <- result_rem$.row_idx
    if (is.na(row_idx) || row_idx < 1L || row_idx > nrow(current_nm_data)) {
      stop("Row index for removal is out of range. Check dataset and model alignment.")
    }
    data <- current_nm_data[-row_idx, ]

    # RMSE
    rmse_val <- result |>
      dplyr::filter(.data$MDV == 0) |>
      dplyr::mutate(rmse = .calculate_nrmse(.data$DV, .data$PRED)) |>
      dplyr::slice(1L) |>
      dplyr::mutate(ITERATION = i) |>
      dplyr::select("rmse", "ITERATION")
    rmse_list[[i]] <- rmse_val

    rem_list[[i]] <- result_rem
  }

  # -- Assemble results ---------------------------------------------------------
  final <- list(
    par  = dplyr::bind_rows(par_list),
    rem  = dplyr::bind_rows(rem_list),
    phi  = dplyr::bind_rows(phi_list),
    rmse = dplyr::bind_rows(rmse_list) |> dplyr::distinct(.data$rmse, .data$ITERATION)
  )

  # Parameter labels parsed from the .mod file (used for plot labels).
  # Use the original file when available; fall back to the rewritten run file.
  final$param_labels <- tryCatch(
    nm_parse_param_labels(if (!is.null(labels_src)) labels_src else new_mod_file),
    error = function(e) list(thetas = character(0), omegas = character(0))
  )

  # -- Stability check ----------------------------------------------------------
  if (stability_check) {
    final$stability <- check_model_stability(final, verbose = FALSE)
  }

  # -- Save results -------------------------------------------------------------
  setwd(original_wd)  # move back before saving (on.exit will be a no-op for wd)
  if (save_results) {
    rds_path <- paste0(run_id, "_results.RDS")
    saveRDS(final, rds_path)
    if (verbose) message(sprintf("Results saved to: %s", rds_path))
  }

  final
}
