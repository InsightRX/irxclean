# =============================================================================
# Iterative Influential Individual Removal
# Based on Karlsson's influential individual approach, extended to iterate
# sequentially (analogous to the CWRES observation-based approach).
# =============================================================================

#' Iteratively detect and remove influential individuals
#'
#' Performs iterative NONMEM-based detection of individuals whose removal most
#' improves overall model fit, as measured by the change in total objective
#' function value (dOFV).  The method extends the influential individual
#' approach of Karlsson by iterating: after each removal the model is
#' re-estimated on the reduced dataset and the next most influential subject is
#' identified from the updated individual OBJ values.
#'
#' At each iteration the algorithm:
#' \enumerate{
#'   \item Runs NONMEM (via PsN \code{execute}) on the current dataset.
#'   \item Reads the individual OBJ values from the \code{.phi} file and
#'     computes \code{sumOFV}.
#'   \item (From iteration 2 onwards) Calculates
#'     \code{dOFV = sumOFV(prev) - sumOFV(curr)} and checks the stopping
#'     criterion.
#'   \item Identifies the subject with the highest individual OBJ (iOFV).
#'   \item Removes all rows for that subject from the dataset.
#'   \item Saves parameter estimates and diagnostics for later plotting.
#' }
#'
#' @param dat File path to the NONMEM-ready \code{.csv} dataset (with or
#'   without the \code{.csv} extension, relative or absolute), \strong{or} an
#'   R data frame that is already loaded in memory.
#' @param mod File path to the NONMEM \code{.mod} control stream (with or
#'   without the \code{.mod} extension, relative or absolute), \strong{or} a
#'   \code{nm_model} list returned by \code{nm_read_model()}.
#' @param run_id Character. Unique run identifier (e.g. \code{"iir_run1"}).
#' @param n Integer or \code{NULL}.  Maximum number of subjects to remove.
#'   Must specify at least one of \code{n} or \code{dofv_threshold}.
#'   If both are supplied the run stops when \emph{either} criterion is met.
#' @param dofv_threshold Numeric or \code{NULL}.  Stop when the change in
#'   sumOFV from removing a subject (\code{sumOFV_before - sumOFV_after})
#'   falls below this value.  A common choice is 3.84 (chi-squared 1 df,
#'   p < 0.05).  When the threshold is not met the last entry in \code{rem}
#'   retains its \code{dOFV} value so the user can inspect it before deciding
#'   whether to apply the removal.  Must specify at least one of \code{n} or
#'   \code{dofv_threshold}.
#' @param id_col Character.  Name of the subject ID column in the dataset.
#'   Default \code{"ID"}.
#' @param verbose Logical.  Print progress messages?  Default \code{TRUE}.
#' @param save_results Logical.  Save an \code{.RDS} results file to disk?
#'   Default \code{TRUE}.
#' @param save_temp_dir Logical.  Copy the temporary working directory to
#'   the current working directory (\code{TRUE}) or delete it after completion
#'   (\code{FALSE}, default)?
#' @param stability_check Logical.  Append model stability metrics to the
#'   returned results?  Default \code{TRUE}.  See \code{check_model_stability()}.
#'
#' @return A list containing:
#' \describe{
#'   \item{par}{Data frame of population parameter estimates per iteration.
#'     ITERATION 1 = baseline (full dataset); ITERATION 2..k+1 = after each
#'     of k subject removals.}
#'   \item{rem}{Data frame of removed subjects with columns: \code{ID},
#'     \code{iOFV} (individual OBJ before removal), \code{sumOFV_before},
#'     \code{sumOFV_after}, \code{dOFV}, and \code{ITERATION} (removal
#'     sequence number, 1-based).}
#'   \item{phi}{Data frame of individual OBJ values per NONMEM run iteration
#'     (columns \code{ID}, \code{OBJ}, \code{ITERATION}).  Removed subjects
#'     are absent from subsequent iterations.}
#'   \item{ofv}{Data frame of \code{sumOFV} and \code{dOFV} per NONMEM run
#'     iteration.  \code{dOFV} at ITERATION 1 is \code{NA} (baseline).}
#'   \item{rmse}{Data frame of normalised RMSE per iteration.}
#'   \item{stability}{(If \code{stability_check = TRUE}) Stability assessment
#'     from \code{check_model_stability()}.}
#'   \item{param_labels}{Named list with elements \code{thetas} and
#'     \code{omegas}: human-readable parameter labels parsed from the
#'     \code{.mod} file comments, used automatically in plots.}
#' }
#'
#' @section Stopping behaviour:
#' When \code{dofv_threshold} is the active stopping criterion, the algorithm
#' runs one additional NONMEM estimation after each removal to obtain the
#' updated \code{sumOFV}.  The iteration that first produces
#' \code{dOFV < dofv_threshold} is recorded in \code{rem} with its
#' \code{dOFV} value; no further removals are performed.  Users should inspect
#' that final entry before deciding whether to apply it.
#'
#' @section Supported input modes:
#' The same three input modes as \code{remove_erroneous_obs()} are supported:
#' bare file name, relative/absolute file path, or in-memory R objects.
#' See \code{remove_erroneous_obs()} for full details.
#'
#' @section PsN requirement:
#' \code{execute} and \code{sumo} must be available on the system \code{PATH}.
#'
#' @seealso \code{remove_erroneous_obs()}, \code{plot_removal_metrics()},
#'   \code{check_model_stability()}, \code{remove_observations_from_data()}
#'
#' @examples
#' \dontrun{
#' # Remove up to 10 subjects
#' results <- remove_influential_individuals(
#'   dat    = "my_data",
#'   mod    = "my_model",
#'   run_id = "iir_run1",
#'   n      = 10
#' )
#'
#' # Stop when dOFV drops below chi-squared threshold (p < 0.05, 1 df)
#' results <- remove_influential_individuals(
#'   dat            = "my_data",
#'   mod            = "my_model",
#'   run_id         = "iir_run1",
#'   dofv_threshold = 3.84
#' )
#'
#' # Both: stop at whichever criterion is hit first
#' results <- remove_influential_individuals(
#'   dat            = "my_data",
#'   mod            = "my_model",
#'   run_id         = "iir_run1",
#'   n              = 20,
#'   dofv_threshold = 3.84
#' )
#'
#' # Inspect removals
#' results$rem
#' results$ofv
#'
#' # Reuse existing plot functions for parameter / RMSE diagnostics
#' plot_removal_metrics(results, metric = "thetas")
#' plot_removal_metrics(results, metric = "rmse")
#' }
#'
#' @export
remove_influential_individuals <- function(
    dat,
    mod,
    run_id,
    n               = NULL,
    dofv_threshold  = NULL,
    id_col          = "ID",
    verbose         = TRUE,
    save_results    = TRUE,
    save_temp_dir   = FALSE,
    stability_check = TRUE
) {
  # -- Input validation ----------------------------------------------------------
  if (!is.character(run_id)) stop("`run_id` must be a character string.")
  if (is.null(n) && is.null(dofv_threshold)) {
    stop("Specify at least one of `n` (subjects to remove) or `dofv_threshold`.")
  }
  if (!is.null(n)) {
    n <- as.integer(n)
    if (is.na(n) || n < 1L) stop("`n` must be a positive integer.")
  }
  if (!is.null(dofv_threshold)) {
    dofv_threshold <- as.numeric(dofv_threshold)
    if (is.na(dofv_threshold)) stop("`dofv_threshold` must be numeric.")
  }

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

  # -- Temporary working directory -----------------------------------------------
  temp_dir <- tempfile(pattern = paste0("irxclean_iir_", run_id, "_"))
  dir.create(temp_dir, recursive = TRUE)

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

  # -- Prepare dataset and model -------------------------------------------------
  if (!is.null(dat_file)) {
    file.copy(from = dat_file, to = basename(dat_file))
    data <- utils::read.csv(basename(dat_file))
  } else {
    data <- dat_obj
  }

  if (!id_col %in% names(data)) {
    stop(sprintf("ID column '%s' not found in dataset.", id_col))
  }

  if (!is.null(mod_file)) {
    mod_base   <- basename(mod_file)
    file.copy(from = mod_file, to = mod_base)
    nm_model   <- nm_read_model(mod_base)
    labels_src <- mod_base
  } else {
    mod_base   <- NULL
    nm_model   <- mod_obj
    labels_src <- NULL
  }

  total_subjects <- length(unique(data[[id_col]]))
  max_removals   <- if (!is.null(n)) n else (total_subjects - 1L)
  if (max_removals < 1L) stop("Not enough subjects to remove.")
  if (max_removals >= total_subjects) {
    stop(sprintf(
      "`n` (%d) must be less than the number of subjects in the dataset (%d).",
      max_removals, total_subjects
    ))
  }

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
  nm_write_model(nm_model, new_mod_file, overwrite = TRUE)
  if (verbose) message("Writing modified model file for iterative runs.")

  # -- Iterative removal loop ----------------------------------------------------
  # iter          : NONMEM run counter (1-based; iter 1 = baseline full-data run)
  # removals_done : subjects removed so far
  # Each entry in rem_list is initially stored after a subject is identified;
  # its sumOFV_after and dOFV are filled in after the subsequent NONMEM run.
  par_list  <- list()
  phi_list  <- list()
  rem_list  <- list()
  rmse_list <- list()

  iter          <- 0L
  removals_done <- 0L
  sumOFV_prev   <- NULL

  repeat {
    iter <- iter + 1L

    # Write current dataset to CSV with IGNORE@ header
    csv_path <- paste0(run_id, ".csv")
    writeLines(paste0("@", paste(names(data), collapse = ",")), csv_path)
    utils::write.table(data, csv_path, append = TRUE, sep = ",",
                       col.names = FALSE, row.names = FALSE, quote = FALSE)

    n_subj_cur <- length(unique(data[[id_col]]))
    if (verbose) {
      if (removals_done == 0L) {
        message(sprintf(
          "\n-- Baseline run (subjects: %d) ----------------------",
          n_subj_cur
        ))
      } else {
        message(sprintf(
          "\n-- Iteration %d/%d  (subjects remaining: %d) -----------",
          removals_done, max_removals, n_subj_cur
        ))
      }
    }

    fit_dir <- sprintf("iteration_%s_%d", run_id, iter - 1L)
    system(
      command       = paste0("execute ", new_mod_file, " --dir=", fit_dir),
      intern        = TRUE,
      ignore.stdout = !verbose,
      ignore.stderr = !verbose
    )
    system(
      command       = paste0("sumo ", run_id, ".lst"),
      intern        = TRUE,
      ignore.stdout = !verbose,
      ignore.stderr = !verbose
    )

    # Read individual OBJ values from phi file
    phifile <- paste0(run_id, ".phi")
    if (!file.exists(phifile)) stop(sprintf("PHI file not found: %s", phifile))
    phi_raw  <- utils::read.table(phifile, skip = 1, header = TRUE)
    sumOFV_i <- sum(phi_raw$OBJ)

    phi_df <- data.frame(
      ID        = phi_raw$ID,
      OBJ       = phi_raw$OBJ,
      ITERATION = iter,
      stringsAsFactors = FALSE
    )
    phi_list[[iter]] <- phi_df

    # Population parameter estimates
    pars   <- nm_read_pars(run_id)
    par_df <- data.frame(pars) |>
      dplyr::mutate(ITERATION = iter)
    par_list[[iter]] <- par_df

    # RMSE from table file
    tabfile <- paste0(run_id, ".tab")
    if (file.exists(tabfile)) {
      result <- tryCatch(
        suppressWarnings(vpc::read_table_nm(tabfile)),
        error = function(e) {
          if (verbose) message("vpc::read_table_nm failed; falling back to read.table.")
          utils::read.table(tabfile, header = TRUE, skip = 1)
        }
      )
      rmse_val <- result |>
        dplyr::filter(.data$MDV == 0) |>
        dplyr::mutate(rmse = .calculate_nrmse(.data$DV, .data$PRED)) |>
        dplyr::slice(1L) |>
        dplyr::mutate(ITERATION = iter) |>
        dplyr::select("rmse", "ITERATION")
      rmse_list[[iter]] <- rmse_val
    }

    # Update dOFV for the most recent removal (now that we have sumOFV_after)
    if (removals_done > 0L) {
      dOFV <- sumOFV_prev - sumOFV_i
      rem_list[[removals_done]]$sumOFV_after <- sumOFV_i
      rem_list[[removals_done]]$dOFV         <- dOFV

      if (verbose) {
        message(sprintf(
          "  dOFV for removal %d (ID=%s): %.3f  [sumOFV: %.2f -> %.2f]",
          removals_done,
          rem_list[[removals_done]]$ID,
          dOFV,
          sumOFV_prev,
          sumOFV_i
        ))
      }

      # Check dOFV stopping criterion
      if (!is.null(dofv_threshold) && dOFV < dofv_threshold) {
        if (verbose) {
          message(sprintf(
            "Stopping: dOFV %.3f < threshold %.3f after removal %d (ID=%s).",
            dOFV, dofv_threshold, removals_done,
            rem_list[[removals_done]]$ID
          ))
        }
        break
      }
    }

    # Check n stopping criterion (all requested removals done + dOFV recorded)
    if (removals_done >= max_removals) {
      if (verbose) {
        message(sprintf(
          "Stopping: reached maximum removals (%d).", max_removals
        ))
      }
      break
    }

    # Identify subject with highest individual OBJ
    worst_idx  <- which.max(phi_raw$OBJ)
    worst_id   <- phi_raw$ID[worst_idx]
    worst_iofv <- phi_raw$OBJ[worst_idx]

    removals_done <- removals_done + 1L
    rem_list[[removals_done]] <- data.frame(
      ID            = worst_id,
      iOFV          = worst_iofv,
      sumOFV_before = sumOFV_i,
      sumOFV_after  = NA_real_,
      dOFV          = NA_real_,
      ITERATION     = removals_done,
      stringsAsFactors = FALSE
    )

    if (verbose) {
      message(sprintf(
        "  Flagging subject ID=%s  (iOFV=%.3f, rank 1 of %d remaining)",
        worst_id, worst_iofv, n_subj_cur
      ))
    }

    # Remove all rows for that subject
    data        <- data[data[[id_col]] != worst_id, ]
    sumOFV_prev <- sumOFV_i
  }

  # -- Assemble results ----------------------------------------------------------
  phi_all  <- dplyr::bind_rows(phi_list)
  par_all  <- dplyr::bind_rows(par_list)
  rem_all  <- dplyr::bind_rows(rem_list)
  rmse_all <- dplyr::bind_rows(rmse_list) |>
    dplyr::distinct(.data$rmse, .data$ITERATION)

  # sumOFV / dOFV summary table (one row per NONMEM run)
  iters_done <- sort(unique(phi_all$ITERATION))
  sumOFV_vec <- vapply(iters_done, function(k) {
    sum(phi_all$OBJ[phi_all$ITERATION == k])
  }, numeric(1L))
  ofv_tab <- data.frame(
    ITERATION = iters_done,
    sumOFV    = sumOFV_vec,
    stringsAsFactors = FALSE
  )
  dOFV_vec        <- c(NA_real_, -diff(sumOFV_vec))
  ofv_tab$dOFV    <- dOFV_vec

  final <- list(
    par  = par_all,
    rem  = rem_all,
    phi  = phi_all,
    ofv  = ofv_tab,
    rmse = rmse_all
  )

  # Parameter labels
  final$param_labels <- tryCatch(
    nm_parse_param_labels(if (!is.null(labels_src)) labels_src else new_mod_file),
    error = function(e) list(thetas = character(0), omegas = character(0))
  )

  # Stability check
  if (stability_check) {
    final$stability <- check_model_stability(final, verbose = FALSE)
  }

  # Save and return
  setwd(original_wd)
  if (save_results) {
    rds_path <- paste0(run_id, "_results.RDS")
    saveRDS(final, rds_path)
    if (verbose) message(sprintf("Results saved to: %s", rds_path))
  }

  final
}


#' Plot dOFV trajectory from iterative influential individual removal
#'
#' Creates a waterfall/trajectory plot of the change in total OFV
#' (\code{dOFV = sumOFV_before - sumOFV_after}) for each subject removed
#' by \code{remove_influential_individuals()}.  An optional threshold line
#' marks the stopping criterion when \code{dofv_threshold} was used.
#'
#' @param results A results list from \code{remove_influential_individuals()}.
#' @param dofv_threshold Numeric or \code{NULL}.  If supplied, a horizontal
#'   reference line is drawn at this value (e.g. 3.84).  Default \code{NULL}.
#' @param label_ids Logical.  Annotate each point with the removed subject ID?
#'   Default \code{TRUE}.
#'
#' @return A \code{ggplot} object.
#'
#' @seealso \code{remove_influential_individuals()}
#'
#' @examples
#' \dontrun{
#' results <- remove_influential_individuals(
#'   dat = "my_data", mod = "my_model",
#'   run_id = "iir_run1", dofv_threshold = 3.84
#' )
#' plot_influential_dofv(results, dofv_threshold = 3.84)
#' }
#'
#' @export
plot_influential_dofv <- function(
    results,
    dofv_threshold = NULL,
    label_ids      = TRUE
) {
  rem <- results$rem
  if (is.null(rem) || nrow(rem) == 0L) {
    stop("`results$rem` is empty - no subjects were removed.")
  }
  if (!"dOFV" %in% names(rem)) {
    stop("`results$rem` does not contain a `dOFV` column.")
  }

  rem$ID <- as.character(rem$ID)

  p <- ggplot2::ggplot(
    rem,
    ggplot2::aes(.data$ITERATION, .data$dOFV)
  ) +
    ggplot2::geom_col(fill = "#003B4F", width = 0.6, alpha = 0.85) +
    ggplot2::geom_point(colour = "#003B4F", size = 2)

  if (!is.null(dofv_threshold)) {
    p <- p + ggplot2::geom_hline(
      yintercept = dofv_threshold,
      linetype   = 2,
      colour     = "#AD0000",
      linewidth  = 0.7
    ) +
      ggplot2::annotate(
        "text",
        x      = max(rem$ITERATION) + 0.3,
        y      = dofv_threshold,
        label  = sprintf("threshold = %.2f", dofv_threshold),
        hjust  = 1,
        vjust  = -0.4,
        size   = 3,
        colour = "#AD0000"
      )
  }

  if (label_ids) {
    p <- p + ggplot2::geom_text(
      ggplot2::aes(label = .data$ID),
      vjust  = -0.6,
      size   = 2.8,
      colour = "#5E5E5E"
    )
  }

  p +
    .clean_theme() +
    ggplot2::ggtitle(
      "dOFV per influential individual removal",
      subtitle = paste0(
        "dOFV = sumOFV(before removal) - sumOFV(after removal)  |  ",
        "Higher = more influential"
      )
    ) +
    ggplot2::xlab("Removal sequence") +
    ggplot2::ylab("\u0394 OFV") +
    ggplot2::scale_x_continuous(breaks = seq_len(nrow(rem)))
}
