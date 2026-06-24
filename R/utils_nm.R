# =============================================================================
# Internal NONMEM utility functions
# These are lightweight parsers / writers that do not depend on any external
# InsightRX package.  Adapted from irxnlme (InsightRX).
# =============================================================================

#' Parse a NONMEM model file into a named list of code blocks
#'
#' @param modelfile Path to the NONMEM `.mod` file.
#' @param as_block If `TRUE` (default `FALSE`), each block is returned as a
#'   single collapsed string rather than a character vector of lines.
#' @param code Character string of NONMEM code (alternative to `modelfile`).
#'
#' @return A named list of class `c("NONMEM", "list")`.  Names correspond to
#'   NONMEM block identifiers without the leading `$` (e.g., `"DATA"`,
#'   `"THETA"`).  If the same block appears multiple times the entries are
#'   concatenated.
#'
#' @keywords internal
nm_read_model <- function(modelfile = NULL, as_block = FALSE, code = NULL) {
  if (is.null(modelfile) && is.null(code)) {
    stop("Please specify a NONMEM modelfile or NONMEM code.")
  }
  nm_txt <- if (!is.null(code)) {
    code
  } else {
    if (!file.exists(modelfile)) {
      stop(paste0("NONMEM modelfile (", modelfile, ") not found."))
    }
    readChar(modelfile, file.info(modelfile)$size)
  }

  nm_lines <- stringr::str_split(nm_txt, "\\n")[[1]]
  # Strip leading whitespace / tabs from each line
  nm_lines <- stringr::str_replace_all(nm_lines, "^[\\s\\t]*", "")

  block_idx <- which(stringr::str_detect(nm_lines, "^\\$"))
  if (length(block_idx) == 0) {
    stop("No NONMEM code blocks detected.")
  }
  block_idx <- c(block_idx, length(nm_lines) + 1L)

  obj <- list()
  for (i in seq_len(length(block_idx) - 1L)) {
    block_id <- stringr::str_replace(
      stringr::str_split(nm_lines[block_idx[i]], "\\s")[[1]][1],
      "^\\$", ""
    )
    block_tmp <- nm_lines[block_idx[i]:(block_idx[i + 1L] - 1L)]
    if (as_block) {
      block_tmp <- stringr::str_c(block_tmp, collapse = "\n")
    }
    if (is.null(obj[[block_id]])) {
      obj[[block_id]] <- block_tmp
    } else {
      obj[[block_id]] <- c(obj[[block_id]], block_tmp)
    }
  }
  class(obj) <- c("NONMEM", "list")
  obj
}


#' Write a NONMEM model object to a `.mod` file
#'
#' @param model A NONMEM model object created by [nm_read_model()].
#' @param modelfile Output file path.
#' @param overwrite Overwrite an existing file?  Default `FALSE`.
#'
#' @return Invisibly returns \code{modelfile} (the output file path).
#'
#' @keywords internal
nm_write_model <- function(model = NULL, modelfile = NULL, overwrite = FALSE) {
  if (is.null(modelfile)) stop("Please specify an output NONMEM modelfile.")
  if (is.null(model))     stop("Please specify a NONMEM model object.")
  if (file.exists(modelfile) && !overwrite) {
    stop("Output file exists and `overwrite = FALSE`.")
  }
  if (!"NONMEM" %in% class(model)) {
    stop("Object does not appear to be a valid NONMEM model. Use nm_read_model().")
  }
  # Enforce standard block ordering ($PROB / $DATA / $INPUT first).
  # nm_read_model strips "$" so the key is the first token after "$".
  # Models using "$PROBLEM" produce key "PROBLEM"; "$PROB" produces "PROB".
  # Both are included so either form is placed first.
  # NONMEM requires $DATA before $INPUT.
  header <- c("PROB", "PROBLEM", "DATA", "INPUT", "ABBR")
  ordered_names <- c(header, setdiff(names(model), header))
  model <- model[ordered_names[ordered_names %in% names(model)]]
  writeLines(unlist(model), con = modelfile)
  invisible(modelfile)
}


# Sentinel value NONMEM writes in the ITERATION column of .ext files to mark
# the row containing the final parameter estimates.
NM_FINAL_ITER <- -1000000000L

#' Read population parameter estimates from a NONMEM `.ext` file
#'
#' Returns the final estimates (iteration \code{-1000000000}).
#'
#' @param model Model name (character, e.g. `"run1"`) or run number (numeric).
#'   The `.ext` suffix is appended automatically if absent.
#'
#' @return A named list of parameter values.
#'
#' @keywords internal
nm_read_pars <- function(model) {
  if (is.numeric(model)) model <- paste0("run", model)
  model <- stringr::str_replace(model, "\\.mod$", "")
  parfile <- paste0(model, ".ext")
  if (!file.exists(parfile)) stop(paste0("Parameter file not found: ", parfile))
  pars <- utils::read.table(parfile, skip = 1, header = TRUE)
  names(pars) <- stringr::str_replace_all(names(pars), "\\.$", "")
  pars <- as.list(pars[pars$ITERATION == NM_FINAL_ITER, , drop = FALSE])
  pars$ITERATION <- NULL
  pars
}


#' Parse THETA and OMEGA parameter labels from NONMEM .mod file comments
#'
#' Extracts human-readable labels from inline comments in the \code{$THETA}
#' and \code{$OMEGA} blocks of a NONMEM control stream.  Labels are used
#' automatically as panel titles in \code{plot_removal_metrics()}.
#'
#' @param mod_file Path to the NONMEM \code{.mod} file.
#'
#' @return A named list with elements \code{thetas} and \code{omegas}, each a
#'   named character vector mapping column names (e.g. \code{"THETA1"},
#'   \code{"OMEGA.3.3"}) to human-readable labels.  Returns
#'   \code{list(thetas = character(0), omegas = character(0))} if no labels
#'   are found.
#'
#' @section Labelling your .mod file for irxclean:
#'
#' To get readable panel titles in irxclean plots, add inline comments to the
#' \code{$THETA} and \code{$OMEGA} blocks of your NONMEM control stream.
#'
#' \strong{Format:}
#' \preformatted{
#'   ; <index>. <Label>          -- single parameter
#'   ; <start>-<end>. <Label>    -- range (e.g. IOV block)
#' }
#'
#' The number(s) before the dot correspond to the 1-based parameter index
#' (the same numbering NONMEM uses in the \code{.ext} file).  The dot and
#' surrounding spaces are flexible -- the parser accepts \code{; 1. CL},
#' \code{; 1.CL}, and \code{; 1 CL}.
#'
#' \strong{THETA example:}
#' \preformatted{
#' $THETA
#'   (0, 11.5555) ; 1. CL
#'   (0,  9.2)    ; 2. V1
#'   (0,  3.1)    ; 3. KA
#' }
#'
#' \strong{OMEGA example (IIV + IOV):}
#' \preformatted{
#' $OMEGA
#'   0.09         ; 1. IIV CL
#'   0.04         ; 2. IIV V1
#'
#' $OMEGA BLOCK(1)
#'   0.03         ; 3. IOV CL
#' $OMEGA BLOCK(1) SAME   ; -- no label needed, auto-excluded from plots
#' $OMEGA BLOCK(1) SAME
#' }
#'
#' \strong{Rules:}
#' \itemize{
#'   \item \strong{THETA}: label every row you want named; unlabelled rows
#'     keep their default name (e.g. \code{THETA4}).
#'   \item \strong{OMEGA}: label only the \emph{first} element of each unique
#'     term.  \code{BLOCK(1) SAME} lines are detected automatically and always
#'     excluded from plots regardless of any comment on that line.
#'   \item Range notation (\code{; 3-7. IOV CL}) is supported; only the first
#'     element (\code{OMEGA.3.3}) is labelled and plotted.
#'   \item The parser is case-insensitive for the \code{SAME} keyword and
#'     tolerates varied spacing around the dash and dot.
#' }
#'
#' @keywords internal
nm_parse_param_labels <- function(mod_file) {
  model <- nm_read_model(mod_file)

  # Helper: extract (start, end, label) from a commented line.
  # Accepts formats:
  #   "; 1. CL"   "; 1.CL"   "; 1 CL"
  #   "; 3-7. IOV CL"   "; 3 - 7 . IOV CL"
  extract_label <- function(line) {
    m <- regexpr(
      ";\\s*(\\d+)(?:\\s*-\\s*(\\d+))?\\s*\\.?\\s*(.+)$",
      line, perl = TRUE
    )
    if (m == -1L) return(NULL)
    matched <- regmatches(line, m)
    cap <- regmatches(
      matched,
      regexec("(\\d+)(?:\\s*-\\s*(\\d+))?\\s*\\.?\\s*(.+)$", matched)
    )[[1L]]
    if (length(cap) < 4L) return(NULL)
    lbl <- trimws(cap[4L])
    if (!nzchar(lbl)) return(NULL)
    list(
      start = as.integer(cap[2L]),
      end   = if (nzchar(cap[3L])) as.integer(cap[3L]) else as.integer(cap[2L]),
      label = lbl
    )
  }

  # Helper: detect a BLOCK(n) SAME header line (case-insensitive).
  is_same_line <- function(line) {
    grepl("\\bSAME\\b", line, ignore.case = TRUE) &&
      grepl("^\\$OMEGA", trimws(line), ignore.case = TRUE)
  }

  # -- THETA labels ------------------------------------------------------------
  theta_labels <- character(0)
  if (!is.null(model$THETA)) {
    for (line in model$THETA) {
      info <- extract_label(line)
      if (!is.null(info)) {
        for (n in info$start:info$end) {
          theta_labels[paste0("THETA", n)] <- info$label
        }
      }
    }
  }

  # -- OMEGA labels (diagonal, first element of each unique term only) ---------
  # BLOCK(1) SAME lines are skipped unconditionally: their diagonal value is
  # identical to the preceding BLOCK definition and would produce a redundant
  # plot panel.  The first element of each unique OMEGA term is labelled using
  # only info$start so that range notation ("3-7. IOV CL") labels OMEGA.3.3
  # and leaves OMEGA.4.4 ... OMEGA.7.7 unlabelled (and thus filtered from
  # plots when param_labels are present).
  omega_labels <- character(0)
  if (!is.null(model$OMEGA)) {
    for (line in model$OMEGA) {
      if (is_same_line(line)) next
      info <- extract_label(line)
      if (!is.null(info)) {
        key <- paste0("OMEGA.", info$start, ".", info$start)
        omega_labels[key] <- info$label
      }
    }
  }

  # -- SIGMA labels (diagonal elements) ----------------------------------------
  sigma_labels <- character(0)
  if (!is.null(model$SIGMA)) {
    for (line in model$SIGMA) {
      info <- extract_label(line)
      if (!is.null(info)) {
        key <- paste0("SIGMA.", info$start, ".", info$start)
        sigma_labels[key] <- info$label
      }
    }
  }

  list(thetas = theta_labels, omegas = omega_labels, sigmas = sigma_labels)
}
