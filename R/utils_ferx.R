# =============================================================================
# Internal ferx utility functions
# Lightweight wrappers around the ferx package for use as an alternative
# estimation engine in remove_erroneous_obs().
# =============================================================================

#' Run a single ferx iteration and return standardized results
#'
#' @param model_path Path to the `.ferx` model file.
#' @param data Data frame to fit.
#' @param run_id Character. Run identifier (used for file naming).
#' @param iteration Integer. Current iteration number.
#' @param method Character. Estimation method passed to `ferx::ferx_fit()`.
#' @param verbose Logical. Print progress?
#'
#' @return A named list with elements: `sdtab`, `theta`, `omega`, `sigma`,
#'   `ofv`, `individual_obj`.
#'
#' @keywords internal
ferx_run_iteration <- function(model_path, data, run_id, iteration,
                               method = "focei", verbose = TRUE) {
  if (!requireNamespace("ferx", quietly = TRUE)) {
    stop("Package 'ferx' is required when engine = \"ferx\". Please install it.")
  }

  if (verbose) {
    message(sprintf("  ferx: fitting iteration %d ...", iteration))
  }

  # ferx::ferx_fit() expects a file path for data, not a data.frame
  data_path <- tempfile(
    pattern = sprintf("%s_iter%d_", run_id, iteration),
    fileext = ".csv"
  )
  ferx_write_data(data, data_path)

  fit <- ferx::ferx_fit(
    model      = model_path,
    data       = data_path,
    method     = method,
    covariance = FALSE,
    verbose    = verbose
  )

  # Individual OBJ: one value per subject from EBE_OFV column in sdtab
  sdtab <- fit$sdtab
  if ("EBE_OFV" %in% names(sdtab) && "ID" %in% names(sdtab)) {
    individual_obj <- stats::aggregate(
      sdtab[, "EBE_OFV", drop = FALSE],
      by = list(ID = sdtab$ID),
      FUN = function(x) x[1L]
    )
    names(individual_obj) <- c("ID", "OBJ")
  } else {
    individual_obj <- data.frame(ID = integer(0), OBJ = numeric(0))
  }

  list(
    sdtab          = sdtab,
    theta          = fit$theta,
    omega          = fit$omega,
    sigma          = fit$sigma,
    ofv            = fit$ofv,
    individual_obj = individual_obj
  )
}


#' Convert ferx fit parameters to a flat named list
#'
#' Normalizes `fit$theta`, `fit$omega`, `fit$sigma` into the same column-name
#' format that [nm_read_pars()] produces (THETA1, THETA2, ..., OMEGA.1.1, ...,
#' SIGMA.1.1, ..., OBJ).
#'
#' @param fit_result List returned by [ferx_run_iteration()].
#'
#' @return A named list of scalar parameter values.
#'
#' @keywords internal
ferx_read_pars <- function(fit_result) {
  pars <- list()

  # Thetas
  theta <- fit_result$theta
  if (is.numeric(theta)) {
    nms <- names(theta)
    for (i in seq_along(theta)) {
      key <- if (!is.null(nms) && nzchar(nms[i])) {
        paste0("THETA", i)
      } else {
        paste0("THETA", i)
      }
      pars[[key]] <- theta[i]
    }
  }

  # Omega (matrix -> diagonal elements)
  omega <- fit_result$omega
  if (is.matrix(omega)) {
    for (i in seq_len(nrow(omega))) {
      for (j in seq_len(i)) {
        key <- paste0("OMEGA.", i, ".", j)
        pars[[key]] <- omega[i, j]
      }
    }
  } else if (is.numeric(omega)) {
    for (i in seq_along(omega)) {
      pars[[paste0("OMEGA.", i, ".", i)]] <- omega[i]
    }
  }

  # Sigma (matrix -> elements)
  sigma <- fit_result$sigma
  if (is.matrix(sigma)) {
    for (i in seq_len(nrow(sigma))) {
      for (j in seq_len(i)) {
        key <- paste0("SIGMA.", i, ".", j)
        pars[[key]] <- sigma[i, j]
      }
    }
  } else if (is.numeric(sigma)) {
    for (i in seq_along(sigma)) {
      pars[[paste0("SIGMA.", i, ".", i)]] <- sigma[i]
    }
  }

  # OBJ
  pars$OBJ <- fit_result$ofv

  pars
}


#' Parse parameter labels from a .ferx model file
#'
#' Reads the `[parameters]` block of a `.ferx` file and extracts named
#' theta/omega/sigma labels.  Returns the same structure as
#' [nm_parse_param_labels()].
#'
#' @param model_path Path to the `.ferx` model file.
#'
#' @return A named list with elements `thetas`, `omegas`, and `sigmas`, each a
#'   named character vector mapping column names to labels.
#'
#' @keywords internal
ferx_parse_param_labels <- function(model_path) {
  if (!file.exists(model_path)) {
    return(list(thetas = character(0), omegas = character(0), sigmas = character(0)))
  }

  lines <- readLines(model_path, warn = FALSE)

  theta_labels <- character(0)
  omega_labels <- character(0)
  sigma_labels <- character(0)

  # Find [parameters] block
  param_start <- grep("^\\s*\\[parameters\\]", lines, ignore.case = TRUE)
  if (length(param_start) == 0L) {
    return(list(thetas = theta_labels, omegas = omega_labels, sigmas = sigma_labels))
  }

  # Block ends at next section header or end of file
  section_headers <- grep("^\\s*\\[", lines)
  next_section <- section_headers[section_headers > param_start[1L]]
  end_line <- if (length(next_section) > 0L) next_section[1L] - 1L else length(lines)

  param_lines <- lines[(param_start[1L] + 1L):end_line]

  theta_idx <- 0L
  omega_idx <- 0L
  sigma_idx <- 0L

  for (line in param_lines) {
    line <- trimws(line)
    if (!nzchar(line) || grepl("^#", line) || grepl("^\\[", line)) next

    # Extract label from inline comment: name = value ; Label
    comment_match <- regmatches(line, regexpr(";\\s*(.+)$", line, perl = TRUE))
    label <- if (length(comment_match) > 0L) {
      trimws(sub("^;\\s*", "", comment_match))
    } else {
      NULL
    }

    # Determine parameter type from name prefix
    name_part <- trimws(sub("\\s*[=;].*", "", line))
    name_lower <- tolower(name_part)

    if (grepl("^theta|^tv|^pop", name_lower, ignore.case = TRUE)) {
      theta_idx <- theta_idx + 1L
      if (!is.null(label) && nzchar(label)) {
        theta_labels[paste0("THETA", theta_idx)] <- label
      }
    } else if (grepl("^omega|^iiv|^eta", name_lower, ignore.case = TRUE)) {
      omega_idx <- omega_idx + 1L
      if (!is.null(label) && nzchar(label)) {
        omega_labels[paste0("OMEGA.", omega_idx, ".", omega_idx)] <- label
      }
    } else if (grepl("^sigma|^err|^eps", name_lower, ignore.case = TRUE)) {
      sigma_idx <- sigma_idx + 1L
      if (!is.null(label) && nzchar(label)) {
        sigma_labels[paste0("SIGMA.", sigma_idx, ".", sigma_idx)] <- label
      }
    }
  }

  list(thetas = theta_labels, omegas = omega_labels, sigmas = sigma_labels)
}


#' Write a data frame as CSV for ferx input
#'
#' @param data Data frame to write.
#' @param path Output file path.
#'
#' @keywords internal
ferx_write_data <- function(data, path) {
  utils::write.csv(data, file = path, row.names = FALSE, quote = FALSE)
}
