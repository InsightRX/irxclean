# =============================================================================
# Shared test helpers and fixtures
# =============================================================================

#' Load the bundled example results (from inst/testdata)
example_results <- function() {
  rds <- system.file("extdata", "busulfan_results.RDS", package = "irxclean")
  if (!nzchar(rds)) {
    rds <- file.path(
      system.file(package = "irxclean"), "..", "..", "inst", "extdata",
      "busulfan_results.RDS"
    )
  }
  # Fall back: read directly from source during development
  if (!file.exists(rds)) {
    rds <- file.path(
      rprojroot::find_package_root_file(), "inst", "extdata",
      "busulfan_results.RDS"
    )
  }
  readRDS(rds)
}

#' Create a minimal synthetic PK dataset for testing
minimal_pk_data <- function(n_subjects = 10, n_obs = 5) {
  set.seed(42)
  n_rows_per_subj <- n_obs + 1L  # 1 dose row + n_obs observation rows
  ids <- rep(seq_len(n_subjects), each = n_rows_per_subj)
  # Dose at time 0, observations at times 1 .. n_obs
  times_per_subj <- c(0, seq_len(n_obs))
  times <- rep(times_per_subj, times = n_subjects)
  evid  <- rep(c(1L, rep(0L, n_obs)), times = n_subjects)
  mdv   <- ifelse(evid == 1L, 1L, 0L)
  dv    <- ifelse(evid == 0L, stats::rexp(length(ids)) * 100, 0)
  data.frame(
    ID   = ids,
    TIME = times,
    DV   = dv,
    EVID = evid,
    MDV  = mdv,
    AMT  = ifelse(evid == 1L, 100, 0),
    AGE  = rep(round(stats::runif(n_subjects, 18, 80)), each = n_rows_per_subj),
    WT   = rep(round(stats::runif(n_subjects, 40, 120), 1), each = n_rows_per_subj),
    SEX  = rep(sample(c(0L, 1L), n_subjects, replace = TRUE), each = n_rows_per_subj),
    stringsAsFactors = FALSE
  )
}

#' Create a minimal synthetic results list for testing (no NONMEM required)
minimal_results <- function(n_iter = 5L) {
  set.seed(1L)
  par_list <- lapply(seq_len(n_iter), function(i) {
    data.frame(
      THETA1 = 10 + stats::rnorm(1, 0, 0.1 * i),
      THETA2 = 50 + stats::rnorm(1, 0, 0.2 * i),
      OMEGA1 = 0.05 + stats::rnorm(1, 0, 0.001 * i),
      OBJ    = -1000 + i * 2,
      ITERATION = i
    )
  })
  rem_list <- lapply(seq_len(n_iter), function(i) {
    data.frame(ID = i, TIME = i * 2.0, DV = 100 - i * 5, CWRES = 4 - i * 0.3, ITERATION = i)
  })
  phi_list <- lapply(seq_len(n_iter), function(i) {
    data.frame(
      ID  = seq_len(20L),
      OBJ = stats::rnorm(20L, -50, 5),
      ITERATION = i
    )
  })
  rmse_list <- lapply(seq_len(n_iter), function(i) {
    data.frame(rmse = 20 - i * 0.5, ITERATION = i)
  })
  list(
    par  = do.call(rbind, par_list),
    rem  = do.call(rbind, rem_list),
    phi  = do.call(rbind, phi_list),
    rmse = do.call(rbind, rmse_list)
  )
}
