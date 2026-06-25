# =============================================================================
# Integration tests for ferx engine
# These tests require the ferx package and exercise the full ferx code path
# in remove_erroneous_obs() and the ferx utility functions.
# =============================================================================

# -- Helper: minimal ferx model and dataset -----------------------------------

#' Create a minimal ferx model file for testing
#' One-compartment IV bolus, IIV on CL and V, proportional error.
write_test_ferx_model <- function(path) {
  writeLines(c(
    "[parameters]",
    "  theta TVCL(4.0, 0.1, 100.0)  ; Clearance",
    "  theta TVV(15.0, 1.0, 500.0)   ; Volume",
    "",
    "  omega ETA_CL ~ 0.09  ; IIV CL",
    "  omega ETA_V  ~ 0.04  ; IIV V",
    "",
    "  sigma PROP_ERR ~ 0.1 (sd)  ; Proportional error",
    "",
    "[individual_parameters]",
    "  CL = TVCL * exp(ETA_CL)",
    "  V  = TVV  * exp(ETA_V)",
    "",
    "[structural_model]",
    "  ode(states=[central])",
    "",
    "[odes]",
    "  d/dt(central) = -(CL/V) * central",
    "",
    "[scaling]",
    "  y = central / V",
    "",
    "[error_model]",
    "  DV ~ proportional(PROP_ERR)",
    "",
    "[fit_options]",
    "  method  = focei",
    "  maxiter = 300"
  ), path)
  invisible(path)
}

#' Create a minimal PK dataset compatible with ferx
#' 5 subjects, IV bolus, 4 obs each -> 25 rows total
make_ferx_test_data <- function() {
  set.seed(123)
  n_subj <- 5L
  n_obs  <- 4L
  rows <- lapply(seq_len(n_subj), function(id) {
    cl <- 4.0 * exp(rnorm(1, 0, 0.3))
    v  <- 15.0 * exp(rnorm(1, 0, 0.2))
    dose_amt <- 100
    obs_times <- c(0.5, 2, 4, 8)
    conc <- (dose_amt / v) * exp(-cl / v * obs_times) *
      (1 + rnorm(n_obs, 0, 0.1))
    conc[conc < 0] <- 0.01
    rbind(
      data.frame(ID = id, TIME = 0, DV = 0, AMT = dose_amt,
                 EVID = 1L, MDV = 1L, CMT = 1L, RATE = 0),
      data.frame(ID = id, TIME = obs_times, DV = round(conc, 4),
                 AMT = 0, EVID = 0L, MDV = 0L, CMT = 1L, RATE = 0)
    )
  })
  do.call(rbind, rows)
}


# =============================================================================
# ferx_run_iteration — full integration
# =============================================================================

test_that("ferx_run_iteration returns expected structure", {
  skip_if_not_installed("ferx")

  model_path <- tempfile(fileext = ".ferx")
  write_test_ferx_model(model_path)
  dat <- make_ferx_test_data()
  on.exit(unlink(model_path))

  result <- ferx_run_iteration(
    model_path = model_path,
    data       = dat,
    run_id     = "test_int",
    iteration  = 1L,
    method     = "focei",
    verbose    = FALSE
  )

  expect_type(result, "list")
  expect_named(result,
    c("sdtab", "theta", "omega", "sigma", "ofv", "individual_obj"),
    ignore.order = TRUE
  )

  # sdtab should have observation rows with CWRES
  expect_s3_class(result$sdtab, "data.frame")
  expect_true("CWRES" %in% names(result$sdtab))
  expect_true("ID" %in% names(result$sdtab))
  expect_true("DV" %in% names(result$sdtab))
  expect_true(nrow(result$sdtab) > 0L)

  # Parameters should be numeric

  expect_true(is.numeric(result$theta))
  expect_true(length(result$theta) >= 2L)
  expect_true(is.numeric(result$ofv))

  # Individual OBJ: one row per subject
  expect_s3_class(result$individual_obj, "data.frame")
  expect_equal(nrow(result$individual_obj), 5L)
  expect_named(result$individual_obj, c("ID", "OBJ"))
})

test_that("ferx_run_iteration with verbose prints progress", {
  skip_if_not_installed("ferx")

  model_path <- tempfile(fileext = ".ferx")
  write_test_ferx_model(model_path)
  dat <- make_ferx_test_data()
  on.exit(unlink(model_path))

  expect_message(
    ferx_run_iteration(model_path, dat, "test_v", 1L, verbose = TRUE),
    "ferx: fitting iteration 1"
  )
})


# =============================================================================
# ferx_read_pars — with real fit output
# =============================================================================

test_that("ferx_read_pars produces correct keys from real fit", {
  skip_if_not_installed("ferx")

  model_path <- tempfile(fileext = ".ferx")
  write_test_ferx_model(model_path)
  dat <- make_ferx_test_data()
  on.exit(unlink(model_path))

  fit <- ferx_run_iteration(model_path, dat, "test_rp", 1L, verbose = FALSE)
  pars <- ferx_read_pars(fit)

  expect_true("THETA1" %in% names(pars))
  expect_true("THETA2" %in% names(pars))
  expect_true("OBJ" %in% names(pars))
  # Should have at least one OMEGA
  omega_keys <- grep("^OMEGA", names(pars), value = TRUE)
  expect_true(length(omega_keys) >= 1L)
  # All values should be numeric scalars
  expect_true(all(vapply(pars, is.numeric, logical(1L))))
})


# =============================================================================
# remove_erroneous_obs with engine = "ferx" — end-to-end
# =============================================================================

test_that("remove_erroneous_obs works end-to-end with ferx engine", {
  skip_if_not_installed("ferx")

  model_path <- tempfile(fileext = ".ferx")
  write_test_ferx_model(model_path)
  dat <- make_ferx_test_data()
  on.exit(unlink(model_path))

  n_remove <- 2L
  results <- remove_erroneous_obs(
    dat             = dat,
    mod             = model_path,
    run_id          = "ferx_e2e",
    n               = n_remove,
    engine          = "ferx",
    ferx_method     = "focei",
    verbose         = FALSE,
    save_results    = FALSE,
    stability_check = TRUE
  )

  # Structure checks
  expect_type(results, "list")
  expect_true(all(c("par", "rem", "phi", "rmse") %in% names(results)))

  # par: one row per iteration
  expect_equal(nrow(results$par), n_remove)
  expect_true("THETA1" %in% names(results$par))
  expect_true("ITERATION" %in% names(results$par))
  expect_true("nOFV" %in% names(results$par))

  # rem: one removed obs per iteration
  expect_equal(nrow(results$rem), n_remove)
  expect_true(all(c("ID", "TIME", "DV", "CWRES") %in% names(results$rem)))

  # phi: individual OBJ per iteration
  expect_true(nrow(results$phi) > 0L)
  expect_true(all(c("ID", "OBJ", "ITERATION") %in% names(results$phi)))

  # rmse: one row per iteration
  expect_equal(nrow(results$rmse), n_remove)
  expect_true("rmse" %in% names(results$rmse))

  # param_labels should be populated from .ferx file
  expect_type(results$param_labels, "list")
  expect_true("thetas" %in% names(results$param_labels))

  # stability should be present

  expect_true("stability" %in% names(results))
})

test_that("remove_erroneous_obs with ferx and data.frame input works", {
  skip_if_not_installed("ferx")

  model_path <- tempfile(fileext = ".ferx")
  write_test_ferx_model(model_path)
  dat <- make_ferx_test_data()
  on.exit(unlink(model_path))

  # Pass data as data.frame (not file path)
  results <- remove_erroneous_obs(
    dat             = dat,
    mod             = model_path,
    run_id          = "ferx_df",
    n               = 1L,
    engine          = "ferx",
    verbose         = FALSE,
    save_results    = FALSE,
    stability_check = FALSE
  )

  expect_equal(nrow(results$rem), 1L)
  expect_true(abs(results$rem$CWRES[1]) > 0)
})

test_that("remove_erroneous_obs ferx rejects non-character mod", {
  skip_if_not_installed("ferx")

  expect_error(
    remove_erroneous_obs(
      dat    = make_ferx_test_data(),
      mod    = list(fake = "model"),
      run_id = "test",
      n      = 1L,
      engine = "ferx"
    ),
    "file path"
  )
})

test_that("remove_erroneous_obs ferx rejects missing model file", {
  skip_if_not_installed("ferx")

  expect_error(
    remove_erroneous_obs(
      dat    = make_ferx_test_data(),
      mod    = "/nonexistent/model.ferx",
      run_id = "test",
      n      = 1L,
      engine = "ferx"
    ),
    "not found"
  )
})

test_that("remove_erroneous_obs ferx with save_results writes RDS", {
  skip_if_not_installed("ferx")

  model_path <- tempfile(fileext = ".ferx")
  write_test_ferx_model(model_path)
  dat <- make_ferx_test_data()
  old_wd <- getwd()
  tmp_out <- tempdir()
  setwd(tmp_out)
  on.exit({ setwd(old_wd); unlink(model_path) })

  results <- remove_erroneous_obs(
    dat          = dat,
    mod          = model_path,
    run_id       = "ferx_save",
    n            = 1L,
    engine       = "ferx",
    verbose      = FALSE,
    save_results = TRUE,
    stability_check = FALSE
  )

  rds_path <- file.path(tmp_out, "ferx_save_results.RDS")
  expect_true(file.exists(rds_path))
  loaded <- readRDS(rds_path)
  expect_equal(nrow(loaded$rem), 1L)
  unlink(rds_path)
})


# =============================================================================
# Downstream functions work with ferx results
# =============================================================================

test_that("plot_removal_metrics works with ferx results", {
  skip_if_not_installed("ferx")

  model_path <- tempfile(fileext = ".ferx")
  write_test_ferx_model(model_path)
  dat <- make_ferx_test_data()
  on.exit(unlink(model_path))

  results <- remove_erroneous_obs(
    dat          = dat,
    mod          = model_path,
    run_id       = "ferx_plot",
    n            = 2L,
    engine       = "ferx",
    verbose      = FALSE,
    save_results = FALSE,
    stability_check = FALSE
  )

  # All metric types should produce ggplot objects without error
  p_pofv   <- plot_removal_metrics(results, metric = "pOFV", verbose = FALSE)
  p_thetas <- plot_removal_metrics(results, metric = "thetas", verbose = FALSE)
  p_rmse   <- plot_removal_metrics(results, metric = "rmse", verbose = FALSE)

  expect_s3_class(p_pofv, "gg")
  expect_s3_class(p_thetas, "gg")
  expect_s3_class(p_rmse, "gg")
})

test_that("check_model_stability works with ferx results", {
  skip_if_not_installed("ferx")

  model_path <- tempfile(fileext = ".ferx")
  write_test_ferx_model(model_path)
  dat <- make_ferx_test_data()
  on.exit(unlink(model_path))

  results <- remove_erroneous_obs(
    dat          = dat,
    mod          = model_path,
    run_id       = "ferx_stab",
    n            = 2L,
    engine       = "ferx",
    verbose      = FALSE,
    save_results = FALSE,
    stability_check = FALSE
  )

  stab <- check_model_stability(results, verbose = FALSE)
  expect_s3_class(stab, "irxclean_stability")
  expect_true("param_summary" %in% names(stab))
})
