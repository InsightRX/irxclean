# =============================================================================
# Tests for ferx utility functions (R/utils_ferx.R)
# =============================================================================

# =============================================================================
# ferx_read_pars
# =============================================================================

test_that("ferx_read_pars converts fit result to flat named list", {
  fit_result <- list(
    theta = c(10.5, 50.2, 3.1),
    omega = matrix(c(0.09, 0, 0, 0.04), nrow = 2),
    sigma = matrix(0.01, nrow = 1),
    ofv   = -1234.5
  )

  pars <- ferx_read_pars(fit_result)

  expect_type(pars, "list")
  expect_equal(pars$THETA1, 10.5)
  expect_equal(pars$THETA2, 50.2)
  expect_equal(pars$THETA3, 3.1)
  expect_equal(pars$OMEGA.1.1, 0.09)
  expect_equal(pars$OMEGA.2.1, 0)
  expect_equal(pars$OMEGA.2.2, 0.04)
  expect_equal(pars$SIGMA.1.1, 0.01)
  expect_equal(pars$OBJ, -1234.5)
})

test_that("ferx_read_pars handles vector omega/sigma", {
  fit_result <- list(
    theta = c(5.0),
    omega = c(0.1, 0.2),
    sigma = c(0.05),
    ofv   = -500
  )

  pars <- ferx_read_pars(fit_result)

  expect_equal(pars$THETA1, 5.0)
  expect_equal(pars$OMEGA.1.1, 0.1)
  expect_equal(pars$OMEGA.2.2, 0.2)
  expect_equal(pars$SIGMA.1.1, 0.05)
})

# =============================================================================
# ferx_parse_param_labels
# =============================================================================

test_that("ferx_parse_param_labels extracts labels from a .ferx file", {
  tmp <- tempfile(fileext = ".ferx")
  writeLines(c(
    "[parameters]",
    "tvCL = 10   ; Clearance",
    "tvV  = 50   ; Volume",
    "omega_CL = 0.09 ; IIV CL",
    "sigma_prop = 0.01 ; Proportional error",
    "",
    "[model]",
    "# model code here"
  ), tmp)
  on.exit(unlink(tmp))

  labels <- ferx_parse_param_labels(tmp)

  expect_type(labels, "list")
  expect_named(labels, c("thetas", "omegas", "sigmas"))
  expect_equal(labels$thetas[["THETA1"]], "Clearance")
  expect_equal(labels$thetas[["THETA2"]], "Volume")
  expect_equal(labels$omegas[["OMEGA.1.1"]], "IIV CL")
  expect_equal(labels$sigmas[["SIGMA.1.1"]], "Proportional error")
})

test_that("ferx_parse_param_labels returns empty when no [parameters] block", {
  tmp <- tempfile(fileext = ".ferx")
  writeLines(c(
    "[model]",
    "# model code only"
  ), tmp)
  on.exit(unlink(tmp))

  labels <- ferx_parse_param_labels(tmp)

  expect_equal(labels$thetas, character(0))
  expect_equal(labels$omegas, character(0))
  expect_equal(labels$sigmas, character(0))
})

test_that("ferx_parse_param_labels returns empty for nonexistent file", {
  labels <- ferx_parse_param_labels("/nonexistent/model.ferx")
  expect_equal(labels$thetas, character(0))
  expect_equal(labels$omegas, character(0))
  expect_equal(labels$sigmas, character(0))
})

# =============================================================================
# ferx_write_data
# =============================================================================

test_that("ferx_write_data writes a CSV file", {
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp))
  df <- data.frame(ID = 1:3, TIME = c(0, 1, 2), DV = c(0, 10, 5))
  ferx_write_data(df, tmp)
  expect_true(file.exists(tmp))
  result <- utils::read.csv(tmp)
  expect_equal(nrow(result), 3L)
  expect_equal(names(result), c("ID", "TIME", "DV"))
})

# =============================================================================
# ferx_run_iteration (requires ferx package)
# =============================================================================

test_that("ferx_run_iteration errors when ferx is not installed", {
  skip_if(requireNamespace("ferx", quietly = TRUE),
          "ferx is installed; skipping missing-package test")
  expect_error(
    ferx_run_iteration(
      model_path = "test.ferx",
      data       = data.frame(ID = 1),
      run_id     = "test",
      iteration  = 1L
    ),
    "ferx"
  )
})

# =============================================================================
# remove_erroneous_obs engine validation
# =============================================================================

test_that("remove_erroneous_obs rejects invalid engine", {
  expect_error(
    remove_erroneous_obs(
      dat    = data.frame(ID = 1),
      mod    = "test",
      run_id = "test",
      n      = 1,
      engine = "invalid"
    ),
    "arg"
  )
})

test_that("remove_erroneous_obs with engine='ferx' errors when ferx not installed", {
  skip_if(requireNamespace("ferx", quietly = TRUE),
          "ferx is installed; skipping missing-package test")
  expect_error(
    remove_erroneous_obs(
      dat    = data.frame(ID = 1),
      mod    = "test.ferx",
      run_id = "test",
      n      = 1,
      engine = "ferx"
    ),
    "ferx"
  )
})
