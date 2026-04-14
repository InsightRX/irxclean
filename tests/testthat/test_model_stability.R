test_that("check_model_stability returns expected structure", {
  res  <- minimal_results(n_iter = 5L)
  stab <- check_model_stability(res, verbose = FALSE)
  expect_s3_class(stab, "irxclean_stability")
  expect_named(
    stab,
    c("param_summary", "parameter_drift", "unstable_params",
      "flagged_iterations", "n_flagged_params", "rmse_trend",
      "threshold_pct"),
    ignore.order = TRUE
  )
})

test_that("check_model_stability flags bouncing parameter as unstable", {
  res <- minimal_results(n_iter = 8L)
  # Alternate THETA1 between 10 and 100: RSE ~88%, 6 reversals -> unstable
  res$par$THETA1 <- rep(c(10, 100), length.out = 8L)
  stab <- check_model_stability(res, threshold_pct = 20, verbose = FALSE)
  expect_true("THETA1" %in% stab$unstable_params)
  expect_true(stab$n_flagged_params >= 1L)
  expect_true(length(stab$flagged_iterations) >= 1L)
})

test_that("check_model_stability reports stable when no drift", {
  # Constant params are treated as fixed and excluded; no unstable params
  res <- minimal_results(n_iter = 4L)
  res$par$THETA1 <- 10
  res$par$THETA2 <- 50
  res$par$OMEGA1 <- 0.05
  stab <- check_model_stability(res, threshold_pct = 20, verbose = FALSE)
  expect_equal(stab$n_flagged_params, 0L)
  expect_equal(length(stab$flagged_iterations), 0L)
})

test_that("check_model_stability detects increasing RMSE", {
  res <- minimal_results(n_iter = 6L)
  res$rmse <- data.frame(
    rmse      = c(10, 15, 20, 25, 30, 35),
    ITERATION = 1:6
  )
  stab <- check_model_stability(res, verbose = FALSE)
  expect_equal(stab$rmse_trend, "increasing")
})

test_that("check_model_stability detects decreasing RMSE", {
  res <- minimal_results(n_iter = 6L)
  res$rmse <- data.frame(
    rmse      = c(30, 25, 20, 15, 10, 5),
    ITERATION = 1:6
  )
  stab <- check_model_stability(res, verbose = FALSE)
  expect_equal(stab$rmse_trend, "decreasing")
})

test_that("print.irxclean_stability outputs text", {
  res  <- minimal_results()
  stab <- check_model_stability(res, verbose = FALSE)
  expect_output(print(stab), "irxclean model stability")
})

test_that("stability_flag_table returns data frame of unstable params", {
  res <- minimal_results(n_iter = 8L)
  res$par$THETA1 <- rep(c(10, 100), length.out = 8L)
  stab <- check_model_stability(res, threshold_pct = 20, verbose = FALSE)
  tbl  <- stability_flag_table(stab)
  expect_true(is.data.frame(tbl))
  expect_true(nrow(tbl) >= 1L)
  expect_true(all(tbl$STATUS == "unstable"))
})

test_that("stability_flag_table errors on wrong class", {
  expect_error(stability_flag_table(list()), "irxclean_stability")
})

test_that("plot.irxclean_stability returns a list with named plots", {
  res <- minimal_results(n_iter = 5L)
  res$par$THETA1[res$par$ITERATION == 3L] <- 999
  stab <- check_model_stability(res, threshold_pct = 20, verbose = FALSE)
  pl   <- plot(stab)
  expect_true(is.list(pl))
  expect_true("parameter_drift" %in% names(pl))
  expect_true("param_rse" %in% names(pl))
})
