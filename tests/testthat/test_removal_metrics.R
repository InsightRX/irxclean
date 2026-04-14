test_that("calculate_pofv_change returns expected columns", {
  res  <- minimal_results(n_iter = 5L)
  pofv <- calculate_pofv_change(res)
  expect_true(all(c("ITERATION", "pofv", "diff_pofv") %in% names(pofv)))
  expect_equal(nrow(pofv), 5L)
})

test_that("calculate_pofv_change first diff_pofv is zero", {
  res  <- minimal_results(n_iter = 4L)
  pofv <- calculate_pofv_change(res)
  expect_equal(pofv$diff_pofv[1], 0)
})

test_that("calculate_pofv_change errors when rem IDs not in phi", {
  res <- minimal_results()
  res$rem$ID <- 99999L  # ID not in phi
  expect_error(calculate_pofv_change(res), "not found in phi")
})

test_that("plot_removal_metrics returns a ggplot for each metric", {
  res <- minimal_results(n_iter = 5L)
  for (metric in c("pOFV", "thetas", "omegas", "rmse")) {
    p <- plot_removal_metrics(res, metric = metric, verbose = FALSE)
    expect_s3_class(p, "gg")
  }
})

test_that("plot_removal_metrics errors on invalid metric", {
  res <- minimal_results()
  expect_error(plot_removal_metrics(res, metric = "invalid"))
})

test_that(".validate_results errors on missing components", {
  res_bad <- list(par = data.frame(), rem = data.frame())
  expect_error(irxclean:::.validate_results(res_bad), "missing components")
})

test_that(".calculate_nrmse is correct", {
  obs  <- c(1, 2, 3, 4, 5)
  pred <- c(1, 2, 3, 4, 5)
  expect_equal(irxclean:::.calculate_nrmse(obs, pred), 0)
  # rmse = 1, mean(obs) = 2, nrmse = 0.5
  expect_equal(irxclean:::.calculate_nrmse(c(2, 2), c(3, 3)), 0.5)
})

test_that("load_removal_results works with RDS fixture", {
  rds_path <- system.file("testdata", "rem_test_results.RDS", package = "irxclean")
  skip_if(!file.exists(rds_path), "Test fixture not available")
  res <- readRDS(rds_path)
  expect_true(all(c("par", "rem", "phi", "rmse") %in% names(res)))
})

test_that("load_removal_results errors when file missing", {
  withr::with_dir(tempdir(), {
    expect_error(load_removal_results("nonexistent_run"), "not found")
  })
})
