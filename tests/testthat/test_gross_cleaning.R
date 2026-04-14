test_that("assess_data_quality returns expected structure", {
  dat <- minimal_pk_data(n_subjects = 15)
  qc  <- assess_data_quality(
    dat,
    covariate_cols = c("AGE", "WT"),
    categorical_cols = "SEX",
    verbose = FALSE
  )
  expect_s3_class(qc, "irxclean_quality")
  expect_true(is.list(qc))
  expect_named(qc, c("covariate_outliers", "concentration_bin_outliers",
                     "conc_scatter_data", "concentration_flags",
                     "missing_summary", "n_subjects", "n_observations",
                     "columns", "covariate_cols", "categorical_cols",
                     "iqr_multiplier", "tad_bin_width", "tad_label"),
               ignore.order = TRUE)
})

test_that("n_subjects and n_observations are correct", {
  dat <- minimal_pk_data(n_subjects = 8)
  qc  <- assess_data_quality(dat, verbose = FALSE)
  expect_equal(qc$n_subjects, 8L)
  # 8 subjects × 5 observation rows
  expect_equal(qc$n_observations, 8L * 5L)
})

test_that("missing_summary has correct columns and no NAs in n_missing", {
  dat <- minimal_pk_data()
  qc  <- assess_data_quality(dat, verbose = FALSE)
  ms  <- qc$missing_summary
  expect_true(all(c("column", "n_missing", "pct_missing") %in% names(ms)))
  expect_false(any(is.na(ms$n_missing)))
})

test_that("covariate outlier detection flags extreme values", {
  dat <- minimal_pk_data(n_subjects = 20)
  # Inject an extreme value for one subject's AGE
  dat$AGE[dat$ID == 1L] <- 999
  qc <- assess_data_quality(
    dat,
    covariate_cols = "AGE",
    iqr_multiplier = 3,
    verbose = FALSE
  )
  expect_false(is.null(qc$covariate_outliers))
  expect_true(1L %in% qc$covariate_outliers$ID)
})

test_that("covariate outlier detection returns NULL when none present", {
  # Homogeneous dataset – no outliers
  dat <- data.frame(
    ID   = 1:20,
    TIME = 0,
    DV   = 0,
    EVID = 1L,
    MDV  = 1L,
    AGE  = rep(40, 20),
    stringsAsFactors = FALSE
  )
  qc <- assess_data_quality(dat, covariate_cols = "AGE", verbose = FALSE)
  expect_null(qc$covariate_outliers)
})

test_that("concentration increase detection flags rise without dose", {
  # Single subject: dose at 0, obs at 1 (low), obs at 2 (very high without dose)
  dat <- data.frame(
    ID   = c(1, 1, 1),
    TIME = c(0, 1, 2),
    DV   = c(0, 10, 100),  # 10× increase, no dose between
    EVID = c(1L, 0L, 0L),
    MDV  = c(1L, 0L, 0L),
    stringsAsFactors = FALSE
  )
  qc <- assess_data_quality(
    dat,
    conc_increase_threshold = 1.5,
    verbose = FALSE
  )
  expect_false(is.null(qc$concentration_flags))
  expect_equal(nrow(qc$concentration_flags), 1L)
  expect_true(qc$concentration_flags$RATIO > 1.5)
})

test_that("concentration increase detection does not flag rise after dose", {
  dat <- data.frame(
    ID   = c(1, 1, 1, 1),
    TIME = c(0, 1, 2, 3),
    DV   = c(0, 10, 0, 100),
    EVID = c(1L, 0L, 1L, 0L),  # second dose at TIME=2
    MDV  = c(1L, 0L, 1L, 0L),
    stringsAsFactors = FALSE
  )
  qc <- assess_data_quality(dat, conc_increase_threshold = 1.5, verbose = FALSE)
  expect_null(qc$concentration_flags)
})

test_that("print.irxclean_quality works without error", {
  dat <- minimal_pk_data()
  qc  <- assess_data_quality(dat, verbose = FALSE)
  expect_output(print(qc), "irxclean gross data quality assessment")
})

test_that("plot.irxclean_quality returns ggplot for missing data", {
  dat <- minimal_pk_data()
  dat$AGE[1] <- NA_real_
  qc <- assess_data_quality(dat, verbose = FALSE)
  pl <- plot(qc)
  # Should be a list with at least missing plot
  if (inherits(pl, "gg")) {
    expect_s3_class(pl, "gg")
  } else {
    expect_true(is.list(pl))
    expect_true(any(vapply(pl, inherits, logical(1), "gg")))
  }
})

test_that("assess_data_quality errors on missing required columns", {
  dat <- data.frame(ID = 1, X = 1)
  expect_error(
    assess_data_quality(dat, verbose = FALSE),
    "not found"
  )
})

test_that("unknown covariate columns produce a warning, not an error", {
  dat <- minimal_pk_data()
  expect_warning(
    assess_data_quality(dat, covariate_cols = c("AGE", "NONEXISTENT"), verbose = FALSE),
    "not found"
  )
})
