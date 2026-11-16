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
                     "covariate_pair_outliers", "covariate_pair_scatter",
                     "missing_summary", "n_subjects", "n_observations",
                     "columns", "covariate_cols", "categorical_cols",
                     "covariate_pairs", "iqr_multiplier", "tad_bin_width",
                     "tad_label"),
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


# =============================================================================
# Covariate pair (bivariate) outlier tests
# =============================================================================

test_that("covariate pair returns NULL when no pairs specified", {
  dat <- minimal_pk_data(n_subjects = 20)
  qc  <- assess_data_quality(dat, verbose = FALSE)
  expect_null(qc$covariate_pair_outliers)
  expect_null(qc$covariate_pair_scatter)
  expect_null(qc$covariate_pairs)
})

test_that("covariate pair flags implausible WT~AGE combination", {
  # Build pediatric dataset with strong WT~AGE correlation
  set.seed(123)
  n <- 30
  ages <- runif(n, 1, 15)
  wts  <- 3 * ages + rnorm(n, 0, 2)  # tight WT ~ 3*AGE
  dat <- data.frame(
    ID   = rep(seq_len(n), each = 2),
    TIME = rep(c(0, 1), n),
    DV   = rep(c(0, 50), n),
    EVID = rep(c(1L, 0L), n),
    MDV  = rep(c(1L, 0L), n),
    AGE  = rep(round(ages, 1), each = 2),
    WT   = rep(round(wts, 1), each = 2),
    stringsAsFactors = FALSE
  )
  # Inject outlier: subject 1 is 10yo but only 5 kg (expected ~30 kg)
  dat$AGE[dat$ID == 1L] <- 10
  dat$WT[dat$ID == 1L]  <- 5
  qc <- assess_data_quality(
    dat,
    covariate_pairs = list(c("WT", "AGE")),
    iqr_multiplier  = 3,
    verbose = FALSE
  )
  expect_false(is.null(qc$covariate_pair_outliers))
  expect_true(1L %in% qc$covariate_pair_outliers$ID)
  expect_equal(qc$covariate_pair_outliers$Y_COL[1], "WT")
  expect_equal(qc$covariate_pair_outliers$X_COL[1], "AGE")
})

test_that("covariate pair returns NULL when no outliers in tight data", {
  set.seed(99)
  n <- 30
  dat <- data.frame(
    ID   = rep(seq_len(n), each = 2),
    TIME = rep(c(0, 1), n),
    DV   = rep(c(0, 50), n),
    EVID = rep(c(1L, 0L), n),
    MDV  = rep(c(1L, 0L), n),
    AGE  = rep(seq(20, 49), each = 2),
    WT   = rep(seq(50, 79), each = 2),
    stringsAsFactors = FALSE
  )
  qc <- assess_data_quality(
    dat,
    covariate_pairs = list(c("WT", "AGE")),
    iqr_multiplier  = 5,
    verbose = FALSE
  )
  expect_null(qc$covariate_pair_outliers)
  # Scatter should still exist
  expect_false(is.null(qc$covariate_pair_scatter))
})

test_that("covariate pair warns and skips when too few subjects", {
  dat <- minimal_pk_data(n_subjects = 5)
  expect_warning(
    qc <- assess_data_quality(
      dat,
      covariate_pairs   = list(c("WT", "AGE")),
      min_subjects_pair = 15L,
      verbose = FALSE
    ),
    "Only 5 subjects"
  )
  expect_null(qc$covariate_pair_outliers)
})

test_that("covariate pair warns on missing columns", {
  dat <- minimal_pk_data(n_subjects = 20)
  expect_warning(
    qc <- assess_data_quality(
      dat,
      covariate_pairs = list(c("NONEXIST", "AGE")),
      verbose = FALSE
    ),
    "not found"
  )
  expect_null(qc$covariate_pair_outliers)
})

test_that("covariate pair handles NA values gracefully", {
  dat <- minimal_pk_data(n_subjects = 25)
  # Set a few WT values to NA
  dat$WT[dat$ID %in% c(3L, 7L)] <- NA_real_
  qc <- assess_data_quality(
    dat,
    covariate_pairs = list(c("WT", "AGE")),
    verbose = FALSE
  )
  # Should not error; NA subjects excluded from scatter
  if (!is.null(qc$covariate_pair_scatter)) {
    expect_false(any(is.na(qc$covariate_pair_scatter$Y_VALUE)))
  }
})

test_that("covariate pair scatter data has correct structure", {
  dat <- minimal_pk_data(n_subjects = 20)
  qc <- assess_data_quality(
    dat,
    covariate_pairs = list(c("WT", "AGE"), c("AGE", "WT")),
    verbose = FALSE
  )
  scat <- qc$covariate_pair_scatter
  expect_false(is.null(scat))
  expect_true(all(c("ID", "Y_COL", "X_COL", "Y_VALUE", "X_VALUE",
                     "PREDICTED", "RESIDUAL", "OUTLIER") %in% names(scat)))
  # Two pairs means Y_COL should have 2 unique values
  expect_equal(length(unique(paste(scat$Y_COL, scat$X_COL))), 2L)
})

test_that("pair_method = 'reference' gives informative error", {
  dat <- minimal_pk_data(n_subjects = 20)
  expect_error(
    assess_data_quality(
      dat,
      covariate_pairs = list(c("WT", "AGE")),
      pair_method = "reference",
      verbose = FALSE
    ),
    "not yet implemented"
  )
})

test_that("covariate_pairs rejects non-list input", {
  dat <- minimal_pk_data(n_subjects = 20)
  expect_error(
    assess_data_quality(
      dat,
      covariate_pairs = c("WT", "AGE"),
      verbose = FALSE
    ),
    "must be a list"
  )
})

test_that("covariate_pairs rejects wrong-length elements", {
  dat <- minimal_pk_data(n_subjects = 20)
  expect_error(
    assess_data_quality(
      dat,
      covariate_pairs = list(c("WT", "AGE", "HT")),
      verbose = FALSE
    ),
    "length-2"
  )
})

test_that("plot with covariate pair panels generates ggplot", {
  set.seed(456)
  n <- 30
  ages <- runif(n, 1, 15)
  wts  <- 3 * ages + rnorm(n, 0, 2)
  dat <- data.frame(
    ID   = rep(seq_len(n), each = 2),
    TIME = rep(c(0, 1), n),
    DV   = rep(c(0, 50), n),
    EVID = rep(c(1L, 0L), n),
    MDV  = rep(c(1L, 0L), n),
    AGE  = rep(round(ages, 1), each = 2),
    WT   = rep(round(wts, 1), each = 2),
    stringsAsFactors = FALSE
  )
  dat$AGE[dat$ID == 1L] <- 10
  dat$WT[dat$ID == 1L]  <- 5
  qc <- assess_data_quality(
    dat,
    covariate_pairs = list(c("WT", "AGE")),
    iqr_multiplier  = 3,
    verbose = FALSE
  )
  pl <- plot(qc)
  if (inherits(pl, "gg")) {
    expect_s3_class(pl, "gg")
  } else {
    expect_true(is.list(pl))
    pair_names <- grep("^pair_", names(pl), value = TRUE)
    expect_true(length(pair_names) >= 1L)
    expect_s3_class(pl[[pair_names[1L]]], "gg")
  }
})

test_that("print includes covariate pair line when pairs specified", {
  dat <- minimal_pk_data(n_subjects = 20)
  qc <- assess_data_quality(
    dat,
    covariate_pairs = list(c("WT", "AGE")),
    verbose = FALSE
  )
  expect_output(print(qc), "Covariate pair outliers")
})

test_that("covariate_pairs = 'all' generates all pairwise combinations", {
  dat <- minimal_pk_data(n_subjects = 20)
  qc <- assess_data_quality(
    dat,
    covariate_cols  = c("AGE", "WT"),
    covariate_pairs = "all",
    verbose = FALSE
  )
  # combn(c("AGE","WT"), 2) gives 1 pair: c("AGE","WT")
  expect_equal(length(qc$covariate_pairs), 1L)
  expect_equal(qc$covariate_pairs[[1]], c("AGE", "WT"))
  # Scatter should exist for that pair
  expect_false(is.null(qc$covariate_pair_scatter))
  expect_true(all(qc$covariate_pair_scatter$Y_COL == "AGE"))
  expect_true(all(qc$covariate_pair_scatter$X_COL == "WT"))
})

test_that("covariate_pairs = 'all' with 3 covariates gives 3 pairs", {
  dat <- minimal_pk_data(n_subjects = 20)
  dat$HT <- rep(round(runif(20, 150, 190), 1), each = 6)
  qc <- assess_data_quality(
    dat,
    covariate_cols  = c("AGE", "WT", "HT"),
    covariate_pairs = "all",
    verbose = FALSE
  )
  # combn(3, 2) = 3 pairs

  expect_equal(length(qc$covariate_pairs), 3L)
})

test_that("covariate_pairs = 'all' errors with fewer than 2 covariate_cols", {
  dat <- minimal_pk_data(n_subjects = 20)
  expect_error(
    assess_data_quality(
      dat,
      covariate_cols  = "AGE",
      covariate_pairs = "all",
      verbose = FALSE
    ),
    "at least 2"
  )
})

test_that("covariate_pairs = 'all' errors with NULL covariate_cols", {
  dat <- minimal_pk_data(n_subjects = 20)
  expect_error(
    assess_data_quality(
      dat,
      covariate_pairs = "all",
      verbose = FALSE
    ),
    "at least 2"
  )
})


# =============================================================================
# Categorical x continuous covariate pair tests
# =============================================================================

test_that("categorical x continuous pair flags outlier within group", {
  set.seed(200)
  n <- 30
  dat <- data.frame(
    ID   = rep(seq_len(n), each = 2),
    TIME = rep(c(0, 1), n),
    DV   = rep(c(0, 50), n),
    EVID = rep(c(1L, 0L), n),
    MDV  = rep(c(1L, 0L), n),
    SEX  = rep(c(rep(0L, 15), rep(1L, 15)), each = 2),
    WT   = rep(c(rnorm(15, 80, 5), rnorm(15, 65, 5)), each = 2),
    stringsAsFactors = FALSE
  )
  # Inject outlier: subject 1 (SEX=0) with WT=5 (way below group mean ~80)
  dat$WT[dat$ID == 1L] <- 5
  qc <- assess_data_quality(
    dat,
    covariate_cols  = "WT",
    categorical_cols = "SEX",
    covariate_pairs = list(c("WT", "SEX")),
    iqr_multiplier  = 3,
    verbose = FALSE
  )
  expect_false(is.null(qc$covariate_pair_outliers))
  expect_true(1L %in% qc$covariate_pair_outliers$ID)
  expect_equal(qc$covariate_pair_outliers$PAIR_TYPE[1], "grouped")
})

test_that("categorical x continuous returns scatter with PAIR_TYPE = grouped", {
  dat <- minimal_pk_data(n_subjects = 20)
  qc <- assess_data_quality(
    dat,
    covariate_cols  = "WT",
    categorical_cols = "SEX",
    covariate_pairs = list(c("WT", "SEX")),
    verbose = FALSE
  )
  expect_false(is.null(qc$covariate_pair_scatter))
  expect_true(all(qc$covariate_pair_scatter$PAIR_TYPE == "grouped"))
})

test_that("categorical column order does not matter in pair", {
  dat <- minimal_pk_data(n_subjects = 20)
  # SEX first or second should both work
  qc1 <- assess_data_quality(
    dat,
    covariate_cols  = "WT",
    categorical_cols = "SEX",
    covariate_pairs = list(c("SEX", "WT")),
    verbose = FALSE
  )
  qc2 <- assess_data_quality(
    dat,
    covariate_cols  = "WT",
    categorical_cols = "SEX",
    covariate_pairs = list(c("WT", "SEX")),
    verbose = FALSE
  )
  # Both should produce scatter with continuous as Y_COL
  expect_equal(qc1$covariate_pair_scatter$Y_COL[1], "WT")
  expect_equal(qc2$covariate_pair_scatter$Y_COL[1], "WT")
})

test_that("covariate_pairs = 'all' includes cross-type pairs", {
  dat <- minimal_pk_data(n_subjects = 20)
  qc <- assess_data_quality(
    dat,
    covariate_cols  = c("AGE", "WT"),
    categorical_cols = "SEX",
    covariate_pairs = "all",
    verbose = FALSE
  )
  # 3 columns -> combn(3,2) = 3 pairs, minus cat x cat = 0
  # AGE~WT (cont x cont), AGE~SEX (cont x cat), WT~SEX (cont x cat)
  expect_equal(length(qc$covariate_pairs), 3L)
  types <- unique(qc$covariate_pair_scatter$PAIR_TYPE)
  expect_true("continuous" %in% types)
  expect_true("grouped" %in% types)
})

test_that("both-categorical pair is skipped with warning", {
  dat <- minimal_pk_data(n_subjects = 20)
  dat$RACE <- rep(sample(1:3, 20, replace = TRUE), each = 6)
  expect_warning(
    qc <- assess_data_quality(
      dat,
      categorical_cols = c("SEX", "RACE"),
      covariate_pairs  = list(c("SEX", "RACE")),
      verbose = FALSE
    ),
    "Both columns are categorical"
  )
  expect_null(qc$covariate_pair_outliers)
})

test_that("plot with grouped pair generates boxplot ggplot", {
  set.seed(201)
  n <- 30
  dat <- data.frame(
    ID   = rep(seq_len(n), each = 2),
    TIME = rep(c(0, 1), n),
    DV   = rep(c(0, 50), n),
    EVID = rep(c(1L, 0L), n),
    MDV  = rep(c(1L, 0L), n),
    SEX  = rep(c(rep(0L, 15), rep(1L, 15)), each = 2),
    WT   = rep(c(rnorm(15, 80, 5), rnorm(15, 65, 5)), each = 2),
    stringsAsFactors = FALSE
  )
  dat$WT[dat$ID == 1L] <- 5
  qc <- assess_data_quality(
    dat,
    covariate_cols  = "WT",
    categorical_cols = "SEX",
    covariate_pairs = list(c("WT", "SEX")),
    iqr_multiplier  = 3,
    verbose = FALSE
  )
  pl <- plot(qc)
  if (inherits(pl, "gg")) {
    expect_s3_class(pl, "gg")
  } else {
    expect_true(is.list(pl))
    pair_names <- grep("^pair_", names(pl), value = TRUE)
    expect_true(length(pair_names) >= 1L)
    expect_s3_class(pl[[pair_names[1L]]], "gg")
  }
})
