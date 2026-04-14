# =============================================================================
# Tests for apply_exclusion_criteria()
# =============================================================================

# Helper: minimal NONMEM dataset with dose + observations and AMT/RATE
make_excl_data <- function(n_subjects = 10, n_obs = 4) {
  set.seed(99)
  n_per_subj <- n_obs + 1L
  ids  <- rep(seq_len(n_subjects), each = n_per_subj)
  evid <- rep(c(1L, rep(0L, n_obs)), times = n_subjects)
  mdv  <- ifelse(evid == 1L, 1L, 0L)
  data.frame(
    ID   = ids,
    TIME = rep(c(0, seq_len(n_obs)), times = n_subjects),
    DV   = ifelse(evid == 0L, stats::rexp(length(ids)) * 100, 0),
    EVID = evid,
    MDV  = mdv,
    AMT  = ifelse(evid == 1L, 100, 0),
    RATE = 0,
    stringsAsFactors = FALSE
  )
}

# =============================================================================
# Return structure
# =============================================================================

test_that("apply_exclusion_criteria returns irxclean_exclusions with expected names", {
  dat  <- make_excl_data()
  excl <- apply_exclusion_criteria(dat, verbose = FALSE)
  expect_s3_class(excl, "irxclean_exclusions")
  expect_named(
    excl,
    c("data_flagged", "data_clean", "exclusion_summary", "criteria_applied"),
    ignore.order = TRUE
  )
})

test_that("data_flagged has same nrow as input and extra excl_ columns", {
  dat  <- make_excl_data()
  excl <- apply_exclusion_criteria(dat, verbose = FALSE)
  expect_equal(nrow(excl$data_flagged), nrow(dat))
  excl_cols <- grep("^excl_", names(excl$data_flagged), value = TRUE)
  expect_gt(length(excl_cols), 0L)
  expect_true(all(vapply(excl$data_flagged[excl_cols], is.logical, logical(1L))))
})

test_that("data_clean has fewer rows than data_flagged when exclusions apply", {
  dat <- make_excl_data()
  # Subject 1 gets an outlier dose
  dat$AMT[dat$ID == 1L & dat$EVID == 1L] <- 99999
  excl <- apply_exclusion_criteria(
    dat, check_dose_outlier = TRUE, dose_iqr_multiplier = 3, verbose = FALSE
  )
  expect_lt(nrow(excl$data_clean), nrow(dat))
})

test_that("exclusion_summary has required columns", {
  dat  <- make_excl_data()
  excl <- apply_exclusion_criteria(dat, verbose = FALSE)
  expect_true(all(c("criterion", "n_records", "n_subjects", "pct_records") %in%
                    names(excl$exclusion_summary)))
})

test_that("criteria_applied matches columns added to data_flagged", {
  dat  <- make_excl_data()
  excl <- apply_exclusion_criteria(dat, verbose = FALSE)
  excl_cols <- grep("^excl_", names(excl$data_flagged), value = TRUE)
  expect_setequal(excl$criteria_applied, excl_cols)
})

# =============================================================================
# excl_no_doses
# =============================================================================

test_that("excl_no_doses flags subjects with no EVID==1 rows", {
  dat <- make_excl_data(n_subjects = 5)
  # Remove dose for subject 3
  dat <- dat[!(dat$ID == 3L & dat$EVID == 1L), ]
  excl <- apply_exclusion_criteria(dat, verbose = FALSE)
  flagged_ids <- unique(excl$data_flagged$ID[excl$data_flagged$excl_no_doses])
  expect_true(3L %in% flagged_ids)
  expect_false(1L %in% flagged_ids)
})

test_that("excl_no_doses is all FALSE when all subjects have doses", {
  dat  <- make_excl_data()
  excl <- apply_exclusion_criteria(dat, check_no_doses = TRUE, verbose = FALSE)
  expect_false(any(excl$data_flagged$excl_no_doses))
})

# =============================================================================
# excl_no_tdm
# =============================================================================

test_that("excl_no_tdm flags subjects with no observation rows", {
  dat <- make_excl_data(n_subjects = 5)
  # Remove all obs rows for subject 2
  dat <- dat[!(dat$ID == 2L & dat$EVID == 0L), ]
  excl <- apply_exclusion_criteria(dat, verbose = FALSE)
  flagged_ids <- unique(excl$data_flagged$ID[excl$data_flagged$excl_no_tdm])
  expect_true(2L %in% flagged_ids)
  expect_false(1L %in% flagged_ids)
})

# =============================================================================
# excl_during_infusion
# =============================================================================

test_that("excl_during_infusion flags obs during ongoing infusion", {
  # Single subject: infusion dose at t=0 (AMT=100, RATE=10 -> ends at t=10)
  # Observations at t=5 (during) and t=15 (after)
  dat <- data.frame(
    ID   = c(1L, 1L, 1L),
    TIME = c(0, 5, 15),
    DV   = c(0, 80, 20),
    EVID = c(1L, 0L, 0L),
    MDV  = c(1L, 0L, 0L),
    AMT  = c(100, 0, 0),
    RATE = c(10, 0, 0),
    stringsAsFactors = FALSE
  )
  excl <- apply_exclusion_criteria(
    dat,
    check_no_doses = FALSE, check_no_tdm = FALSE,
    check_conc_increase = FALSE, check_dose_outlier = FALSE,
    verbose = FALSE
  )
  expect_true(excl$data_flagged$excl_during_infusion[2L])   # t=5 flagged
  expect_false(excl$data_flagged$excl_during_infusion[3L])  # t=15 not flagged
})

test_that("excl_during_infusion is all FALSE when RATE=0 (bolus)", {
  dat  <- make_excl_data()  # all RATE = 0
  excl <- apply_exclusion_criteria(
    dat, check_during_infusion = TRUE, verbose = FALSE
  )
  expect_false(any(excl$data_flagged$excl_during_infusion))
})

# =============================================================================
# excl_conc_increase
# =============================================================================

test_that("excl_conc_increase flags concentration rise without dose", {
  dat <- data.frame(
    ID   = c(1L, 1L, 1L),
    TIME = c(0, 1, 2),
    DV   = c(0, 10, 100),
    EVID = c(1L, 0L, 0L),
    MDV  = c(1L, 0L, 0L),
    AMT  = c(100, 0, 0),
    RATE = c(0, 0, 0),
    stringsAsFactors = FALSE
  )
  excl <- apply_exclusion_criteria(
    dat,
    check_no_doses = FALSE, check_no_tdm = FALSE,
    check_during_infusion = FALSE, check_dose_outlier = FALSE,
    conc_increase_threshold = 1.5,
    verbose = FALSE
  )
  # DV goes from 10 to 100 (10x) with no intervening dose -> flagged
  flagged_times <- dat$TIME[excl$data_flagged$excl_conc_increase]
  expect_true(2 %in% flagged_times)
})

# =============================================================================
# excl_dose_outlier
# =============================================================================

test_that("excl_dose_outlier flags extreme AMT values", {
  dat <- make_excl_data(n_subjects = 20)
  # Inject an extreme dose for subject 1
  dat$AMT[dat$ID == 1L & dat$EVID == 1L] <- 100000
  excl <- apply_exclusion_criteria(
    dat, check_dose_outlier = TRUE, dose_iqr_multiplier = 3, verbose = FALSE
  )
  flagged_ids <- unique(dat$ID[excl$data_flagged$excl_dose_outlier])
  expect_true(1L %in% flagged_ids)
})

test_that("excl_dose_outlier is all FALSE when doses are uniform", {
  dat  <- make_excl_data()  # all AMT = 100
  excl <- apply_exclusion_criteria(
    dat, check_dose_outlier = TRUE, verbose = FALSE
  )
  expect_false(any(excl$data_flagged$excl_dose_outlier))
})

# =============================================================================
# Custom criteria
# =============================================================================

test_that("custom_criteria adds excl_<name> column", {
  dat  <- make_excl_data(n_subjects = 5)
  flag <- dat$ID == 2L
  excl <- apply_exclusion_criteria(
    dat,
    custom_criteria = list(my_flag = flag),
    verbose = FALSE
  )
  expect_true("excl_my_flag" %in% names(excl$data_flagged))
  expect_equal(excl$data_flagged$excl_my_flag, flag)
})

test_that("custom_criteria wrong length throws error", {
  dat <- make_excl_data()
  expect_error(
    apply_exclusion_criteria(
      dat,
      custom_criteria = list(bad = rep(TRUE, 3L)),
      verbose = FALSE
    ),
    "length"
  )
})

# =============================================================================
# Toggle checks off
# =============================================================================

test_that("all checks off produces no excl_ columns in criteria_applied", {
  dat  <- make_excl_data()
  excl <- apply_exclusion_criteria(
    dat,
    check_no_doses        = FALSE,
    check_no_tdm          = FALSE,
    check_during_infusion = FALSE,
    check_conc_increase   = FALSE,
    check_dose_outlier    = FALSE,
    verbose               = FALSE
  )
  expect_length(excl$criteria_applied, 0L)
  expect_equal(nrow(excl$data_clean), nrow(dat))
})

# =============================================================================
# Input validation
# =============================================================================

test_that("missing required column raises error", {
  dat <- data.frame(ID = 1, X = 1)
  expect_error(
    apply_exclusion_criteria(dat, verbose = FALSE),
    "not found"
  )
})

test_that("missing columns list key raises error", {
  dat <- make_excl_data()
  expect_error(
    apply_exclusion_criteria(dat, columns = list(ID = "ID"), verbose = FALSE),
    "missing keys"
  )
})

# =============================================================================
# print method
# =============================================================================

test_that("print.irxclean_exclusions prints without error", {
  dat  <- make_excl_data()
  excl <- apply_exclusion_criteria(dat, verbose = FALSE)
  expect_output(print(excl), "irxclean exclusion criteria summary")
})

# =============================================================================
# pct_records reflects obs rows only
# =============================================================================

test_that("pct_records in summary is based on observation rows not total rows", {
  # 1 subject with 1 dose + 1 obs; flag everything
  dat <- data.frame(
    ID = 1L, TIME = c(0, 1), DV = c(0, 50),
    EVID = c(1L, 0L), MDV = c(1L, 0L),
    AMT = c(100, 0), RATE = c(0, 0),
    stringsAsFactors = FALSE
  )
  # Custom flag for only the obs row
  excl <- apply_exclusion_criteria(
    dat,
    check_no_doses        = FALSE,
    check_no_tdm          = FALSE,
    check_during_infusion = FALSE,
    check_conc_increase   = FALSE,
    check_dose_outlier    = FALSE,
    custom_criteria       = list(test = c(FALSE, TRUE)),
    verbose               = FALSE
  )
  summary_row <- excl$exclusion_summary[excl$exclusion_summary$criterion == "excl_test", ]
  # 1 obs row flagged / 1 total obs row = 100%
  expect_equal(summary_row$pct_records, 100)
})
