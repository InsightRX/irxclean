test_that("remove_observations_from_data removes correct rows", {
  dat <- minimal_pk_data(n_subjects = 5, n_obs = 5)
  # Manually craft a results object pointing to an existing obs
  obs_row <- dat[dat$EVID == 0, ][1, ]
  fake_rem <- data.frame(
    ID        = obs_row$ID,
    TIME      = obs_row$TIME,
    DV        = obs_row$DV,
    CWRES     = 4.2,
    ITERATION = 1L
  )
  fake_results <- minimal_results(n_iter = 1L)
  fake_results$rem <- fake_rem

  cleaned <- remove_observations_from_data(
    data    = dat,
    results = fake_results,
    n       = 1L,
    columns = list(ID = "ID", TIME = "TIME", DV = "DV")
  )
  expect_equal(nrow(cleaned), nrow(dat) - 1L)
})

test_that("remove_observations_from_data n=NULL removes all rem rows", {
  res <- minimal_results(n_iter = 3L)
  dat <- minimal_pk_data(n_subjects = 20, n_obs = 10)
  # Inject some valid rows into rem
  obs_rows <- dat[dat$EVID == 0, ][1:3, ]
  res$rem <- data.frame(
    ID = obs_rows$ID, TIME = obs_rows$TIME, DV = obs_rows$DV,
    CWRES = c(4, 3, 2), ITERATION = 1:3
  )
  cleaned <- remove_observations_from_data(dat, res, n = NULL)
  expect_equal(nrow(cleaned), nrow(dat) - 3L)
})

test_that("remove_observations_from_data warns on row count mismatch", {
  res <- minimal_results(n_iter = 1L)
  dat <- minimal_pk_data(n_subjects = 5)
  # rem points to a non-existent row
  res$rem <- data.frame(ID = 999, TIME = 999, DV = 999, CWRES = 5, ITERATION = 1L)
  expect_warning(
    remove_observations_from_data(dat, res),
    "Removed 0 rows but requested 1"
  )
})

test_that("remove_observations_from_data errors on bad n", {
  res <- minimal_results()
  dat <- minimal_pk_data()
  expect_error(remove_observations_from_data(dat, res, n = 0L), "positive integer")
  expect_error(remove_observations_from_data(dat, res, n = -1L), "positive integer")
})
