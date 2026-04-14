# Apply structured exclusion criteria to a NONMEM-format dataset

Applies a battery of model-independent exclusion checks and returns a
flagged dataset together with a per-criterion summary. This function is
designed for the *initial* data-cleaning step run once after data
wrangling, before a stable base model is available.

## Usage

``` r
apply_exclusion_criteria(
  data,
  columns = list(ID = "ID", TIME = "TIME", DV = "DV", EVID = "EVID", MDV = "MDV", AMT =
    "AMT", RATE = "RATE"),
  check_no_doses = TRUE,
  check_no_tdm = TRUE,
  check_during_infusion = TRUE,
  check_conc_increase = TRUE,
  conc_increase_threshold = 1.5,
  conc_floor = 0,
  check_dose_outlier = TRUE,
  dose_iqr_multiplier = 3,
  custom_criteria = NULL,
  verbose = TRUE
)
```

## Arguments

- data:

  A data frame in NONMEM format.

- columns:

  Named list mapping logical column roles to actual column names.
  Required keys: `ID`, `TIME`, `DV`, `EVID`, `MDV`. Optional: `AMT`,
  `RATE` (needed for dose-outlier and infusion checks).

- check_no_doses:

  Logical. Flag all rows for subjects who have no dose records (EVID ==
  1)? Default `TRUE`.

- check_no_tdm:

  Logical. Flag all rows for subjects who have no TDM observations (EVID
  == 0, MDV == 0)? Default `TRUE`.

- check_during_infusion:

  Logical. Flag observation rows collected during an ongoing infusion?
  Requires `AMT` and `RATE` columns. Default `TRUE`.

- check_conc_increase:

  Logical. Flag observations where concentration rises by more than
  `conc_increase_threshold`-fold without an intervening dose? Default
  `TRUE`.

- conc_increase_threshold:

  Numeric. Fold-increase threshold for flagging concentration
  elevations. Default `1.5`.

- conc_floor:

  Numeric. Minimum preceding concentration for evaluating fold-increases
  (pairs where DV_prev \<= `conc_floor` are skipped). Default `0`.

- check_dose_outlier:

  Logical. Flag dose records whose AMT is an IQR-based outlier across
  all dose records? Requires `AMT` column. Default `TRUE`.

- dose_iqr_multiplier:

  Numeric. IQR multiplier for dose-outlier detection. Default `3`.

- custom_criteria:

  Named list of row-level logical vectors (one per criterion, length
  equal to `nrow(data)`). Each vector is added as an `excl_<name>`
  column. `NULL` skips custom criteria.

- verbose:

  Logical. Print a summary to the console? Default `TRUE`.

## Value

An object of class `irxclean_exclusions` (a named list):

- data_flagged:

  The original data frame with `excl_*` logical columns appended.

- data_clean:

  Rows from `data` where no `excl_*` column is `TRUE`.

- exclusion_summary:

  Data frame with columns `criterion`, `n_records`, `n_subjects`, and
  `pct_records` (percentage of total observation rows affected).

- criteria_applied:

  Character vector of `excl_*` column names added.

## Details

Each criterion adds a logical `excl_*` column to the data. Rows where
*any* exclusion column is `TRUE` are removed to produce `data_clean`.

## See also

\[assess_data_quality()\], \[remove_erroneous_obs()\]

## Examples

``` r
if (FALSE) { # \dontrun{
dat <- read.csv("my_data.csv")
excl <- apply_exclusion_criteria(
  dat,
  check_dose_outlier = TRUE,
  dose_iqr_multiplier = 3
)
print(excl)
head(excl$data_clean)
} # }
```
