# Assess data quality before model-dependent cleaning

Performs a series of model-independent checks on a NONMEM-format
dataset:

## Usage

``` r
assess_data_quality(
  data,
  columns = list(ID = "ID", TIME = "TIME", DV = "DV", EVID = "EVID", MDV = "MDV"),
  covariate_cols = NULL,
  categorical_cols = NULL,
  iqr_multiplier = 3,
  tad_col = NULL,
  tad_bin_width = NULL,
  conc_increase_threshold = 1.5,
  conc_floor = 0,
  verbose = TRUE
)
```

## Arguments

- data:

  A data frame in NONMEM format.

- columns:

  Named list of column labels. Must include at minimum: \`ID\`,
  \`TIME\`, \`DV\`, \`EVID\`, \`MDV\`.

- covariate_cols:

  Character vector of continuous covariate column names to screen for
  outliers. \`NULL\` skips this check.

- categorical_cols:

  Character vector of categorical covariate column names. Currently used
  only for completeness in the missing-data summary.

- iqr_multiplier:

  Numeric. IQR multiplier for outlier bounds, applied to both covariate
  outlier detection and concentration bin outlier detection. Default
  \`3\`.

- tad_col:

  Character or `NULL`. Name of a pre-computed time-after-dose column.
  `NULL` (default) triggers auto-detection: looks for columns named
  `"TAD"`, `"TAFD"`, or `"TSFD"`; if none found, computes TAD from dose
  records (EVID == 1 / AMT \> 0); if no dose records exist, falls back
  to `TIME`.

- tad_bin_width:

  Numeric. Width of TAD bins used when checking for concentration
  outliers within each TAD window. `NULL` (default) auto-selects a round
  value from the TAD range. Bins with fewer than 4 observations are
  skipped.

- conc_increase_threshold:

  Numeric. Minimum fold-increase in concentration (without intervening
  dose) to flag as suspicious. Default \`1.5\`.

- conc_floor:

  Numeric. Minimum value for the preceding concentration when evaluating
  fold-increases. Pairs where the earlier DV is at or below this value
  are skipped (ratio is undefined / noise-dominated near the assay LOQ).
  Set to your assay LLOQ. Default \`0\` (no floor).

- verbose:

  Logical. Print a summary of findings? Default \`TRUE\`.

## Value

An object of class `irxclean_quality` (a named list) with elements:

- covariate_outliers:

  Data frame of flagged covariate values, or `NULL` if none detected /
  not requested.

- concentration_bin_outliers:

  Data frame of observations whose DV is extreme relative to other
  observations in the same TAD bin, or `NULL` if none detected. Columns:
  ID, TIME, TAD, TAD_BIN, DV, LOWER_BOUND, UPPER_BOUND, N_IN_BIN.

- conc_scatter_data:

  Data frame of all observation rows with columns `ID`, `TIME`, `TAD`,
  `DV`, and `OUTLIER` (logical). Used for the DV vs TAD scatter plot.

- concentration_flags:

  Data frame of suspicious concentration elevations without an
  intervening dose, or `NULL` if none detected.

- missing_summary:

  Data frame with columns `column`, `n_missing`, and `pct_missing`.

- n_subjects:

  Number of unique subjects.

- n_observations:

  Number of observation rows (EVID==0, MDV==0).

## Details

1.  \*\*Covariate outliers\*\* - flags subjects whose covariate values
    fall outside \\Q_1 - k \cdot IQR\\ or \\Q_3 + k \cdot IQR\\ (where
    \\k\\ = \`iqr_multiplier\`).

2.  \*\*Concentration TAD-bin outliers\*\* - observations whose DV falls
    outside IQR-based bounds within their time-after-dose (TAD) bin. TAD
    is read from a TAD column (auto-detected or specified via
    \`tad_col\`); if absent it is computed from dose records (EVID ==
    1).

3.  \*\*Concentration elevations without a preceding dose\*\* - flags
    consecutive observation pairs within a subject where the
    concentration rises by more than \`conc_increase_threshold\`-fold
    with no intervening dose event.

4.  \*\*Missing data summary\*\* - tabulates \`NA\` rates by column.

## See also

\[plot.irxclean_quality()\], \[remove_erroneous_obs()\]

## Examples

``` r
if (FALSE) { # \dontrun{
dat <- read.csv("my_data.csv")
qc  <- assess_data_quality(
  dat,
  covariate_cols   = c("AGE", "WT", "HT"),
  categorical_cols = "SEX"
)
print(qc)
plot(qc)
} # }
```
