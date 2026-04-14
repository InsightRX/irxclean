# Compare demographics of subjects with removed observations vs. whole dataset

After iterative removal, this function checks whether subjects whose
observations were removed are demographically representative of the full
dataset, and whether removals cluster at a particular position on the PK
curve (by time after dose). Systematic differences may indicate that the
removal algorithm is biased against a specific patient population or a
specific part of the concentration-time profile.

## Usage

``` r
plot_demographic_comparison(
  data,
  results,
  continuous_cols = NULL,
  categorical_cols = NULL,
  id_col = "ID",
  time_col = "TIME",
  evid_col = "EVID",
  tad_col = NULL,
  tad_bin_width = NULL,
  alpha = 0.05
)
```

## Arguments

- data:

  A data frame containing the original NONMEM dataset (one row per
  record).

- results:

  A results list from \[remove_erroneous_obs()\] or
  \[load_removal_results()\].

- continuous_cols:

  Character vector of continuous covariate column names (e.g.
  `c("AGE", "WT")`).

- categorical_cols:

  Character vector of categorical covariate column names (e.g. `"SEX"`).

- id_col:

  Character. Subject identifier column. Default `"ID"`.

- time_col:

  Character. Time column used to compute TAD when no explicit TAD column
  is found. Default `"TIME"`.

- evid_col:

  Character. EVID column used to distinguish observations (EVID == 0)
  from dose records (EVID == 1). Default `"EVID"`. If not present, all
  rows are treated as observations.

- tad_col:

  Character or `NULL`. Name of a pre-computed time-after- dose column.
  `NULL` (default) triggers auto-detection: the function looks for
  columns named `"TAD"`, `"TAFD"`, or `"TSFD"`; if none are found it
  computes TAD from the dose records in `data`; if dosing records cannot
  be identified it falls back to `time_col`.

- tad_bin_width:

  Numeric. Width of TAD bins for the curve-position analysis. `NULL`
  (default) auto-selects a round value from the TAD range.

- alpha:

  Numeric. Significance level. Default `0.05`.

## Value

An object of class `irxclean_demographics` (a named list):

- plots:

  Named list of `ggplot` objects, one per covariate.

- panels:

  List of `patchwork` grids, each containing up to 9 covariate plots.
  Significant covariates have red bold titles.

- tad_plot:

  A `ggplot` showing the proportion of removed observations per TAD bin,
  with Bonferroni-corrected per-bin binomial tests (elevated bins in
  red).

- stats:

  Data frame: `covariate`, `test`, `p_value`, `significant`.

- removed_ids:

  IDs with removed observations.

- n_removed_subjects:

  Number of subjects with any removal.

## Details

Statistical tests applied:

- Continuous covariates: two-sided Wilcoxon rank-sum test.

- Categorical covariates: chi-squared test of proportions.

- TAD bins: per-bin one-sided binomial test vs. the overall removal
  rate, Bonferroni-corrected.

## See also

\[plot_demographic_comparison()\], \[remove_erroneous_obs()\]

## Examples

``` r
if (FALSE) { # \dontrun{
original_data <- read.csv("my_data.csv")
results       <- load_removal_results("err_rem1")
demo <- plot_demographic_comparison(
  data             = original_data,
  results          = results,
  continuous_cols  = c("AGE", "WT", "HT"),
  categorical_cols = "SEX"
)
demo$panels[[1]]
demo$tad_plot
demo$stats
} # }
```
