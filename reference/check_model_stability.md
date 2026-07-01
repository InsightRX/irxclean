# Check model stability across removal iterations

Evaluates whether model parameters reached stable estimates as
observations are iteratively removed. Rather than flagging every
iteration where a parameter drifted from baseline (which conflates
expected convergence to a new value with genuine instability), this
function classifies each parameter as:

## Usage

``` r
check_model_stability(
  results,
  threshold_pct = 20,
  late_window_fraction = 1/3,
  verbose = TRUE
)
```

## Arguments

- results:

  A results list from
  [`remove_erroneous_obs()`](https://insightrx.github.io/irxclean/reference/remove_erroneous_obs.md)
  or
  [`load_removal_results()`](https://insightrx.github.io/irxclean/reference/load_removal_results.md).

- threshold_pct:

  Numeric. RSE % threshold above which a parameter is considered to have
  high variability. Also used (as `threshold_pct / 4`) to assess
  late-iteration stability. Default 20 (i.e. 20%).

- late_window_fraction:

  Numeric. Fraction of iterations (from the end) used to assess
  late-stage CV stability. Default `1/3` (final third, minimum 2
  iterations).

- verbose:

  Logical. Print a stability summary to the console? Default TRUE.

## Value

An object of class `irxclean_stability` (a named list):

- param_summary:

  Data frame with one row per parameter containing `PARAMETER`, `MEAN`,
  `FIRST`, `LAST`, `RSE_PCT` (SD / \|mean\| x 100), `N_REVERSALS` (trend
  direction changes), `CV_LATE_PCT` (CV of the final third of
  iterations), `DRIFT_PCT` (\|first - last\| / \|first\| x 100), and
  `STATUS` (stable / drifted / unstable).

- parameter_drift:

  Data frame of per-iteration values with `PCT_CHANGE` from baseline and
  `STATUS`, used for the trajectory plot.

- unstable_params:

  Character vector of parameters classified as unstable.

- flagged_iterations:

  Integer vector of iterations associated with unstable parameters.

- n_flagged_params:

  Number of unstable parameters.

- rmse_trend:

  Character: `"decreasing"`, `"stable"`, or `"increasing"` based on a
  linear trend of nRMSE over iterations.

- threshold_pct:

  The threshold used.

## Details

- stable:

  Low variability across iterations.

- drifted:

  Monotonic shift to a new value (few reversals); the parameter moved
  but settled. This is expected behaviour when true outliers are
  removed.

- unstable:

  RSE \> 50% with at least three trend reversals, or CV of the final
  third of iterations \> 10%. Warrants investigation.

## See also

\[plot.irxclean_stability()\], \[stability_flag_table()\],
\[remove_erroneous_obs()\]

## Examples

``` r
if (FALSE) { # \dontrun{
results   <- load_removal_results("err_rem1")
stability <- check_model_stability(results)
print(stability)
plot(stability)
stability_flag_table(stability)
} # }
```
