# Iteratively detect and remove influential individuals

Performs iterative NONMEM-based detection of individuals whose removal
most improves overall model fit, as measured by the change in total
objective function value (dOFV). The method extends the influential
individual approach of Karlsson by iterating: after each removal the
model is re-estimated on the reduced dataset and the next most
influential subject is identified from the updated individual OBJ
values.

## Usage

``` r
remove_influential_individuals(
  dat,
  mod,
  run_id,
  n = NULL,
  dofv_threshold = NULL,
  id_col = "ID",
  verbose = TRUE,
  save_results = TRUE,
  save_temp_dir = FALSE,
  stability_check = TRUE
)
```

## Arguments

- dat:

  File path to the NONMEM-ready `.csv` dataset (with or without the
  `.csv` extension, relative or absolute), **or** an R data frame that
  is already loaded in memory.

- mod:

  File path to the NONMEM `.mod` control stream (with or without the
  `.mod` extension, relative or absolute), **or** a `nm_model` list
  returned by
  [`nm_read_model()`](https://insightrx.github.io/irxclean/reference/nm_read_model.md).

- run_id:

  Character. Unique run identifier (e.g. `"iir_run1"`).

- n:

  Integer or `NULL`. Maximum number of subjects to remove. Must specify
  at least one of `n` or `dofv_threshold`. If both are supplied the run
  stops when *either* criterion is met.

- dofv_threshold:

  Numeric or `NULL`. Stop when the change in sumOFV from removing a
  subject (`sumOFV_before - sumOFV_after`) falls below this value. A
  common choice is 3.84 (chi-squared 1 df, p \< 0.05). When the
  threshold is not met the last entry in `rem` retains its `dOFV` value
  so the user can inspect it before deciding whether to apply the
  removal. Must specify at least one of `n` or `dofv_threshold`.

- id_col:

  Character. Name of the subject ID column in the dataset. Default
  `"ID"`.

- verbose:

  Logical. Print progress messages? Default `TRUE`.

- save_results:

  Logical. Save an `.RDS` results file to disk? Default `TRUE`.

- save_temp_dir:

  Logical. Copy the temporary working directory to the current working
  directory (`TRUE`) or delete it after completion (`FALSE`, default)?

- stability_check:

  Logical. Append model stability metrics to the returned results?
  Default `TRUE`. See
  [`check_model_stability()`](https://insightrx.github.io/irxclean/reference/check_model_stability.md).

## Value

A list containing:

- par:

  Data frame of population parameter estimates per iteration. ITERATION
  1 = baseline (full dataset); ITERATION 2..k+1 = after each of k
  subject removals.

- rem:

  Data frame of removed subjects with columns: `ID`, `iOFV` (individual
  OBJ before removal), `sumOFV_before`, `sumOFV_after`, `dOFV`, and
  `ITERATION` (removal sequence number, 1-based).

- phi:

  Data frame of individual OBJ values per NONMEM run iteration (columns
  `ID`, `OBJ`, `ITERATION`). Removed subjects are absent from subsequent
  iterations.

- ofv:

  Data frame of `sumOFV` and `dOFV` per NONMEM run iteration. `dOFV` at
  ITERATION 1 is `NA` (baseline).

- rmse:

  Data frame of normalised RMSE per iteration.

- stability:

  (If `stability_check = TRUE`) Stability assessment from
  [`check_model_stability()`](https://insightrx.github.io/irxclean/reference/check_model_stability.md).

- param_labels:

  Named list with elements `thetas` and `omegas`: human-readable
  parameter labels parsed from the `.mod` file comments, used
  automatically in plots.

## Details

At each iteration the algorithm:

1.  Runs NONMEM (via PsN `execute`) on the current dataset.

2.  Reads the individual OBJ values from the `.phi` file and computes
    `sumOFV`.

3.  (From iteration 2 onwards) Calculates
    `dOFV = sumOFV(prev) - sumOFV(curr)` and checks the stopping
    criterion.

4.  Identifies the subject with the highest individual OBJ (iOFV).

5.  Removes all rows for that subject from the dataset.

6.  Saves parameter estimates and diagnostics for later plotting.

## Stopping behaviour

When `dofv_threshold` is the active stopping criterion, the algorithm
runs one additional NONMEM estimation after each removal to obtain the
updated `sumOFV`. The iteration that first produces
`dOFV < dofv_threshold` is recorded in `rem` with its `dOFV` value; no
further removals are performed. Users should inspect that final entry
before deciding whether to apply it.

## Supported input modes

The same three input modes as
[`remove_erroneous_obs()`](https://insightrx.github.io/irxclean/reference/remove_erroneous_obs.md)
are supported: bare file name, relative/absolute file path, or in-memory
R objects. See
[`remove_erroneous_obs()`](https://insightrx.github.io/irxclean/reference/remove_erroneous_obs.md)
for full details.

## PsN requirement

`execute` and `sumo` must be available on the system `PATH`.

## See also

[`remove_erroneous_obs()`](https://insightrx.github.io/irxclean/reference/remove_erroneous_obs.md),
[`plot_removal_metrics()`](https://insightrx.github.io/irxclean/reference/plot_removal_metrics.md),
[`check_model_stability()`](https://insightrx.github.io/irxclean/reference/check_model_stability.md),
[`remove_observations_from_data()`](https://insightrx.github.io/irxclean/reference/remove_observations_from_data.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Remove up to 10 subjects
results <- remove_influential_individuals(
  dat    = "my_data",
  mod    = "my_model",
  run_id = "iir_run1",
  n      = 10
)

# Stop when dOFV drops below chi-squared threshold (p < 0.05, 1 df)
results <- remove_influential_individuals(
  dat            = "my_data",
  mod            = "my_model",
  run_id         = "iir_run1",
  dofv_threshold = 3.84
)

# Both: stop at whichever criterion is hit first
results <- remove_influential_individuals(
  dat            = "my_data",
  mod            = "my_model",
  run_id         = "iir_run1",
  n              = 20,
  dofv_threshold = 3.84
)

# Inspect removals
results$rem
results$ofv

# Reuse existing plot functions for parameter / RMSE diagnostics
plot_removal_metrics(results, metric = "thetas")
plot_removal_metrics(results, metric = "rmse")
} # }
```
