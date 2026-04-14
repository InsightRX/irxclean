# Iteratively detect potentially erroneous pharmacokinetic observations

Performs iterative NONMEM-based detection of potentially erroneous
observations in a pharmacokinetic dataset. The ranked list of flagged
observations is intended for user review before any exclusions are
applied to the final analysis dataset (see
[`remove_observations_from_data()`](https://insightrx.github.io/irxclean/reference/remove_observations_from_data.md)).
At each iteration the algorithm:

1.  Runs NONMEM (via PsN \`execute\`) on the current dataset.

2.  Identifies the observation with the highest absolute CWRES.

3.  Removes that observation from the dataset.

4.  Saves parameter estimates and diagnostics for later plotting.

## Usage

``` r
remove_erroneous_obs(
  dat,
  mod,
  run_id,
  n,
  columns = list(ID = "ID", TIME = "TIME", DV = "DV"),
  verbose = TRUE,
  save_results = TRUE,
  save_temp_dir = FALSE,
  stability_check = TRUE
)
```

## Arguments

- dat:

  File path to the NONMEM-ready \`.csv\` dataset (with or without the
  \`.csv\` extension, relative or absolute), \*\*or\*\* an R data frame
  that is already loaded in memory.

- mod:

  File path to the NONMEM \`.mod\` control stream (with or without the
  \`.mod\` extension, relative or absolute), \*\*or\*\* a `nm_model`
  list returned by
  [`nm_read_model()`](https://insightrx.github.io/irxclean/reference/nm_read_model.md).

- run_id:

  Character. Unique run identifier (e.g. \`"err_rem1"\`).

- n:

  Integer. Number of observations to remove.

- columns:

  Named list specifying column labels in the dataset. Default: \`list(ID
  = "ID", TIME = "TIME", DV = "DV")\`.

- verbose:

  Logical. Print progress messages? Default \`TRUE\`.

- save_results:

  Logical. Save an \`.RDS\` results file to disk? Default \`TRUE\`.

- save_temp_dir:

  Logical. Copy the temporary working directory to the current working
  directory (\`TRUE\`) or delete it after completion (\`FALSE\`,
  default)?

- stability_check:

  Logical. Append model stability metrics to the returned results?
  Default \`TRUE\`. See \[check_model_stability()\].

## Value

A list containing:

- par:

  Data frame of population parameter estimates per iteration.

- rem:

  Data frame of removed observations with their iteration.

- phi:

  Data frame of individual objective function (OBJ) values per
  iteration.

- rmse:

  Data frame of RMSE per iteration.

- stability:

  (If `stability_check = TRUE`) Stability assessment from
  [`check_model_stability()`](https://insightrx.github.io/irxclean/reference/check_model_stability.md).

- param_labels:

  Named list with elements `thetas` and `omegas`: human-readable
  parameter labels parsed from the `.mod` file comments, used
  automatically in plots.

## Specifying inputs – three supported modes

**Mode 1 – bare file name (files in the working directory)**

Pass the stem of the filename without the extension. The files must
exist in the current working directory. This was the original interface
and remains the simplest option when your working directory is already
set to the folder that contains the data and model.

    results <- remove_erroneous_obs(
      dat    = "my_data",    # reads my_data.csv from getwd()
      mod    = "my_model",   # reads my_model.mod from getwd()
      run_id = "clean_run1",
      n      = 10
    )

**Mode 2 – file path (relative or absolute, extension optional)**

Pass a relative or absolute path. The `.csv` / `.mod` extension is
stripped and re-appended automatically, so you can include or omit it.
Useful when files live in a different directory from your R session.

    results <- remove_erroneous_obs(
      dat    = "/projects/busulfan/data/bu_data.csv",
      mod    = "/projects/busulfan/models/bu_base",
      run_id = "clean_run1",
      n      = 10
    )

**Mode 3 – R objects already loaded in memory**

Pass a `data.frame` for `dat` and/or a `nm_model` list (from
[`nm_read_model()`](https://insightrx.github.io/irxclean/reference/nm_read_model.md))
for `mod`. This is useful when you have already read and pre-processed
your data in R (e.g. after applying
[`apply_exclusion_criteria()`](https://insightrx.github.io/irxclean/reference/apply_exclusion_criteria.md)),
or when you want to modify the model programmatically before running.
NONMEM still runs from temporary files on disk – the objects are written
there automatically.

    pk_clean <- apply_exclusion_criteria(...)$data_clean
    nm_mod   <- nm_read_model("/projects/busulfan/models/bu_base.mod")

    results <- remove_erroneous_obs(
      dat    = pk_clean,   # data.frame -- no file needed
      mod    = nm_mod,     # nm_model list -- no file copy needed
      run_id = "clean_run1",
      n      = 10
    )

The three modes can be mixed: you may supply an in-memory data frame for
`dat` while pointing `mod` at a file path, or vice versa.

## PsN requirement

\`execute\` and \`sumo\` must be available on the system \`PATH\`.

## See also

\[plot_removal_metrics()\], \[check_model_stability()\],
\[remove_observations_from_data()\], \[generate_report()\]

## Examples

``` r
if (FALSE) { # \dontrun{
# Mode 1: bare file names (files in working directory)
results <- remove_erroneous_obs(
  dat    = "my_data",
  mod    = "my_model",
  run_id = "err_rem1",
  n      = 5
)

# Mode 2: full file paths (extension optional)
results <- remove_erroneous_obs(
  dat    = "/projects/pk/my_data.csv",
  mod    = "/projects/pk/my_model.mod",
  run_id = "err_rem1",
  n      = 5
)

# Mode 3: R objects already in memory
pk_data  <- read.csv("/projects/pk/my_data.csv")
nm_model <- nm_read_model("/projects/pk/my_model.mod")
results  <- remove_erroneous_obs(
  dat    = pk_data,
  mod    = nm_model,
  run_id = "err_rem1",
  n      = 5
)

# Mix and match: in-memory data + file path for model
pk_clean <- apply_exclusion_criteria(pk_data)$data_clean
results  <- remove_erroneous_obs(
  dat    = pk_clean,
  mod    = "/projects/pk/my_model",
  run_id = "err_rem1",
  n      = 5
)

plot_removal_metrics(results, metric = "pOFV")
} # }
```
