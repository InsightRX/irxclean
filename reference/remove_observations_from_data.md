# Apply user-reviewed observation exclusions to produce a cleaned dataset

After running
[`remove_erroneous_obs()`](https://insightrx.github.io/irxclean/reference/remove_erroneous_obs.md)
and reviewing the flagged observations, use this function to apply
exclusions and produce a cleaned dataset suitable for final analysis.
The top \`n\` flagged observations (by iteration order) are excluded.

## Usage

``` r
remove_observations_from_data(
  data,
  results,
  n = NULL,
  columns = list(ID = "ID", TIME = "TIME", DV = "DV")
)
```

## Arguments

- data:

  A data frame containing the original NONMEM dataset.

- results:

  A results list from \[remove_erroneous_obs()\] or
  \[load_removal_results()\].

- n:

  Integer. Number of observations to remove (taken from the first \`n\`
  iterations in \`results\$rem\`). Default: all removed observations.

- columns:

  Named list specifying column labels for matching. Default: \`list(ID =
  "ID", TIME = "TIME", DV = "DV")\`.

## Value

The input \`data\` frame with the identified observations removed.

## User review

Flagged observations should be investigated before calling this
function. A flagged observation is not necessarily erroneous; it may
represent a genuine patient population subgroup or an unusual but valid
measurement. Final exclusion decisions rest with the analyst.

Matching is performed on rounded (\`3\` decimal places) \`ID\`,
\`TIME\`, and \`DV\` values. A warning is raised if the number of rows
removed does not equal \`n\`.

## Examples

``` r
if (FALSE) { # \dontrun{
original_data <- read.csv("my_data.csv")
results       <- load_removal_results("err_rem1")
cleaned       <- remove_observations_from_data(original_data, results, n = 5)
} # }
```
