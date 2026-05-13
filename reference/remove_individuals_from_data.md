# Apply user-reviewed subject exclusions to produce a cleaned dataset

After running
[`remove_influential_individuals()`](https://insightrx.github.io/irxclean/reference/remove_influential_individuals.md)
and reviewing the flagged subjects, use this function to apply
exclusions at the *subject* level and produce a cleaned dataset suitable
for final analysis. All rows belonging to the top `n` flagged subjects
(by iteration order) are removed.

## Usage

``` r
remove_individuals_from_data(data, results, n = NULL, id_col = "ID")
```

## Arguments

- data:

  A data frame containing the original NONMEM dataset.

- results:

  A results list from
  [`remove_influential_individuals()`](https://insightrx.github.io/irxclean/reference/remove_influential_individuals.md).

- n:

  Integer or `NULL`. Number of subjects to remove (taken from the first
  `n` entries in `results$rem`, ordered by `ITERATION`). Default: all
  removed subjects.

- id_col:

  Character. Name of the subject ID column. Default `"ID"`.

## Value

The input `data` frame with all rows for the excluded subjects removed.
A warning is raised if the resulting number of unique IDs does not
decrease by exactly `n`.

## Details

This is the subject-level analogue of
[`remove_observations_from_data()`](https://insightrx.github.io/irxclean/reference/remove_observations_from_data.md).

## User review

Flagged subjects should be investigated before calling this function. A
high individual OFV may reflect a genuine patient subgroup or an unusual
but valid pharmacokinetic profile, not necessarily erroneous data. Final
exclusion decisions rest with the analyst.

## See also

[`remove_influential_individuals()`](https://insightrx.github.io/irxclean/reference/remove_influential_individuals.md),
[`remove_observations_from_data()`](https://insightrx.github.io/irxclean/reference/remove_observations_from_data.md)

## Examples

``` r
if (FALSE) { # \dontrun{
original_data <- read.csv("my_data.csv")
results       <- remove_influential_individuals(
  dat    = original_data,
  mod    = "my_model",
  run_id = "iir_run1",
  n      = 10
)
# Remove the top 5 most influential subjects
cleaned <- remove_individuals_from_data(original_data, results, n = 5)
} # }
```
