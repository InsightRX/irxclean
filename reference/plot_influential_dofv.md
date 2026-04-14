# Plot dOFV trajectory from iterative influential individual removal

Creates a waterfall/trajectory plot of the change in total OFV
(`dOFV = sumOFV_before - sumOFV_after`) for each subject removed by
[`remove_influential_individuals()`](https://insightrx.github.io/irxclean/reference/remove_influential_individuals.md).
An optional threshold line marks the stopping criterion when
`dofv_threshold` was used.

## Usage

``` r
plot_influential_dofv(results, dofv_threshold = NULL, label_ids = TRUE)
```

## Arguments

- results:

  A results list from
  [`remove_influential_individuals()`](https://insightrx.github.io/irxclean/reference/remove_influential_individuals.md).

- dofv_threshold:

  Numeric or `NULL`. If supplied, a horizontal reference line is drawn
  at this value (e.g. 3.84). Default `NULL`.

- label_ids:

  Logical. Annotate each point with the removed subject ID? Default
  `TRUE`.

## Value

A `ggplot` object.

## See also

[`remove_influential_individuals()`](https://insightrx.github.io/irxclean/reference/remove_influential_individuals.md)

## Examples

``` r
if (FALSE) { # \dontrun{
results <- remove_influential_individuals(
  dat = "my_data", mod = "my_model",
  run_id = "iir_run1", dofv_threshold = 3.84
)
plot_influential_dofv(results, dofv_threshold = 3.84)
} # }
```
