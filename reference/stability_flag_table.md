# Summarise unstable parameters in a tidy table

Returns the rows of the parameter summary for parameters classified as
`"unstable"` (i.e. those that did not reach a stable estimate).

## Usage

``` r
stability_flag_table(stability)
```

## Arguments

- stability:

  An \`irxclean_stability\` object.

## Value

A data frame of unstable parameters with stability metrics.
