# Plot method for irxclean_quality objects

Generates a panel of up to three diagnostic plots: covariate
distributions with outlier bounds, flagged concentration profiles, and
missing-data rates.

## Usage

``` r
# S3 method for class 'irxclean_quality'
plot(x, ...)
```

## Arguments

- x:

  An \`irxclean_quality\` object from \[assess_data_quality()\].

- ...:

  Ignored.

## Value

A \`ggplot\` object (or a list of \`ggplot\` objects).
