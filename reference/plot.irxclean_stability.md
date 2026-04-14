# Plot method for irxclean_stability objects

Produces two panels:

1.  Parameter trajectories (% change from baseline) coloured by
    stability status (stable = teal, drifted = grey, unstable = red).

2.  RSE bar chart per parameter, sorted by RSE and coloured by status.

## Usage

``` r
# S3 method for class 'irxclean_stability'
plot(x, ...)
```

## Arguments

- x:

  An \`irxclean_stability\` object from \[check_model_stability()\].

- ...:

  Ignored.

## Value

A list of `ggplot` objects: `$parameter_drift` and `$param_rse`.
