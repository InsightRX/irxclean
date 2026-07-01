# Convert ferx fit parameters to a flat named list

Normalizes \`fit\$theta\`, \`fit\$omega\`, \`fit\$sigma\` into the same
column-name format that \[nm_read_pars()\] produces (THETA1, THETA2,
..., OMEGA.1.1, ..., SIGMA.1.1, ..., OBJ).

## Usage

``` r
ferx_read_pars(fit_result)
```

## Arguments

- fit_result:

  List returned by \[ferx_run_iteration()\].

## Value

A named list of scalar parameter values.
