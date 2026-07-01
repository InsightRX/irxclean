# Read population parameter estimates from a NONMEM \`.ext\` file

Returns the final estimates (iteration `-1000000000`).

## Usage

``` r
nm_read_pars(model)
```

## Arguments

- model:

  Model name (character, e.g. \`"run1"\`) or run number (numeric). The
  \`.ext\` suffix is appended automatically if absent.

## Value

A named list of parameter values.
