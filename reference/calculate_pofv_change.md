# Calculate the pseudo-objective function value (pOFV) change across iterations

The pOFV is defined as the sum of individual \\-2\log L_i\\ (OBJ) values
restricted to subjects that had \*\*no\*\* observations removed across
all iterations. The change in pOFV from one iteration to the next
reflects the genuine model-fit improvement attributable to each removal.

## Usage

``` r
calculate_pofv_change(results)
```

## Arguments

- results:

  A results list from \[remove_erroneous_obs()\] or
  \[load_removal_results()\].

## Value

A data frame with columns \`ITERATION\`, \`pofv\`, and \`diff_pofv\`.

## Examples

``` r
if (FALSE) { # \dontrun{
results <- load_removal_results("err_rem1")
pofv_df <- calculate_pofv_change(results)
} # }
```
