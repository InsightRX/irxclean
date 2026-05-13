# Plot metrics from the observation removal process

Creates diagnostic plots tracking the effect of iterative observation
removal on model fit and parameter stability.

## Usage

``` r
plot_removal_metrics(
  results,
  metric = c("pOFV", "thetas", "omegas", "sigmas", "rmse"),
  verbose = TRUE
)
```

## Arguments

- results:

  A results list from \[remove_erroneous_obs()\] or
  \[load_removal_results()\].

- metric:

  Character. Which metric to plot:

  \`"pOFV"\`

  :   Change in pseudo-OFV at each iteration.

  \`"thetas"\`

  :   THETA estimates normalised to the final estimate.

  \`"omegas"\`

  :   OMEGA estimates normalised to the final estimate.

  \`"sigmas"\`

  :   SIGMA estimates normalised to the final estimate. Returns `NULL`
      with a message if all SIGMA parameters are fixed.

  \`"rmse"\`

  :   Root mean square error over iterations.

- verbose:

  Logical. Print progress messages? Default \`TRUE\`.

## Value

A \`ggplot\` object.

## Parameter labels

When the results object contains a \`param_labels\` element (populated
automatically by \[remove_erroneous_obs()\]), the THETA and OMEGA plot
panels are labelled with the human-readable names parsed from the `.mod`
file.

To enable labelling, add inline comments to your `$THETA` and `$OMEGA`
blocks following the format `; <index>. <Label>`:


    $THETA
      (0, 11.6) ; 1. CL
      (0,  9.2) ; 2. V1

    $OMEGA
      0.09      ; 1. IIV CL
      0.04      ; 2. IIV V1

    $OMEGA BLOCK(1)
      0.03      ; 3. IOV CL
    $OMEGA BLOCK(1) SAME

`BLOCK(1) SAME` entries are detected automatically and excluded from the
OMEGA plot regardless of whether a comment is present. Range notation
(`; 3-7. IOV CL`) is also supported; only the first diagonal element is
plotted.

## Examples

``` r
if (FALSE) { # \dontrun{
results <- load_removal_results("err_rem1")
plot_removal_metrics(results, metric = "pOFV")
plot_removal_metrics(results, metric = "thetas")
} # }
```
