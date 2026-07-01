# Run a single ferx iteration and return standardized results

Run a single ferx iteration and return standardized results

## Usage

``` r
ferx_run_iteration(
  model_path,
  data,
  run_id,
  iteration,
  method = "focei",
  verbose = TRUE,
  ...
)
```

## Arguments

- model_path:

  Path to the \`.ferx\` model file.

- data:

  Data frame to fit.

- run_id:

  Character. Run identifier (used for file naming).

- iteration:

  Integer. Current iteration number.

- method:

  Character. Estimation method passed to \`ferx::ferx_fit()\`.

- verbose:

  Logical. Print progress?

- ...:

  Additional arguments passed to \`ferx::ferx_fit()\`.

## Value

A named list with elements: \`sdtab\`, \`theta\`, \`omega\`, \`sigma\`,
\`ofv\`, \`individual_obj\`.
