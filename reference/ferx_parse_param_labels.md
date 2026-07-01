# Parse parameter labels from a .ferx model file

Reads the \`\[parameters\]\` block of a \`.ferx\` file and extracts
named theta/omega/sigma labels. Returns the same structure as
\[nm_parse_param_labels()\].

## Usage

``` r
ferx_parse_param_labels(model_path)
```

## Arguments

- model_path:

  Path to the \`.ferx\` model file.

## Value

A named list with elements \`thetas\`, \`omegas\`, and \`sigmas\`, each
a named character vector mapping column names to labels.
