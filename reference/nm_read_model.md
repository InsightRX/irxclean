# Parse a NONMEM model file into a named list of code blocks

Parse a NONMEM model file into a named list of code blocks

## Usage

``` r
nm_read_model(modelfile = NULL, as_block = FALSE, code = NULL)
```

## Arguments

- modelfile:

  Path to the NONMEM \`.mod\` file.

- as_block:

  If \`TRUE\` (default \`FALSE\`), each block is returned as a single
  collapsed string rather than a character vector of lines.

- code:

  Character string of NONMEM code (alternative to \`modelfile\`).

## Value

A named list of class \`c("NONMEM", "list")\`. Names correspond to
NONMEM block identifiers without the leading \`\$\` (e.g., \`"DATA"\`,
\`"THETA"\`). If the same block appears multiple times the entries are
concatenated.
