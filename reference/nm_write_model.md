# Write a NONMEM model object to a \`.mod\` file

Write a NONMEM model object to a \`.mod\` file

## Usage

``` r
nm_write_model(model = NULL, modelfile = NULL, overwrite = FALSE)
```

## Arguments

- model:

  A NONMEM model object created by \[nm_read_model()\].

- modelfile:

  Output file path.

- overwrite:

  Overwrite an existing file? Default \`FALSE\`.

## Value

Invisibly returns `modelfile` (the output file path).
