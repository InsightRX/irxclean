# Parse THETA and OMEGA parameter labels from NONMEM .mod file comments

Extracts human-readable labels from inline comments in the `$THETA` and
`$OMEGA` blocks of a NONMEM control stream. Labels are used
automatically as panel titles in
[`plot_removal_metrics()`](https://insightrx.github.io/irxclean/reference/plot_removal_metrics.md).

## Usage

``` r
nm_parse_param_labels(mod_file)
```

## Arguments

- mod_file:

  Path to the NONMEM `.mod` file.

## Value

A named list with elements `thetas` and `omegas`, each a named character
vector mapping column names (e.g. `"THETA1"`, `"OMEGA.3.3"`) to
human-readable labels. Returns
`list(thetas = character(0), omegas = character(0))` if no labels are
found.

## Labelling your .mod file for irxclean

To get readable panel titles in irxclean plots, add inline comments to
the `$THETA` and `$OMEGA` blocks of your NONMEM control stream.

**Format:**

      ; <index>. <Label>          -- single parameter
      ; <start>-<end>. <Label>    -- range (e.g. IOV block)

The number(s) before the dot correspond to the 1-based parameter index
(the same numbering NONMEM uses in the `.ext` file). The dot and
surrounding spaces are flexible – the parser accepts `; 1. CL`,
`; 1.CL`, and `; 1 CL`.

**THETA example:**

    $THETA
      (0, 11.5555) ; 1. CL
      (0,  9.2)    ; 2. V1
      (0,  3.1)    ; 3. KA

**OMEGA example (IIV + IOV):**

    $OMEGA
      0.09         ; 1. IIV CL
      0.04         ; 2. IIV V1

    $OMEGA BLOCK(1)
      0.03         ; 3. IOV CL
    $OMEGA BLOCK(1) SAME   ; -- no label needed, auto-excluded from plots
    $OMEGA BLOCK(1) SAME

**Rules:**

- **THETA**: label every row you want named; unlabelled rows keep their
  default name (e.g. `THETA4`).

- **OMEGA**: label only the *first* element of each unique term.
  `BLOCK(1) SAME` lines are detected automatically and always excluded
  from plots regardless of any comment on that line.

- Range notation (`; 3-7. IOV CL`) is supported; only the first element
  (`OMEGA.3.3`) is labelled and plotted.

- The parser is case-insensitive for the `SAME` keyword and tolerates
  varied spacing around the dash and dot.
