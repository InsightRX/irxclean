# Generate an automated HTML cleaning report

Renders an R Markdown report summarising the iterative outlier-removal
workflow. The report is structured into two sections:

1.  **Key Results** — always visible: observations removed, pOFV change,
    error-parameter stability (if detected), and nRMSE trend.

2.  **Detailed Results** — collapsible: full THETA, OMEGA, and SIGMA
    stability facet plots plus the model stability assessment.

An optional demographic comparison section is rendered when `data` and
covariate columns are supplied.

## Usage

``` r
generate_report(
  results,
  data = NULL,
  prop_error_param = NULL,
  add_error_param = NULL,
  continuous_cov_cols = NULL,
  categorical_cov_cols = NULL,
  id_col = "ID",
  time_col = "TIME",
  evid_col = "EVID",
  tad_col = NULL,
  tad_bin_width = NULL,
  stability_threshold_pct = 20,
  title = "irxclean Data Cleaning Report",
  output_file = "irxclean_report.html",
  output_dir = ".",
  open_report = FALSE,
  quiet = TRUE
)
```

## Arguments

- results:

  A results list from \[remove_erroneous_obs()\] or
  \[load_removal_results()\]. \*\*Required.\*\*

- data:

  A data frame containing the original NONMEM dataset. Required for
  demographic comparison; \`NULL\` skips that section.

- prop_error_param:

  Character or `NULL`. Column name of the proportional residual error
  parameter to highlight in the Key Results section. Accepts the column
  name as stored in `results$par` (e.g. `"THETA7"` or `"SIGMA.1.1"`), or
  NONMEM index notation such as `"THETA(7)"` or `"SIGMA(1)"` which is
  normalised automatically. `NULL` triggers auto-detection from
  `results$param_labels`; the detected column is printed when
  `quiet = FALSE`.

- add_error_param:

  Character or `NULL`. Column name of the additive residual error
  parameter. Accepts the same formats as `prop_error_param`. `NULL`
  triggers auto-detection.

- continuous_cov_cols:

  Character vector of continuous covariate column names for demographic
  comparison. \`NULL\` skips continuous covariates.

- categorical_cov_cols:

  Character vector of categorical covariate column names. \`NULL\` skips
  categorical covariates.

- id_col:

  Character. Subject ID column name. Default \`"ID"\`.

- time_col:

  Character. Time column for TAD computation. Default `"TIME"`.

- evid_col:

  Character. EVID column for identifying dose records. Default `"EVID"`.

- tad_col:

  Character or `NULL`. Pre-computed TAD column name; `NULL` triggers
  auto-detection / computation. Default `NULL`.

- tad_bin_width:

  Numeric or `NULL`. TAD bin width for the curve-position plot. `NULL`
  auto-selects. Default `NULL`.

- stability_threshold_pct:

  Numeric. Parameter drift threshold (%) passed to
  \[check_model_stability()\]. Default \`20\`.

- title:

  Character. Report title. Default \`"irxclean Data Cleaning Report"\`.

- output_file:

  Character. Output file name. Default \`"irxclean_report.html"\`.

- output_dir:

  Character. Directory for the output file. Default: current working
  directory.

- open_report:

  Logical. Open the report in the default browser after rendering?
  Default \`FALSE\`.

- quiet:

  Logical. Suppress rmarkdown rendering messages? Default \`TRUE\`.

## Value

Invisibly returns the path to the rendered HTML file.

## Requirements

Requires the \`rmarkdown\` and \`knitr\` packages (listed in
\`Suggests\`).

## Examples

``` r
if (FALSE) { # \dontrun{
results <- load_removal_results("err_rem1")
data    <- read.csv("my_data.csv")

generate_report(
  results              = results,
  data                 = data,
  continuous_cov_cols  = c("AGE", "WT", "HT"),
  categorical_cov_cols = "SEX",
  title                = "Busulfan PK Cleaning Report"
)
} # }
```
