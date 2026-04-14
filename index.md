# irxclean

Model-informed detection and iterative removal of potentially erroneous
pharmacokinetic observations in clinical data.

## Description

`irxclean` provides tools for flagging potentially erroneous
pharmacokinetic (PK) observations for user investigation and review. The
package supports two complementary approaches:

- **Gross data quality assessment**: model-independent checks for
  covariate distribution outliers and concentration anomalies
  (e.g. negative values, implausible timing, magnitude outliers).
- **Model-based iterative flagging**: iterative CWRES-based
  identification of potentially erroneous observations using
  Perl-speaks-NONMEM (PsN) and NONMEM. Includes model stability
  monitoring across iterations and demographic comparison of subjects
  with flagged observations versus the full dataset.

An automated HTML report summarising all findings can be generated with
[`generate_report()`](https://insightrx.github.io/irxclean/reference/generate_report.md).
Final exclusion decisions remain with the user after reviewing flagged
observations.

## Installation

``` r
devtools::install_github("InsightRX/irxclean")
```

## Requirements

The model-based iterative flagging functionality
([`remove_erroneous_obs()`](https://insightrx.github.io/irxclean/reference/remove_erroneous_obs.md))
requires working installations of
[NONMEM](https://www.iconplc.com/solutions/technologies/nonmem/) and
[PsN](https://uupharmacometrics.github.io/PsN/).

## Contributing

We welcome input from the community:

- If you think you have encountered a bug, please [submit an
  issue](https://github.com/InsightRX/irxclean/issues) on the GitHub
  page. Please include a reproducible example of the unexpected
  behavior.

- Please [open a pull
  request](https://github.com/InsightRX/irxclean/pulls) if you have a
  fix or updates that would improve the package. If you’re not sure if
  your proposed changes are useful or within scope of the package, feel
  free to contact one of the authors of this package.

## Disclaimer

The functionality in this R package is provided “as is”. While its
authors adhere to software development best practices, the software may
still contain unintended errors.

InsightRX Inc. and the authors of this package can not be held liable
for any damages resulting from any use of this software. By the use of
this software package, the user waives all warranties, expressed or
implied, including any warranties to the accuracy, quality or
suitability of InsightRX for any particular purpose, either medical or
non-medical.

© ![InsightRX logo](reference/figures/insightrx_logo_color.png)
