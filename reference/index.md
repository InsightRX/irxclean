# Package index

## Gross data quality

- [`assess_data_quality()`](https://insightrx.github.io/irxclean/reference/assess_data_quality.md)
  : Assess data quality before model-dependent cleaning
- [`plot(`*`<irxclean_quality>`*`)`](https://insightrx.github.io/irxclean/reference/plot.irxclean_quality.md)
  : Plot method for irxclean_quality objects
- [`print(`*`<irxclean_quality>`*`)`](https://insightrx.github.io/irxclean/reference/print.irxclean_quality.md)
  : Print method for irxclean_quality objects

## Exclusion criteria

- [`apply_exclusion_criteria()`](https://insightrx.github.io/irxclean/reference/apply_exclusion_criteria.md)
  : Apply structured exclusion criteria to a NONMEM-format dataset
- [`print(`*`<irxclean_exclusions>`*`)`](https://insightrx.github.io/irxclean/reference/print.irxclean_exclusions.md)
  : Print method for irxclean_exclusions objects

## Model-based iterative cleaning

- [`remove_erroneous_obs()`](https://insightrx.github.io/irxclean/reference/remove_erroneous_obs.md)
  : Iteratively detect potentially erroneous pharmacokinetic
  observations
- [`remove_influential_individuals()`](https://insightrx.github.io/irxclean/reference/remove_influential_individuals.md)
  : Iteratively detect and remove influential individuals
- [`load_removal_results()`](https://insightrx.github.io/irxclean/reference/load_removal_results.md)
  : Load results from a previous \`remove_erroneous_obs()\` run
- [`calculate_pofv_change()`](https://insightrx.github.io/irxclean/reference/calculate_pofv_change.md)
  : Calculate the pseudo-objective function value (pOFV) change across
  iterations

## Model stability

- [`check_model_stability()`](https://insightrx.github.io/irxclean/reference/check_model_stability.md)
  : Check model stability across removal iterations
- [`stability_flag_table()`](https://insightrx.github.io/irxclean/reference/stability_flag_table.md)
  : Summarise unstable parameters in a tidy table
- [`plot(`*`<irxclean_stability>`*`)`](https://insightrx.github.io/irxclean/reference/plot.irxclean_stability.md)
  : Plot method for irxclean_stability objects
- [`print(`*`<irxclean_stability>`*`)`](https://insightrx.github.io/irxclean/reference/print.irxclean_stability.md)
  : Print method for irxclean_stability objects

## Dataset filtering

- [`remove_observations_from_data()`](https://insightrx.github.io/irxclean/reference/remove_observations_from_data.md)
  : Apply user-reviewed observation exclusions to produce a cleaned
  dataset
- [`remove_individuals_from_data()`](https://insightrx.github.io/irxclean/reference/remove_individuals_from_data.md)
  : Apply user-reviewed subject exclusions to produce a cleaned dataset

## Visualization & reporting

- [`plot_removal_metrics()`](https://insightrx.github.io/irxclean/reference/plot_removal_metrics.md)
  : Plot metrics from the observation removal process
- [`plot_influential_dofv()`](https://insightrx.github.io/irxclean/reference/plot_influential_dofv.md)
  : Plot dOFV trajectory from iterative influential individual removal
- [`plot_demographic_comparison()`](https://insightrx.github.io/irxclean/reference/plot_demographic_comparison.md)
  : Compare demographics of subjects with removed observations vs. whole
  dataset
- [`print(`*`<irxclean_demographics>`*`)`](https://insightrx.github.io/irxclean/reference/print.irxclean_demographics.md)
  : Print method for irxclean_demographics objects
- [`generate_report()`](https://insightrx.github.io/irxclean/reference/generate_report.md)
  : Generate an automated HTML cleaning report

## NONMEM utilities

- [`nm_read_model()`](https://insightrx.github.io/irxclean/reference/nm_read_model.md)
  : Parse a NONMEM model file into a named list of code blocks
- [`nm_write_model()`](https://insightrx.github.io/irxclean/reference/nm_write_model.md)
  : Write a NONMEM model object to a \`.mod\` file
- [`nm_read_pars()`](https://insightrx.github.io/irxclean/reference/nm_read_pars.md)
  : Read population parameter estimates from a NONMEM \`.ext\` file
- [`nm_parse_param_labels()`](https://insightrx.github.io/irxclean/reference/nm_parse_param_labels.md)
  : Parse THETA and OMEGA parameter labels from NONMEM .mod file
  comments

## Datasets

- [`busulfan_sim`](https://insightrx.github.io/irxclean/reference/busulfan_sim.md)
  : Simulated busulfan pharmacokinetic dataset
