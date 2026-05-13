# Changelog

## irxclean (development version)

## irxclean 0.1.0

- Initial release.
- Gross data quality assessment: covariate distribution outlier
  detection and concentration anomaly checks
  ([`assess_data_quality()`](https://insightrx.github.io/irxclean/reference/assess_data_quality.md)).
- Model-based iterative CWRES flagging via PsN/NONMEM
  ([`remove_erroneous_obs()`](https://insightrx.github.io/irxclean/reference/remove_erroneous_obs.md)).
- Model stability monitoring across iterations
  ([`check_model_stability()`](https://insightrx.github.io/irxclean/reference/check_model_stability.md)).
- Exclusion criteria application and tracking
  ([`apply_exclusion_criteria()`](https://insightrx.github.io/irxclean/reference/apply_exclusion_criteria.md)).
- Demographic comparison of flagged vs. retained subjects
  ([`plot_demographic_comparison()`](https://insightrx.github.io/irxclean/reference/plot_demographic_comparison.md)).
- Automated HTML report generation
  ([`generate_report()`](https://insightrx.github.io/irxclean/reference/generate_report.md)).
