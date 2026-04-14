# irxclean 0.1.0

- Initial release.
- Gross data quality assessment: covariate distribution outlier detection and concentration anomaly checks (`assess_data_quality()`).
- Model-based iterative CWRES flagging via PsN/NONMEM (`remove_erroneous_obs()`).
- Model stability monitoring across iterations (`check_model_stability()`).
- Exclusion criteria application and tracking (`apply_exclusion_criteria()`).
- Demographic comparison of flagged vs. retained subjects (`plot_demographic_comparison()`).
- Automated HTML report generation (`generate_report()`).
