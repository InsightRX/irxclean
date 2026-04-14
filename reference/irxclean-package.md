# irxclean: Model-Informed Detection of Potentially Erroneous Pharmacokinetic Observations

\`irxclean\` provides a suite of tools for detecting and flagging
potentially erroneous pharmacokinetic (PK) observations in real-world
clinical datasets for user investigation and review. Because assessment
of potentially erroneous data is inherently model-dependent in
population PK/PD analyses, the package centres on an iterative approach
that uses NONMEM (via Perl-speaks-NONMEM, PsN) to re-estimate the model
after each observation is removed, surfacing the observations that most
distort the fitted model.

\## Workflow

A typical \`irxclean\` analysis proceeds in four stages:

1\. \*\*Gross cleaning\*\* - model-independent assessment of raw data
quality: covariate distribution outlier screening, detection of
concentration elevations without an intervening dose, and missing-data
summaries. See \[assess_data_quality()\].

2\. \*\*Iterative flagging\*\* - CWRES-based iterative identification of
the worst-fitting observation at each NONMEM estimation, producing a
ranked list of potentially erroneous observations for user review. See
\[remove_erroneous_obs()\].

3\. \*\*Stability & metrics\*\* - track parameter stability, pOFV, RMSE,
and flag potential model instabilities across iterations. See
\[check_model_stability()\] and \[plot_removal_metrics()\].

4\. \*\*Demographic comparison\*\* - compare subjects with flagged
observations against the full dataset to screen for demographic bias.
See \[plot_demographic_comparison()\].

After reviewing flagged observations,
\[remove_observations_from_data()\] applies the user-approved exclusions
to produce the final cleaned dataset. A single \[generate_report()\]
call produces an HTML summary of all findings.

\## NONMEM / PsN requirement

\[remove_erroneous_obs()\] requires a working PsN installation
(\`execute\`, \`sumo\`) accessible on the system \`PATH\`. All other
functions work on pre-computed results and do not require NONMEM or PsN.

## See also

Useful links:

- <https://github.com/InsightRX/irxclean>

- Report bugs at <https://github.com/InsightRX/irxclean/issues>

## Author

**Maintainer**: Jordan Brooks <jordan.brooks@insight-rx.com>

Other contributors:

- InsightRX \[copyright holder\]
