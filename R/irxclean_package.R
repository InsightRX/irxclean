#' irxclean: Model-Informed Detection of Potentially Erroneous Pharmacokinetic Observations
#'
#' @description
#' `irxclean` provides a suite of tools for detecting and flagging potentially
#' erroneous pharmacokinetic (PK) observations in real-world clinical datasets
#' for user investigation and review. Because assessment of potentially erroneous
#' data is inherently model-dependent in population PK/PD analyses, the package
#' centres on an iterative approach that uses NONMEM (via Perl-speaks-NONMEM,
#' PsN) to re-estimate the model after each observation is removed, surfacing
#' the observations that most distort the fitted model.
#'
#' ## Workflow
#'
#' A typical `irxclean` analysis proceeds in four stages:
#'
#' 1. **Gross cleaning** - model-independent assessment of raw data quality:
#'    covariate distribution outlier screening, detection of concentration
#'    elevations without an intervening dose, and missing-data summaries.
#'    See [assess_data_quality()].
#'
#' 2. **Iterative flagging** - CWRES-based iterative identification of the
#'    worst-fitting observation at each NONMEM estimation, producing a ranked
#'    list of potentially erroneous observations for user review.
#'    See [remove_erroneous_obs()].
#'
#' 3. **Stability & metrics** - track parameter stability, pOFV, RMSE, and
#'    flag potential model instabilities across iterations.  See
#'    [check_model_stability()] and [plot_removal_metrics()].
#'
#' 4. **Demographic comparison** - compare subjects with flagged observations
#'    against the full dataset to screen for demographic bias.
#'    See [plot_demographic_comparison()].
#'
#' After reviewing flagged observations, [remove_observations_from_data()]
#' applies the user-approved exclusions to produce the final cleaned dataset.
#' A single [generate_report()] call produces an HTML summary of all findings.
#'
#' ## NONMEM / PsN requirement
#'
#' [remove_erroneous_obs()] requires a working PsN installation (`execute`,
#' `sumo`) accessible on the system `PATH`.  All other functions work on
#' pre-computed results and do not require NONMEM or PsN.
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom dplyr arrange bind_rows case_when filter first lag mutate select
#'   slice distinct across where
#' @importFrom ggplot2 aes element_blank element_line element_text facet_wrap
#'   geom_boxplot geom_bar geom_density geom_hline geom_line geom_point
#'   geom_violin geom_jitter geom_vline ggplot ggtitle labs scale_fill_manual
#'   scale_colour_manual scale_x_continuous expansion annotate stat_smooth
#'   theme ylab xlab theme_bw position_dodge geom_text geom_col coord_flip
#' @importFrom patchwork wrap_plots
#' @importFrom tidyr pivot_longer
#' @importFrom stringr str_detect str_replace str_replace_all str_split str_c
#' @importFrom rlang .data
#' @importFrom stats median quantile wilcox.test chisq.test binom.test sd IQR setNames
#' @importFrom utils read.csv read.table write.csv tail
#' @importFrom vpc read_table_nm
#' @importFrom methods is
## usethis namespace: end
NULL
