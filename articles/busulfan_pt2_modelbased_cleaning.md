# Simulated Busulfan PK: Part 2 — Model-Informed Iterative Cleaning

## Overview

> **Prerequisites:** This vignette continues from **Part 1**
> (`busulfan_pt1_gc_exclusion`). Run Part 1 first and save
> `busulfan_pt1_ready.csv` to your working directory before proceeding
> here. Part 1 covers dataset loading, concentration visualisation,
> gross QC with
> [`assess_data_quality()`](https://insightrx.github.io/irxclean/reference/assess_data_quality.md),
> covariate anomaly investigation, and structured exclusion criteria.

This vignette walks through the **model-informed** cleaning stages of an
`irxclean` analysis using the simulated busulfan pharmacokinetic
dataset. It covers Phase 2 of the full workflow: NONMEM-based iterative
observation flagging, parameter stability assessment, demographic
comparison, and HTML report generation.

### Intentional errors remaining after Phase 1

After Part 1 removes the subject with a weight decimal error (ID 18),
the following errors remain in the dataset for the model to detect:

| FLAG bit | Error type | Count |
|----|----|----|
| 2 | Concentration IQR outlier (\>3\*IQR above Q3 per TAD bin) | 2 observations |
| 4 | Timing error (+/- 5, 10, 30, or 60 minutes from true time) | 10 observations |
| 8 | Concentration magnitude error (multipliers 0.25–1.75x) | 10 observations |

------------------------------------------------------------------------

## 1. Load the Phase 1 dataset

Load `busulfan_pt1_ready.csv` saved at the end of Part 1. If that file
is not present in the working directory (e.g. when rendering from the
package without running Part 1 first), the raw simulated dataset is used
as a fallback — note that in this case one extra subject (ID 18) with a
covariate error remains in the data.

``` r

pt1_file <- "busulfan_pt1_ready.csv"

if (file.exists(pt1_file)) {
  dat <- read.csv(pt1_file)
  message("Loaded Phase 1 output: ", pt1_file)
} else {
  message(
    "busulfan_pt1_ready.csv not found -- falling back to raw dataset.\n",
    "Run busulfan_pt1_gc_exclusion first to generate the cleaned input."
  )
  dat <- busulfan_sim
}
#> busulfan_pt1_ready.csv not found -- falling back to raw dataset.
#> Run busulfan_pt1_gc_exclusion first to generate the cleaned input.

cat(sprintf(
  "Patients: %d  |  Rows: %d  |  Observations: %d\n",
  length(unique(dat$ID)), nrow(dat), sum(dat$EVID == 0)
))
#> Patients: 100  |  Rows: 1600  |  Observations: 1200
```

------------------------------------------------------------------------

## 2. Iterative observation flagging (NONMEM required)

[`remove_erroneous_obs()`](https://insightrx.github.io/irxclean/reference/remove_erroneous_obs.md)
fits NONMEM iteratively, removing the observation with the largest
\|CWRES\| at each step. Three input modes are supported — see
[`?remove_erroneous_obs`](https://insightrx.github.io/irxclean/reference/remove_erroneous_obs.md)
for details. Here we demonstrate **Mode 3**: pass the dataset as an
in-memory R data frame (no separate file copy needed) and point `mod` at
the bundled model file path.

Results are cached to `busulfan_results.RDS` in the working directory so
subsequent renders skip the NONMEM run. If neither a cache file nor a
PsN/NONMEM installation is present, the 5-iteration pre-computed result
bundled with the package is used as a fallback.

``` r

cache_file <- "busulfan_results.RDS"

if (file.exists(cache_file)) {
  results <- readRDS(cache_file)
  message("Loaded cached results from ", cache_file)
} else if (nzchar(Sys.which("execute"))) {
  # Mode 3: pass data frame directly; model given as a file path
  results <- remove_erroneous_obs(
    dat    = dat,   # Phase 1 cleaned data frame
    mod    = system.file("extdata", "pk_busulfan_ucsf_2025.mod",
                         package = "irxclean"),
    run_id = "busulfan_removal",
    n      = 20
  )
  saveRDS(results, cache_file)
  message("Results saved to ", cache_file)
} else {
  message("NONMEM/PsN not found and no cache -- loading package testdata.")
  results <- readRDS(system.file("extdata", "busulfan_results.RDS",
                                 package = "irxclean"))
}
#> NONMEM/PsN not found and no cache -- loading package testdata.

cat("Iterations completed:", nrow(results$rmse), "\n")
#> Iterations completed: 20
cat("Flagged observations (ranked by |CWRES|):\n")
#> Flagged observations (ranked by |CWRES|):
print(results$rem[, c("ID", "TIME", "DV", "CWRES")])
#>    ID  TIME      DV   CWRES
#> 1  12  5.18  534.31 -8.1999
#> 2  39 51.45 4853.00  7.4235
#> 3  39  6.30 3925.00  6.3175
#> 4  13 28.35 1532.10 -5.0007
#> 5  18 10.40  104.09 -4.3844
#> 6  39 52.67 2968.50  4.0617
#> 7  73 29.73 2311.90  4.0013
#> 8  54 27.31 7712.70  3.9900
#> 9  90 60.65  700.97  3.3856
#> 10 79 54.71 1790.30  3.4172
#> 11 39  4.55 4493.50  3.1948
#> 12 73 56.59 1257.80  3.1732
#> 13 81  3.17 5075.90  3.3176
#> 14  3  7.01 1641.80  3.1707
#> 15 52  7.37 1222.60  3.1883
#> 16 28 52.58 2957.40  2.9034
#> 17 63 51.66 3498.40  2.9089
#> 18 72 51.67 3474.70  2.7608
#> 19 16 56.82 1220.30  2.7564
#> 20  8 59.63  644.55  2.7565
```

------------------------------------------------------------------------

## 3. Removal metrics

``` r

plot_removal_metrics(results, metric = "pOFV", verbose = FALSE)
```

![Change in pseudo-OFV per iteration. A sharp early drop followed by a
plateau suggests meaningful early removals with diminishing
returns.](busulfan_pt2_modelbased_cleaning_files/figure-html/metrics-pofv-1.png)

Change in pseudo-OFV per iteration. A sharp early drop followed by a
plateau suggests meaningful early removals with diminishing returns.

``` r

plot_removal_metrics(results, metric = "thetas", verbose = FALSE)
```

![THETA estimates normalised to final value. Free axes per panel;
phantom reference lines at 0.75 and 1.5 ensure consistent scale even
when individual parameters are
stable.](busulfan_pt2_modelbased_cleaning_files/figure-html/metrics-thetas-1.png)

THETA estimates normalised to final value. Free axes per panel; phantom
reference lines at 0.75 and 1.5 ensure consistent scale even when
individual parameters are stable.

``` r

plot_removal_metrics(results, metric = "rmse", verbose = FALSE)
```

![nRMSE across
iterations.](busulfan_pt2_modelbased_cleaning_files/figure-html/metrics-rmse-1.png)

nRMSE across iterations.

------------------------------------------------------------------------

## 4. Error parameter trajectories

Proportional and additive residual error parameters are auto-detected
from the `param_labels` element of the results (parsed from `$THETA`
comments in the model file). For this model, `THETA7` is proportional
error and `THETA8` is additive. The detection source is printed when
`quiet = FALSE` in
[`generate_report()`](https://insightrx.github.io/irxclean/reference/generate_report.md).

``` r

detected <- irxclean:::.find_error_params(
  results$param_labels, results$par, verbose = TRUE
)
#> irxclean: Auto-detected error parameters (source: THETA label):
#>   Proportional : "THETA7" (label: "PROP error")
#>   Additive     : "THETA8" (label: "ADD error")
#>   Override via prop_error_param / add_error_param (e.g. "THETA(7)" or "SIGMA(1)").
```

``` r

irxclean:::.plot_param_trajectory(
  results,
  detected$prop_col,
  title = "Proportional residual error (THETA7) across iterations"
)
```

![Proportional residual error (THETA7) stability across
iterations.](busulfan_pt2_modelbased_cleaning_files/figure-html/traj-prop-1.png)

Proportional residual error (THETA7) stability across iterations.

``` r

irxclean:::.plot_param_trajectory(
  results,
  detected$add_col,
  title = "Additive residual error (THETA8) across iterations"
)
```

![Additive residual error (THETA8) stability across
iterations.](busulfan_pt2_modelbased_cleaning_files/figure-html/traj-add-1.png)

Additive residual error (THETA8) stability across iterations.

------------------------------------------------------------------------

## 5. Model stability

``` r

stability <- check_model_stability(results, threshold_pct = 20)
#> 
#> irxclean :: Model stability check (threshold: 20%)
#> --------------------------------------------------
#>   Parameters assessed : 17  (fixed excluded)
#>     stable     : 4
#>     drifted    : 3
#>     unstable   : 10
#>   Unstable: MAT-MAG, ADD error, DROP, SHAPE, Q, V2, IIV V1, IOV CL, IOV V1, IIV V2
#>   nRMSE trend: decreasing
print(stability)
#> irxclean model stability (threshold: 20%)
#>   stable     : 4 parameter(s)
#>   drifted    : 3 parameter(s)
#>   unstable   : 10 parameter(s)
#> 
#> Parameter summary:
#>        PARAMETER RSE_PCT N_REVERSALS CV_LATE_PCT DRIFT_PCT   STATUS
#>            TH_CL   37.32          12       8.065  5.08e+01  drifted
#>             TH_V   34.43          10       0.326  1.92e+01  drifted
#>          MAT-MAG 7870.87           9      16.079  1.36e+01 unstable
#>            K_MAT   37.68           7       3.686  5.21e+00   stable
#>       PROP error    6.10           4       0.797  1.97e+01  drifted
#>        ADD error   78.03          11       0.270  4.77e-01 unstable
#>             DROP   88.68          12      13.318  5.86e+02 unstable
#>            SHAPE  443.64          10      15.904  1.07e+04 unstable
#>     allo ffm exp    4.61           8       0.984  6.59e+00   stable
#>  Sex effect on V    1.27          11       0.149  2.02e+00   stable
#>                Q  119.05           8       4.428  7.14e+01 unstable
#>               V2   98.14          10       2.221  5.37e+01 unstable
#>           IIV CL    9.84           9       0.268  1.96e+00   stable
#>           IIV V1  137.25          10       0.721  2.30e+01 unstable
#>           IOV CL  135.55           9       0.000  9.78e+01 unstable
#>           IOV V1  370.66           3       0.000  0.00e+00 unstable
#>           IIV V2   28.55           9      10.694  3.16e+01 unstable
#> 
#> nRMSE trend: decreasing
```

``` r

flag_tbl <- stability_flag_table(stability)
if (nrow(flag_tbl) > 0) knitr::kable(flag_tbl, digits = 2) else
  cat("All parameters stable.\n")
```

|  | PARAMETER | MEAN | FIRST | LAST | RSE_PCT | N_REVERSALS | CV_LATE_PCT | DRIFT_PCT | STATUS |
|:---|:---|---:|---:|---:|---:|---:|---:|---:|:---|
| 3 | MAT-MAG | 0.01 | 0.10 | 0.12 | 7870.87 | 9 | 16.08 | 13.62 | unstable |
| 6 | ADD error | 26.36 | 32.77 | 32.62 | 78.03 | 11 | 0.27 | 0.48 | unstable |
| 7 | DROP | 0.85 | 0.23 | 1.57 | 88.68 | 12 | 13.32 | 586.12 | unstable |
| 8 | SHAPE | -0.01 | 0.00 | 0.00 | 443.64 | 10 | 15.90 | 10737.04 | unstable |
| 11 | Q | 6.42 | 10.92 | 3.13 | 119.05 | 8 | 4.43 | 71.35 | unstable |
| 12 | V2 | 7.01 | 9.07 | 4.20 | 98.14 | 10 | 2.22 | 53.65 | unstable |
| 14 | IIV V1 | 0.10 | 0.07 | 0.05 | 137.25 | 10 | 0.72 | 22.99 | unstable |
| 15 | IOV CL | 0.00 | 0.00 | 0.00 | 135.55 | 9 | 0.00 | 97.81 | unstable |
| 16 | IOV V1 | 0.03 | 0.00 | 0.00 | 370.66 | 3 | 0.00 | 0.00 | unstable |
| 17 | IIV V2 | 0.06 | 0.11 | 0.07 | 28.55 | 9 | 10.69 | 31.63 | unstable |

------------------------------------------------------------------------

## 6. Demographic comparison

``` r

plot_demographic_comparison(
  data             = dat,
  results          = results,
  continuous_cols  = c("AGE", "WT", "HT"),
  categorical_cols = "SEX"
)
#> TAD column not found - computed from dose records (EVID == 1 / AMT > 0).
#> irxclean demographic comparison (16 subjects with removed obs)
#> 
#> Statistical tests:
#>  covariate              test   p_value significant
#>        AGE Wilcoxon rank-sum 0.4351332       FALSE
#>         WT Wilcoxon rank-sum 0.2774897       FALSE
#>         HT Wilcoxon rank-sum 0.2691470       FALSE
#>        SEX       Chi-squared 0.4710967       FALSE
```

------------------------------------------------------------------------

## 7. Apply final exclusions

After reviewing the ranked flagged observations, apply user-approved
exclusions:

``` r

cleaned <- remove_observations_from_data(
  data    = dat,
  results = results,
  n       = 10  # remove top 10ranked observations
)

write.csv(cleaned, "busulfan_cleaned.csv", row.names = FALSE)
cat("Cleaned dataset:", nrow(cleaned), "rows\n")
```

------------------------------------------------------------------------

## 8. Generate HTML report

[`generate_report()`](https://insightrx.github.io/irxclean/reference/generate_report.md)
compiles all findings into a structured, self-contained HTML report.
Proportional and additive error parameters are auto-detected from
`results$param_labels`; override explicitly if needed using either the
column name (`"THETA7"`) or NONMEM index notation (`"THETA(7)"`).

``` r

generate_report(
  results              = results,
  data                 = dat,
  # Auto-detected; shown here for illustration -- override if needed:
  prop_error_param     = "THETA(7)",
  add_error_param      = "THETA(8)",
  continuous_cov_cols  = c("AGE", "WT", "HT"),
  categorical_cov_cols = "SEX",
  title                = "Busulfan PK -- Model-Informed Cleaning Report",
  output_file          = "busulfan_irxclean_report.html",
  quiet                = FALSE   # prints auto-detection summary
)
```

------------------------------------------------------------------------

## Session information

``` r

sessionInfo()
#> R version 4.6.0 (2026-04-24)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.4 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] irxclean_0.1.0.9000
#> 
#> loaded via a namespace (and not attached):
#>  [1] Matrix_1.7-5       gtable_0.3.6       jsonlite_2.0.0     dplyr_1.2.1       
#>  [5] compiler_4.6.0     tidyselect_1.2.1   stringr_1.6.0      tidyr_1.3.2       
#>  [9] jquerylib_0.1.4    splines_4.6.0      systemfonts_1.3.2  scales_1.4.0      
#> [13] textshaping_1.0.5  yaml_2.3.12        fastmap_1.2.0      lattice_0.22-9    
#> [17] ggplot2_4.0.3      R6_2.6.1           labeling_0.4.3     vpc_1.2.4         
#> [21] generics_0.1.4     patchwork_1.3.2    knitr_1.51         tibble_3.3.1      
#> [25] desc_1.4.3         bslib_0.10.0       pillar_1.11.1      RColorBrewer_1.1-3
#> [29] rlang_1.2.0        stringi_1.8.7      cachem_1.1.0       xfun_0.57         
#> [33] fs_2.1.0           sass_0.4.10        S7_0.2.2           cli_3.6.6         
#> [37] mgcv_1.9-4         withr_3.0.2        pkgdown_2.2.0      magrittr_2.0.5    
#> [41] digest_0.6.39      grid_4.6.0         nlme_3.1-169       lifecycle_1.0.5   
#> [45] vctrs_0.7.3        evaluate_1.0.5     glue_1.8.1         farver_2.1.2      
#> [49] ragg_1.5.2         purrr_1.2.2        rmarkdown_2.31     tools_4.6.0       
#> [53] pkgconfig_2.0.3    htmltools_0.5.9
```
