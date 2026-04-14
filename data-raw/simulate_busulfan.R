# =============================================================================
# simulate_busulfan.R
# Simulate a busulfan PK dataset for the irxclean vignette.
#
# Model  : pk_busulfan_ucsf_2025 (2-cmt IV, allometric FFM + maturation,
#          time-varying clearance, IOV on CL and V)
# Tools  : PKPDsim (simulation), PKPDmap (MAP Bayesian dose adaptation)
# Dosing : 4 x Q24 h IV infusions (3 h), TDM at 4 random times per dose
#          (windows: 10-45 min, 1-3 h, 3-6 h, 6-12 h post-infusion end)
#          for each of doses 1-3.  Doses 2-4 are adapted to target
#          cAUC = 90 mg*h/L
#          (= 90 000 ng/mL*h cumulative) via individual MAP Bayesian AUC
#          prediction using TDM collected after the previous dose.
#
# Output : inst/extdata/busulfan_sim.csv  (NONMEM-formatted, 100 patients)
#
# Prerequisites:
#   - PKPDsim (>= 1.5.0), PKPDmap (>= 1.1.6), MASS installed
#   - pkbusulfanucsf2025 model package installed; if not present, run:
#       PKPDsim::model_from_api(
#         model = "pk_busulfan_ucsf_2025",
#         url   = system.file("models", "pk_busulfan_ucsf_2025.json5",
#                             package = "irxclean"),
#         to_package = TRUE, force = TRUE, install = TRUE
#       )
#
# Run from the package root:
#   source("data-raw/simulate_busulfan.R")
# =============================================================================

library(PKPDsim)
library(PKPDmap)
library(pkbusulfanucsf2025)
library(MASS)
library(dplyr)

set.seed(2025)
outfile <- "inst/extdata/busulfan_sim.csv"

# -----------------------------------------------------------------------------
# 1. Load model components from the installed package
# -----------------------------------------------------------------------------
mod  <- model()
pars <- parameters()      # population parameters (structural + IOV kappas)
om   <- omega_matrix()    # lower triangle of IIV omega (3 params: CL, V, V2)
iov_ <- iov()             # IOV spec: CL CV=11.6%, V CV=11.2%, 5 occasions
ruv_ <- ruv()             # prop=0.069, add=28.6 ng/mL

om_full <- PKPDsim::triangle_to_full(om)   # full 3x3 IIV omega matrix

# Parameters always fixed in MAP (structural + IOV kappas for occasions 2-5)
# Only occasion-1 kappas (kappa_CL_1, kappa_V_1) are estimated after dose 1.
fixed_structural <- names(pars)[!grepl("^kappa_", names(pars))]
fixed_iov_2to5   <- c(paste0("kappa_CL_", 2:5), paste0("kappa_V_", 2:5))
map_fixed   <- c(fixed_structural, fixed_iov_2to5)
map_as_eta  <- c("kappa_CL_1", "kappa_V_1")

message("Model loaded: pkbusulfanucsf2025")
message("  Parameters: ", paste(names(pars), collapse = ", "))

# -----------------------------------------------------------------------------
# 2. Population parameters & variability
# -----------------------------------------------------------------------------
target_cauc <- 90000   # ng/mL*h cumulative (4 doses)
t_inf       <- 3       # infusion duration (h)

# TDM sampling windows: (lower, upper) hours relative to dose start (TAD).
# Each window corresponds to the requested sampling time post-infusion end
# (infusion ends at t_inf = 3 h):
#   window 1:  10-45 min post infusion  -> TAD 3.17 - 3.75 h
#   window 2:   1-3 h after infusion    -> TAD 4.00 - 6.00 h
#   window 3:   3-6 h after infusion    -> TAD 6.00 - 9.00 h
#   window 4:  6-12 h after infusion    -> TAD 9.00 - 15.00 h
tdm_windows <- list(
  c(t_inf + 10/60, t_inf + 45/60),
  c(t_inf + 1,     t_inf + 3),
  c(t_inf + 3,     t_inf + 6),
  c(t_inf + 6,     t_inf + 12)
)

# IOV variance (CV -> variance on log scale)
iov_sd_cl <- sqrt(log(1 + iov_$cv$CL^2))
iov_sd_v  <- sqrt(log(1 + iov_$cv$V^2))

# -----------------------------------------------------------------------------
# 3. Helper functions
# -----------------------------------------------------------------------------

# Draw one random TAD per window (rounded to 2 dp, ~36 s precision)
.sample_tad <- function() {
  vapply(tdm_windows, function(w) round(runif(1, w[1], w[2]), 2), numeric(1))
}

# Build covariate list for PKPDsim
.covs <- function(age, wt, ht, sex) {
  list(
    AGE = new_covariate(age),
    WT  = new_covariate(wt),
    HT  = new_covariate(ht),
    SEX = new_covariate(sex)
  )
}

# Predict AUC for one dose using current MAP-estimated parameters
.predict_dose_auc <- function(dose, t_start, params_map, cov_i, A_init) {
  reg <- new_regimen(
    amt      = dose,
    n        = 1,
    interval = 24,
    type     = "infusion",
    t_inf    = t_inf
  )
  s <- sim(
    ode        = mod,
    parameters = params_map,
    regimen    = reg,
    covariates = cov_i,
    t_obs      = 24,         # end of dose interval
    A_init     = A_init,
    only_obs   = FALSE
  )
  list(
    auc   = max(s$y[s$comp == 3]),   # cumulative AUC in comp 3 (ng/mL*h)
    A_end = s$y[s$comp != "obs" & s$t == max(s$t)]
  )
}

# Simulate TDM observations for one dose, return conc (ng/mL) and state.
# obs_tad: TAD times (relative to dose start) at which TDM is collected.
# t_obs passed to sim() must be RELATIVE to the regimen start (TAD), not
# absolute clock time.  We also simulate to t=24 (end of dosing interval)
# so that A_end reflects the true carry-over state for the next dose.
.sim_dose_tdm <- function(dose, t_start, params_true, cov_i, A_init, obs_tad) {
  reg <- new_regimen(
    amt      = dose,
    n        = 1,
    interval = 24,
    type     = "infusion",
    t_inf    = t_inf
  )
  # Simulate at TDM TADs + end of dosing interval (for A_end handoff)
  t_sim <- sort(unique(c(obs_tad, 24)))
  s <- sim(
    ode        = mod,
    parameters = params_true,
    regimen    = reg,
    covariates = cov_i,
    t_obs      = t_sim,   # relative times (TAD); NOT absolute clock times
    A_init     = A_init,
    only_obs   = FALSE
  )
  # TDM concentrations at TAD time points
  conc  <- s$y[s$comp == "obs" & s$t %in% obs_tad]
  # Compartment state at end of dosing interval (t=24h relative to dose)
  A_end <- s$y[s$comp != "obs" & s$t == 24]
  list(conc = conc, t = t_start + obs_tad, A_end = A_end)
}

# Run MAP after one dose and adapt the next dose
.map_adapt_dose <- function(dv_obs, t_obs, dose_used, cov_i, occ_kappas) {
  reg_map <- new_regimen(
    amt      = dose_used,
    n        = 1,
    interval = 24,
    type     = "infusion",
    t_inf    = t_inf
  )
  # Set occasion kappas in population parameters
  pars_map <- pars
  pars_map[names(occ_kappas)] <- occ_kappas

  fit <- get_map_estimates(
    model      = mod,
    data       = data.frame(t = t_obs, y = dv_obs),
    parameters = pars_map,
    covariates = cov_i,
    omega      = om_full,
    error      = ruv_,
    regimen    = reg_map,
    fixed      = map_fixed,
    as_eta     = map_as_eta,
    verbose    = FALSE
  )
  fit$parameters
}

# -----------------------------------------------------------------------------
# 4. Generate patient covariates
# -----------------------------------------------------------------------------
n       <- 100
n_ped   <- 80
n_adult <- 20

# Pediatric (1-18 yr)
age_ped  <- round(runif(n_ped,  1, 18), 1)
sex_ped  <- rbinom(n_ped, 1, 0.40)
# Weight: approx 3*age+6 kg for children, capped at ~65 kg
wt_ped_mean <- pmin(3 * age_ped + 6, 65)
wt_ped   <- pmax(rlnorm(n_ped, log(wt_ped_mean), 0.15), 5)
ht_ped   <- pmax(75 + 5.5 * age_ped + rnorm(n_ped, 0, 6), 60)

# Adult (18-70 yr)
age_adult <- round(runif(n_adult, 18, 70), 0)
sex_adult <- rbinom(n_adult, 1, 0.40)
wt_adult  <- pmax(rnorm(n_adult, 72, 16), 40)
ht_adult  <- pmax(rnorm(n_adult, 170, 10), 150)

covariates <- data.frame(
  ID  = seq_len(n),
  AGE = c(age_ped,  age_adult),
  WT  = round(c(wt_ped,  wt_adult), 1),
  HT  = round(c(ht_ped,  ht_adult), 0),
  SEX = c(sex_ped,  sex_adult)    # 0 = female, 1 = male
)

# -----------------------------------------------------------------------------
# 5. Sample IIV from omega matrix
# -----------------------------------------------------------------------------
etas_iiv <- mvrnorm(n, mu = c(0, 0, 0), Sigma = om_full)
colnames(etas_iiv) <- c("CL", "V", "V2")

# Sample IOV (5 occasions) for each patient
# IOV kappas are independent per occasion, zero-mean normal on log scale
iov_kappas <- lapply(seq_len(n), function(i) {
  list(
    kappa_CL_1 = rnorm(1, 0, iov_sd_cl),
    kappa_CL_2 = rnorm(1, 0, iov_sd_cl),
    kappa_CL_3 = rnorm(1, 0, iov_sd_cl),
    kappa_CL_4 = rnorm(1, 0, iov_sd_cl),
    kappa_CL_5 = rnorm(1, 0, iov_sd_cl),
    kappa_V_1  = rnorm(1, 0, iov_sd_v),
    kappa_V_2  = rnorm(1, 0, iov_sd_v),
    kappa_V_3  = rnorm(1, 0, iov_sd_v),
    kappa_V_4  = rnorm(1, 0, iov_sd_v),
    kappa_V_5  = rnorm(1, 0, iov_sd_v)
  )
})

# -----------------------------------------------------------------------------
# 6. Simulate PK for each patient
# -----------------------------------------------------------------------------
message("Simulating ", n, " patients...")
all_rows <- vector("list", n)

for (i in seq_len(n)) {
  cv  <- covariates[i, ]
  cov <- .covs(cv$AGE, cv$WT, cv$HT, cv$SEX)

  # True individual parameters (IIV + IOV)
  params_true <- pars
  params_true[["kappa_CL_1"]] <- etas_iiv[i, "CL"] + iov_kappas[[i]][["kappa_CL_1"]]
  params_true[["kappa_CL_2"]] <- etas_iiv[i, "CL"] + iov_kappas[[i]][["kappa_CL_2"]]
  params_true[["kappa_CL_3"]] <- etas_iiv[i, "CL"] + iov_kappas[[i]][["kappa_CL_3"]]
  params_true[["kappa_CL_4"]] <- etas_iiv[i, "CL"] + iov_kappas[[i]][["kappa_CL_4"]]
  params_true[["kappa_CL_5"]] <- etas_iiv[i, "CL"] + iov_kappas[[i]][["kappa_CL_5"]]
  params_true[["kappa_V_1"]]  <- etas_iiv[i, "V"]  + iov_kappas[[i]][["kappa_V_1"]]
  params_true[["kappa_V_2"]]  <- etas_iiv[i, "V"]  + iov_kappas[[i]][["kappa_V_2"]]
  params_true[["kappa_V_3"]]  <- etas_iiv[i, "V"]  + iov_kappas[[i]][["kappa_V_3"]]
  params_true[["kappa_V_4"]]  <- etas_iiv[i, "V"]  + iov_kappas[[i]][["kappa_V_4"]]
  params_true[["kappa_V_5"]]  <- etas_iiv[i, "V"]  + iov_kappas[[i]][["kappa_V_5"]]

  # ------------------------------------------------------------------
  # Dose 1: weight-based empirical starting dose (3.2 mg/kg, round to 5 mg)
  dose1  <- pmax(round(3.2 * cv$WT / 5) * 5, 5)
  A0     <- rep(0, attr(mod, "size"))

  tad1   <- .sample_tad()
  res1   <- .sim_dose_tdm(dose1, 0, params_true, cov, A0, tad1)
  prop1  <- rnorm(length(tad1), 0, ruv_$prop)
  add1   <- rnorm(length(tad1), 0, ruv_$add)
  dv1    <- pmax(res1$conc * (1 + prop1) + add1, 0.01)

  # ------------------------------------------------------------------
  # MAP after dose 1 -> adapt dose 2
  # Estimate kappa_CL_1 and kappa_V_1 from dose-1 TDM
  # (occasion-2+ kappas remain at population value 0)
  occ1_pars <- list(kappa_CL_1=0, kappa_CL_2=0, kappa_CL_3=0, kappa_CL_4=0, kappa_CL_5=0,
                    kappa_V_1=0,  kappa_V_2=0,  kappa_V_3=0,  kappa_V_4=0,  kappa_V_5=0)
  params_map1 <- .map_adapt_dose(dv1, tad1, dose1, cov, occ1_pars)

  # Predict dose-1 AUC with MAP parameters (A_init=0, predict over 24h)
  pred1      <- .predict_dose_auc(dose1, 0, params_map1, cov, A0)
  auc1_pred  <- pred1$auc
  dose2      <- pmax(round((target_cauc / 4) * dose1 / auc1_pred / 5) * 5, 5)

  # ------------------------------------------------------------------
  # Dose 2 (t=24)
  tad2  <- .sample_tad()
  res2  <- .sim_dose_tdm(dose2, 24, params_true, cov, res1$A_end, tad2)
  prop2 <- rnorm(length(tad2), 0, ruv_$prop)
  add2  <- rnorm(length(tad2), 0, ruv_$add)
  dv2   <- pmax(res2$conc * (1 + prop2) + add2, 0.01)

  # MAP after doses 1-2 -> adapt dose 3
  # Update occ-2 kappas; kappa_CL_2 = MAP kappa_CL_1 (IIV estimate carries over)
  occ2_pars <- occ1_pars
  occ2_pars[["kappa_CL_1"]] <- params_map1[["kappa_CL_1"]]
  occ2_pars[["kappa_V_1"]]  <- params_map1[["kappa_V_1"]]
  params_map2 <- .map_adapt_dose(dv2, tad2, dose2, cov, occ2_pars)

  pred2     <- .predict_dose_auc(dose2, 0, params_map2, cov, A0)
  dose3     <- pmax(round((target_cauc / 4) * dose2 / pred2$auc / 5) * 5, 5)

  # ------------------------------------------------------------------
  # Dose 3 (t=48)
  tad3  <- .sample_tad()
  res3  <- .sim_dose_tdm(dose3, 48, params_true, cov, res2$A_end, tad3)
  prop3 <- rnorm(length(tad3), 0, ruv_$prop)
  add3  <- rnorm(length(tad3), 0, ruv_$add)
  dv3   <- pmax(res3$conc * (1 + prop3) + add3, 0.01)

  # MAP after doses 1-3 -> adapt dose 4
  occ3_pars <- occ2_pars
  occ3_pars[["kappa_CL_1"]] <- params_map2[["kappa_CL_1"]]
  occ3_pars[["kappa_V_1"]]  <- params_map2[["kappa_V_1"]]
  params_map3 <- .map_adapt_dose(dv3, tad3, dose3, cov, occ3_pars)

  pred3  <- .predict_dose_auc(dose3, 0, params_map3, cov, A0)
  dose4  <- pmax(round((target_cauc / 4) * dose3 / pred3$auc / 5) * 5, 5)

  # ------------------------------------------------------------------
  # Build NONMEM-format rows for this patient
  build_rows <- function(dose_k, time_dose, dv_k, times_obs) {
    dose_row <- data.frame(
      ID = i, TIME = time_dose, AMT = dose_k, RATE = dose_k / t_inf,
      DV = 0, MDV = 1, EVID = 1, CMT = 1,
      AGE = cv$AGE, WT = cv$WT, HT = cv$HT, SEX = cv$SEX,
      FLAG = 0L, FLAG_LABEL = ""
    )
    obs_rows <- data.frame(
      ID = i, TIME = times_obs, AMT = 0, RATE = 0,
      DV = dv_k, MDV = 0, EVID = 0, CMT = 1,
      AGE = cv$AGE, WT = cv$WT, HT = cv$HT, SEX = cv$SEX,
      FLAG = 0L, FLAG_LABEL = ""
    )
    rbind(dose_row, obs_rows)
  }

  all_rows[[i]] <- rbind(
    build_rows(dose1, 0,  dv1, tad1),
    build_rows(dose2, 24, dv2, tad2 + 24),
    build_rows(dose3, 48, dv3, tad3 + 48),
    data.frame(
      ID = i, TIME = 72, AMT = dose4, RATE = dose4 / t_inf, DV = 0,
      MDV = 1, EVID = 1, CMT = 1,
      AGE = cv$AGE, WT = cv$WT, HT = cv$HT, SEX = cv$SEX,
      FLAG = 0L, FLAG_LABEL = ""
    )
  )
}

dat <- dplyr::bind_rows(all_rows) |>
  dplyr::arrange(ID, TIME)

message("Simulation complete: ", nrow(dat), " rows, ",
        sum(dat$EVID == 0), " observations across ", n, " patients.")

# -----------------------------------------------------------------------------
# 7. Introduce intentional errors
# -----------------------------------------------------------------------------
obs_idx  <- which(dat$EVID == 0 & dat$MDV == 0)

# -- 7a. Weight decimal error (1 patient) ------------------------------------
wt_candidates <- which(covariates$AGE < 15 & covariates$WT >= 20 &
                         covariates$WT < 45)
pt_wt_err <- sample(wt_candidates, 1)
wt_true   <- covariates$WT[pt_wt_err]
wt_error  <- round(wt_true / 10, 1)   # e.g. 28.4 -> 2.8
dat$WT[dat$ID == pt_wt_err] <- wt_error
dat$FLAG[dat$ID == pt_wt_err] <- bitwOr(dat$FLAG[dat$ID == pt_wt_err], 1L)
dat$FLAG_LABEL[dat$ID == pt_wt_err] <- ifelse(
  nchar(dat$FLAG_LABEL[dat$ID == pt_wt_err]) > 0,
  paste0(dat$FLAG_LABEL[dat$ID == pt_wt_err], "|covariate_weight"),
  "covariate_weight"
)
message(sprintf("  WT error: patient %d, WT %.1f -> %.1f",
                pt_wt_err, wt_true, wt_error))

# -- 7b. Concentration IQR outliers (2 observations) -------------------------
# Target observations in the later two sampling windows (TAD >= 6 h), one
# from each window (windows 3 and 4), from different patients.
dat_tad <- dat$TIME %% 24

set.seed(7)
iqr_err_rows <- integer(0)
for (win_bounds in list(c(6, 9), c(9, 15))) {
  cands <- obs_idx[
    dat_tad[obs_idx] >= win_bounds[1] &
    dat_tad[obs_idx] <  win_bounds[2] &
    !(dat$ID[obs_idx] %in% c(pt_wt_err, dat$ID[iqr_err_rows]))
  ]
  if (length(cands) > 0) {
    pick     <- sample(cands, 1)
    win_dv   <- dat$DV[obs_idx[dat_tad[obs_idx] >= win_bounds[1] &
                                dat_tad[obs_idx] <  win_bounds[2]]]
    upper_thr <- quantile(win_dv, 0.75) + 3 * IQR(win_dv)
    dat$DV[pick] <- upper_thr * runif(1, 1.5, 2.5)
    iqr_err_rows <- c(iqr_err_rows, pick)
  }
}
dat$FLAG[iqr_err_rows] <- bitwOr(dat$FLAG[iqr_err_rows], 2L)
dat$FLAG_LABEL[iqr_err_rows] <- ifelse(
  nchar(dat$FLAG_LABEL[iqr_err_rows]) > 0,
  paste0(dat$FLAG_LABEL[iqr_err_rows], "|conc_iqr_outlier"),
  "conc_iqr_outlier"
)
message(sprintf("  IQR outliers: rows %s", paste(iqr_err_rows, collapse = ", ")))

# -- 7c. Timing errors (10 observations) -------------------------------------
timing_offsets_hr <- c(-5, 5, -10, 10, -30, 30, -60, 30, -60, 60) / 60
already_flagged   <- c(obs_idx[dat$ID[obs_idx] == pt_wt_err], iqr_err_rows)
timing_pool       <- setdiff(obs_idx, already_flagged)
timing_sel        <- sample(timing_pool, 10)
for (k in seq_along(timing_sel)) {
  row_k  <- timing_sel[k]
  offset <- timing_offsets_hr[k]
  dat$TIME[row_k] <- dat$TIME[row_k] + offset
  label_k <- sprintf("timing_%dmin", round(abs(offset * 60)))
  dat$FLAG[row_k]       <- bitwOr(dat$FLAG[row_k], 4L)
  dat$FLAG_LABEL[row_k] <- ifelse(
    nchar(dat$FLAG_LABEL[row_k]) > 0,
    paste0(dat$FLAG_LABEL[row_k], "|", label_k),
    label_k
  )
}
message(sprintf("  Timing errors: %d observations adjusted", length(timing_sel)))

# -- 7d. Concentration magnitude errors (10 observations) --------------------
mult_factors     <- c(0.40, 0.25, 0.60, 0.80, 0.95, 1.05, 1.20, 1.40, 1.75, 1.60)
already_flagged2 <- c(already_flagged, timing_sel)
mag_pool         <- setdiff(obs_idx, already_flagged2)
mag_sel          <- sample(mag_pool, 10)
for (k in seq_along(mag_sel)) {
  row_k  <- mag_sel[k]
  mult   <- mult_factors[k]
  dat$DV[row_k]         <- dat$DV[row_k] * mult
  label_k               <- sprintf("magnitude_%.2fx", mult)
  dat$FLAG[row_k]       <- bitwOr(dat$FLAG[row_k], 8L)
  dat$FLAG_LABEL[row_k] <- ifelse(
    nchar(dat$FLAG_LABEL[row_k]) > 0,
    paste0(dat$FLAG_LABEL[row_k], "|", label_k),
    label_k
  )
}
message(sprintf("  Magnitude errors: %d observations modified", length(mag_sel)))

# -----------------------------------------------------------------------------
# 8. Final formatting and export
# -----------------------------------------------------------------------------
dat$DV  <- round(dat$DV, 4)
dat$WT  <- round(dat$WT, 1)
dat$AMT <- round(dat$AMT, 1)

n_cov_err <- length(unique(dat$ID[dat$ID == pt_wt_err]))
message("\nError summary:")
message(sprintf("  Covariate (WT)  : %d patient  (%d rows affected)",
                n_cov_err, sum(dat$ID == pt_wt_err)))
message(sprintf("  IQR outliers    : %d observations", length(iqr_err_rows)))
message(sprintf("  Timing errors   : %d observations", length(timing_sel)))
message(sprintf("  Magnitude errors: %d observations", length(mag_sel)))
message(sprintf("  Rows with FLAG>0: %d (obs only)",
                sum(dat$FLAG > 0 & dat$EVID == 0)))

write.csv(dat, outfile, row.names = FALSE, quote = FALSE)
message("\nDataset saved: ", outfile)
message("Rows: ", nrow(dat), "  |  Patients: ", n,
        "  |  Observations: ", sum(dat$EVID == 0))
