#' Simulated busulfan pharmacokinetic dataset
#'
#' A NONMEM-formatted dataset simulated from a two-compartment IV busulfan
#' population PK model (ADVAN3 TRANS4) with allometric fat-free mass scaling,
#' maturation function, time-varying clearance, and inter-individual /
#' inter-occasion variability.  100 patients (80 paediatric, 20 adult) each
#' received 4 MAP-adapted doses with 4 therapeutic drug monitoring samples per
#' dose (1,200 observations total).
#'
#' Intentional errors were introduced to support demonstration and testing of
#' \code{irxclean} workflows:
#' \describe{
#'   \item{FLAG bit 1}{Covariate weight decimal-place shift (1 patient)}
#'   \item{FLAG bit 2}{Concentration IQR outlier >3*IQR above Q3 per TAD bin (2 observations)}
#'   \item{FLAG bit 4}{Timing error +/- 5-60 minutes from true time (10 observations)}
#'   \item{FLAG bit 8}{Concentration magnitude error (multipliers 0.25-1.75x; 10 observations)}
#' }
#'
#' @format A data frame with 1,600 rows and 14 columns:
#' \describe{
#'   \item{ID}{Subject identifier (integer)}
#'   \item{TIME}{Nominal time since first dose (hours)}
#'   \item{AMT}{Dose amount (mg); 0 for observation rows}
#'   \item{RATE}{Infusion rate (mg/h); 0 for observation rows}
#'   \item{DV}{Observed busulfan concentration (ng/mL); 0 for dose rows}
#'   \item{MDV}{Missing dependent variable flag (1 = dose/missing, 0 = observation)}
#'   \item{EVID}{Event ID (1 = dose, 0 = observation)}
#'   \item{CMT}{Compartment number}
#'   \item{AGE}{Age (years)}
#'   \item{WT}{Body weight (kg)}
#'   \item{HT}{Height (cm)}
#'   \item{SEX}{Sex (1 = male, 2 = female)}
#'   \item{FLAG}{Bitmask indicating intentionally introduced error type(s) (0 = clean)}
#'   \item{FLAG_LABEL}{Character description of error type(s); empty string for clean records}
#' }
#'
#' @source Simulated using PKPDsim with the \code{pkbusulfanucsf2025} model.
#'   See \code{data-raw/simulate_busulfan.R} for the full simulation script.
#'
#' @keywords datasets
"busulfan_sim"
