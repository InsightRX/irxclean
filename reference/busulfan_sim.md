# Simulated busulfan pharmacokinetic dataset

A NONMEM-formatted dataset simulated from a two-compartment IV busulfan
population PK model (ADVAN3 TRANS4) with allometric fat-free mass
scaling, maturation function, time-varying clearance, and
inter-individual / inter-occasion variability. 100 patients (80
paediatric, 20 adult) each received 4 MAP-adapted doses with 4
therapeutic drug monitoring samples per dose (1,200 observations total).

## Usage

``` r
busulfan_sim
```

## Format

A data frame with 1,600 rows and 14 columns:

- ID:

  Subject identifier (integer)

- TIME:

  Nominal time since first dose (hours)

- AMT:

  Dose amount (mg); 0 for observation rows

- RATE:

  Infusion rate (mg/h); 0 for observation rows

- DV:

  Observed busulfan concentration (ng/mL); 0 for dose rows

- MDV:

  Missing dependent variable flag (1 = dose/missing, 0 = observation)

- EVID:

  Event ID (1 = dose, 0 = observation)

- CMT:

  Compartment number

- AGE:

  Age (years)

- WT:

  Body weight (kg)

- HT:

  Height (cm)

- SEX:

  Sex (1 = male, 2 = female)

- FLAG:

  Bitmask indicating intentionally introduced error type(s) (0 = clean)

- FLAG_LABEL:

  Character description of error type(s); empty string for clean records

## Source

Simulated using PKPDsim with the `pkbusulfanucsf2025` model. See
`data-raw/simulate_busulfan.R` for the full simulation script.

## Details

Intentional errors were introduced to support demonstration and testing
of `irxclean` workflows:

- FLAG bit 1:

  Covariate weight decimal-place shift (1 patient)

- FLAG bit 2:

  Concentration IQR outlier \>3\*IQR above Q3 per TAD bin (2
  observations)

- FLAG bit 4:

  Timing error +/- 5-60 minutes from true time (10 observations)

- FLAG bit 8:

  Concentration magnitude error (multipliers 0.25-1.75x; 10
  observations)
