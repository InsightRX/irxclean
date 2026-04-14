$PROB   pk_busulfan_ucsf_2025
;       Two-compartment IV busulfan model (UCSF/InsightRX 2025)
;       CL : allometric FFM scaling, maturation function, time-varying
;       IIV: CL, V1, V2   (3x3 BLOCK omega)
;       IOV: CL, V1       (5 occasions, 24-h bins)
;       RUV: proportional + additive
;
;       Dataset: busulfan_sim.csv  (cols: ID TIME AMT DV MDV EVID CMT
;                                         AGE WT HT SEX FLAG FLAG_LABEL)
;       Note : FLAG_LABEL is a character column -- listed as DROP below.
;              ID=0 in IGNORE removes the CSV column-header row.

$DATA   busulfan_sim.csv IGNORE@ IGNORE=(ID.EQ.0) COMMA

$INPUT  ID TIME AMT RATE DV MDV EVID CMT AGE WT HT SEX FLAG DROP

$SUBROUTINE ADVAN3 TRANS4
$PK
; -- Estimated Coefficients -- ;
TH_CL   = THETA(1)
TH_V    = THETA(2)
MAT_MAG = THETA(3)
K_MAT   = THETA(4)
TH_REGI = THETA(5)
TH_DAY  = THETA(6)
TH_Q = THETA(13)
TH_V2 = THETA(14)

; -- Size Descriptors -- ;
BMI = WT/((HT*HT)/10000);
IF(SEX.EQ.0) THEN
  FFM = (1.11 + ((1-1.11)/(1+((AGE/7.1)**(-1.1))))) * ((9270 * WT)/(8780 + (244 * BMI)))
ELSE
  FFM = (0.88 + ((1-0.88)/(1+((AGE/13.4)**(-12.7))))) * ((9270 * WT)/(6680 + (216 * BMI)))
ENDIF

SWTCL = FFM
SWTV = FFM

FSIZECL = (SWTCL/55) ** THETA(11)
FSIZEV = SWTV/55

; -- Sex Effect -- ;
FSEX = 1
IF(SEX.EQ.0) FSEX = THETA(12) ; Female

; -- Maturation Function -- ;
FMAT = (MAT_MAG + (1-MAT_MAG) * (1-EXP(-AGE*K_MAT)))

; -- Clearance Decrease over Time Effect -- ;
DROP = THETA(9)
SHAPE = THETA(10)
FTIME = (1+ DROP * EXP(-SHAPE * TIME))

; -- Inter-occasion variability -- ;
BOVCL = ETA(3)
IF (TIME.GE.24.AND.TIME.LT.48) BOVCL = ETA(4)
IF (TIME.GE.48.AND.TIME.LT.72) BOVCL = ETA(5)
IF (TIME.GE.72.AND.TIME.LT.96) BOVCL = ETA(6)
IF (TIME.GE.96.AND.TIME.LT.9999) BOVCL = ETA(7)
BOVV = ETA(8)
IF (TIME.GE.24.AND.TIME.LT.48) BOVV = ETA(9)
IF (TIME.GE.48.AND.TIME.LT.72) BOVV = ETA(10)
IF (TIME.GE.72.AND.TIME.LT.96) BOVV = ETA(11)
IF (TIME.GE.96.AND.TIME.LT.9999) BOVV = ETA(12)

; -- Population Parameters -- ;
TVV = TH_V * FSIZEV * FSEX
TVCL = TH_CL * FMAT * FSIZECL * FTIME
TVQ = TH_Q * FSIZECL
TVV2 = TH_V2 * FSIZEV * FSEX

; -- Individual Parameters -- ;
CL = TVCL * EXP(ETA(1) + BOVCL)
V1  = TVV  * EXP(ETA(2) + BOVV)
Q = TVQ
V2  = TVV2  * EXP(ETA(13))

S1 = V1/1000 ; convert mg to mcg

$ERROR
CONC  = A(1)/S1
IPRED = CONC
PROP  = CONC * THETA(7)
ADD   = THETA(8)
SD    = SQRT(PROP*PROP + ADD*ADD)
Y     = CONC + SD*ERR(1)

TH_CL   = THETA(1)
TH_V    = THETA(2)
MAT_MAG = THETA(3)
HL      = THETA(4)
TH_REGI = THETA(5)
TH_DAY  = THETA(6)

$THETA  11.5555 ; 1. TH_CL
 50.1056 ; 2. TH_V
 0.446221 ; 3. MAT-MAG
 1.1199 ; 4. K_MAT
 0.100036 FIX ; 5. H_REGI
 -0.135 FIX ; 6. TH_DAY
 0.0690238 ; 7. PROP error
 28.5474 ; 8. ADD error
 0.201437 ; 9. DROP
 0.0414926 ; 10. SHAPE
 0.725589 ; 11. allo ffm exp
 1.11433 ; 12. Sex effect on V
 3.79717 ; 13. Q
 4.76 ; 14. V2
$OMEGA  BLOCK(2)
 0.0555324  ;  1. IIV CL
 0.0401747 0.0418246  ;  2. IIV V1
$OMEGA  BLOCK(1)
 0.0134219  ; 3-7. IOV CL
$OMEGA  BLOCK(1) SAME
$OMEGA  BLOCK(1) SAME
$OMEGA  BLOCK(1) SAME
$OMEGA  BLOCK(1) SAME
$OMEGA  BLOCK(1)
 0.0125156  ; 8-12. IOV V1
$OMEGA  BLOCK(1) SAME
$OMEGA  BLOCK(1) SAME
$OMEGA  BLOCK(1) SAME
$OMEGA  BLOCK(1) SAME
$OMEGA  0.209624  ;     13. IIV V2
$SIGMA  1  FIX

$ESTIMATION MAXEVAL=9999 METHOD=1 PRINT=5 SIG=4 NOABORT ; FOCEI
$COVARIANCE PRINT=E

$TABLE ID TIME AMT RATE DV MDV EVID CMT AGE WT HT SEX FLAG IPRED IRES IWRES CWRES ETA1 ETA2 ETA3
  ONEHEADER NOPRINT FILE=tab_busulfan
