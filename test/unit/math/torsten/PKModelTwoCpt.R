## FIX ME - Put simulations for all tests in PKModelOneCpt_test.cpp
## in this file.

## Simulate solutions for unit tests in Torsten
rm(list = ls())
gc()

# .libPaths("~/svn-StanPmetrics/script/lib")
library(mrgsolve)
library(dplyr)

## Simulate for One Compartment Model
code <- '
$PARAM CL = 5, Q = 8, V2 = 35, V3 = 105, KA = 1.2
$CMT GUT CENT PERI
$GLOBAL
#define CP (CENT/V2)

$PKMODEL ncmt = 2, depot = TRUE
$SIGMA 0.01  // variance
$TABLE
capture DV = CP * exp(EPS(1));
$CAPTURE CP
'

mod <- mread("acum", tempdir(), code)
e1 <- ev(amt = 1200, rate = 1200, addl = 14, ii = 12)
mod %>% ev(e1) %>% mrgsim(end = 500) %>% plot # plot data

## save some data for unit tests (see amounts at t = 1 hour, with no noise)
time <- seq(from = 0.25, to = 2, by = 0.25)
time <- c(time, 4)
xdata <- mod %>% ev(e1) %>% mrgsim(Req = "GUT, CENT, PERI",
                                   end = -1, add = time,
                                   rescort = 3) %>% as.data.frame
xdata
# ID time       GUT      CENT        PERI
# 1   1 0.00   0.00000   0.00000   0.0000000
# 2   1 0.25 259.18178  39.55748   0.7743944
# 3   1 0.50 451.18836 139.65573   5.6130073
# 4   1 0.75 593.43034 278.43884  17.2109885
# 5   1 1.00 698.80579 440.32663  37.1629388
# 6   1 1.25 517.68806 574.76950  65.5141658
# 7   1 1.50 383.51275 653.13596  99.2568509
# 8   1 1.75 284.11323 692.06145 135.6122367
# 9   1 2.00 210.47626 703.65965 172.6607082
# 10  1 4.00  19.09398 486.11014 406.6342765





