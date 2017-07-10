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
$PARAM CL = 10, V = 80, KA = 1.2
$CMT GUT CENT
$GLOBAL
#define CP (CENT/V)

$PKMODEL ncmt = 1, depot = TRUE
$SIGMA 0.01  // variance
$TABLE
capture DV = CP * exp(EPS(1));
$CAPTURE CP
'

## Input 1: multiple truncated infusions
mod <- mread("acum", tempdir(), code)
e1 <- ev(amt = 1200, rate = 1200, addl = 14, ii = 12)
mod %>% ev(e1) %>% mrgsim(end = 500) %>% plot # plot data

## save some data for unit tests (see amounts at t = 1 hour, with no noise)
time <- seq(from = 0.25, to = 2, by = 0.25)
time <- c(time, 4)
xdata <- mod %>% ev(e1) %>% mrgsim(Req = "GUT, CENT",
                                   end = -1, add = time,
                                   rescort = 3) %>% as.data.frame
xdata
# ID time       GUT      CENT
# 1   1 0.00   0.00000   0.00000
# 2   1 0.25 259.18178  40.38605
# 3   1 0.50 451.18836 145.61440
# 4   1 0.75 593.43034 296.56207
# 5   1 1.00 698.80579 479.13371
# 6   1 1.25 517.68806 642.57025
# 7   1 1.50 383.51275 754.79790
# 8   1 1.75 284.11323 829.36134
# 9   1 2.00 210.47626 876.28631
# 10  1 4.00  19.09398 844.11769

## Input 2: multiple truncated infusions. Steady state.
e1 <- ev(amt = 1200, rate = 150, addl = 14, ii = 6, ss = 1)
mod %>% ev(e1) %>% mrgsim(end = 100) %>% plot # plot data

## save some data for unit tests (see amounts at t = 1 hour, with no noise)
time <- seq(from = 0.25, to = 2, by = 0.25)
time <- c(time, 4)
xdata <- mod %>% ev(e1) %>% mrgsim(Req = "GUT, CENT",
                                   end = -1, add = time,
                                   rescort = 3) %>% as.data.frame

xdata
# ID time       GUT       CENT
# 1   1 0.00   0.00000   0.000000
# 2   1 0.25  32.39772   5.048256
# 3   1 0.50  56.39855  18.201800
# 4   1 0.75  74.17879  37.070259
# 5   1 1.00  87.35072  59.891714
# 6   1 1.25  97.10873  85.369537
# 7   1 1.50 104.33764 112.551538
# 8   1 1.75 109.69295 140.740426
# 9   1 2.00 113.66026 169.427503
# 10  1 4.00 123.97128 388.679360


## Input 3: multiple drug infusion, where the duration of the infusion
## is longer than the interdose interval.
e1 <- ev(amt = 1200, rate = 75, ii = 12, ss = 1, addl = 14)
mod %>% ev(e1) %>% mrgsim(end = 100) %>% plot # plot data

## save some data for unit tests (see amounts at t = 1 hour, with no noise)
time <- seq(from = 0.25, to = 2, by = 0.25)
time <- c(time, 4)
xdata <- mod %>% ev(e1) %>% mrgsim(Req = "GUT, CENT",
                                   end = -1, add = time,
                                   rescort = 3) %>% as.data.frame
xdata
# ID time       GUT     CENT
# 1   1 0.00  62.50420 724.7889
# 2   1 0.25  78.70197 723.4747
# 3   1 0.50  90.70158 726.3310
# 4   1 0.75  99.59110 732.1591
# 5   1 1.00 106.17663 740.0744
# 6   1 1.25 111.05530 749.4253
# 7   1 1.50 114.66951 759.7325
# 8   1 1.75 117.34699 770.6441
# 9   1 2.00 119.33051 781.9027
# 10  1 4.00 124.48568 870.0308


## Input 4: constant rate infusion, steady state
## Doesn't look like I can use ss = 1 and ii <= 0. I'll hack it with ii = amt / rate.
e1 <- ev(amt = 1200, rate = 150, ss = 1, ii = 8)
mod %>% ev(e1) %>% mrgsim(end = 100) %>% plot

## save some data for unit tests (see amounts at t = 1 hour, with no noise)
time <- seq(from = 0, to = 22.5, by = 2.5)
xdata <- mod %>% ev(e1) %>% mrgsim(Req = "GUT, CENT",
                                   end = -1, add = time,
                                   rescort = 3) %>% as.data.frame

xdata
# ID time          GUT      CENT
# 1   1  0.0 0.000000e+00    0.0000
# 2   1  0.0 1.250000e+02 1200.0000
# 3   1  2.5 1.250000e+02 1200.0000
# 4   1  5.0 1.250000e+02 1200.0000
# 5   1  7.5 1.250000e+02 1200.0000
# 6   1 10.0 1.133974e+01 1030.5725
# 7   1 12.5 5.645726e-01  762.6137
# 8   1 15.0 2.810842e-02  558.3698
# 9   1 17.5 1.399436e-03  408.5335
# 10  1 20.0 6.967380e-05  298.8906
# 11  1 22.5 3.468854e-06  218.6731
