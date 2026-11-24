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
## According to Bill, HACK DOESN'T WORK.
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

## Input 5: multiple truncated infusion, steady state
e1 <- ev(amt = 1200, rate = 150, ii = 16, ss = 1)
mod %>% ev(e1) %>% mrgsim(end = 100) %>% plot

time <- seq(from = 0, to = 22.5, by = 2.5)
xdata <- mod %>% ev(e1) %>% mrgsim(Req = "GUT, CENT",
                                   end = -1, add = time,
                                   rescort = 3) %>% as.data.frame

xdata
# ID time          GUT     CENT
# 1   1  0.0 0.000000e+00   0.0000
# 2   1  0.0 8.465519e-03 360.2470
# 3   1  2.5 1.187770e+02 490.4911
# 4   1  5.0 1.246902e+02 676.1759
# 5   1  7.5 1.249846e+02 816.5263
# 6   1 10.0 1.133898e+01 750.0054
# 7   1 12.5 5.645344e-01 557.3459
# 8   1 15.0 2.810651e-02 408.1926
# 9   1 17.5 1.399341e-03 298.6615
# 10  1 20.0 6.966908e-05 218.5065
# 11  1 22.5 3.468619e-06 159.8628

## Input 6: multiple doses IV (not absorption from the GUT)
## (value of ka should not matter)
e1 <- ev(amt = 1000, cmt = 2, ii = 12, addl = 14)
mod %>% ev(e1) %>% mrgsim(end = 250) %>% plot

time <- seq(from = 0, to = 2, by = 0.25)
time <- c(time, 4.0)

xdata <- mod %>% ev(e1) %>% mrgsim(Req = "GUT, CENT",
                                   end = -1, add = time,
                                   rescort = 3) %>% as.data.frame
xdata
# ID time GUT      CENT
# 1   1 0.00   0    0.0000
# 2   1 0.00   0 1000.0000
# 3   1 0.25   0  969.2332
# 4   1 0.50   0  939.4131
# 5   1 0.75   0  910.5104
# 6   1 1.00   0  882.4969
# 7   1 1.25   0  855.3453
# 8   1 1.50   0  829.0291
# 9   1 1.75   0  803.5226
# 10  1 2.00   0  778.8008
# 11  1 4.00   0  606.5307


## Case 7: Model with lag times
## Simulate for One Compartment Model with lag time
code <- '
$PARAM CL = 10, V = 80, KA = 1.2, ALG = 0
$CMT GUT CENT
$GLOBAL
#define CP (CENT/V)

$PKMODEL ncmt = 1, depot = TRUE

$MAIN
_ALAG(1) = ALG;

$SIGMA 0.01  // variance
$TABLE
capture DV = CP * exp(EPS(1));
$CAPTURE CP
'
mod <- mread("acum", tempdir(), code)
e1 <- ev(amt = 1200, addl = 14, ii = 12, ALG = 5)
mod %>% ev(e1) %>% mrgsim(end = 50) %>% plot # plot data

time <- seq(from = 0, to = 10, by = 1)
xdata <- mod %>% ev(e1) %>% mrgsim(Req = "GUT, CENT",
                                   end = -1, add = time,
                                   rescort = 3) %>% as.data.frame
xdata
# ID time        GUT     CENT
# 1   1    0   0.000000   0.0000
# 2   1    0   0.000000   0.0000
# 3   1    1   0.000000   0.0000
# 4   1    2   0.000000   0.0000
# 5   1    3   0.000000   0.0000
# 6   1    4   0.000000   0.0000
# 7   1    5   0.000000   0.0000
# 8   1    6 698.805788 479.1337
# 9   1    7 210.476259 876.2863
# 10  1    8  63.394231 909.8972
# 11  1    9  19.093975 844.1177
# 12  1   10   5.750995 757.3213
