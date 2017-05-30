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





