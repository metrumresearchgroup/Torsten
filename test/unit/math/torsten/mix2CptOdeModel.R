## Simulate solutions for unit tests in Torsten
rm(list = ls())
gc()

# .libPaths("~/svn-StanPmetrics/script/lib")
library(mrgsolve)
library(dplyr)

## Simulate ME-2 plasma concentrations and ANC values
## for a simplified Friberg-Karlsson model
code <- '
$PARAM CL = 10, Q = 15, VC = 35, VP = 105, KA = 2.0, MTT = 125, 
Circ0 = 5, alpha = 3E-4, gamma = 0.17

$SET delta=0.1 // simulation grid

$CMT GUT CENT PERI PROL TRANSIT CIRC

$MAIN

// Reparametrization
double k10 = CL / VC;
double k12 = Q / VC;
double k21 = Q / VP;
double ktr = 4/MTT;

$ODE 
dxdt_GUT = -KA * GUT;
dxdt_CENT = KA * GUT - (k10 + k12) * CENT + k21 * PERI;
dxdt_PERI = k12 * CENT - k21 * PERI;
dxdt_PROL = ktr * (PROL + Circ0) * ((1 - alpha * CENT/VC) * pow(Circ0/(CIRC + Circ0),gamma) - 1);
dxdt_TRANSIT = ktr * (PROL - TRANSIT);
dxdt_CIRC = ktr * (TRANSIT - CIRC);
'

mod <- mread("acum", tempdir(), code)
e1 <- ev(amt = 10000)
mod %>% ev(e1) %>% mrgsim(end = 500) %>% plot # plot data

## save some data for unit tests (see amounts at t = 1 hour, with no noise)
time <- seq(from = 0.25, to = 2, by = 0.25)
time <- c(time, 4)
xdata <- mod %>% ev(e1) %>% mrgsim(Req = "GUT, CENT, PERI, PROL, TRANSIT, CIRC",
                                   end = -1, add = time,
                                   rescort = 3) %>% as.data.frame
xdata
# ID time          GUT     CENT      PERI          PROL       TRANSIT          CIRC
# 1   1 0.00 10000.000000    0.000    0.0000  0.0000000000  0.000000e+00  0.000000e+00
# 2   1 0.25  6065.306597 3579.304  212.1623 -0.0006874417 -1.933282e-06 -3.990995e-09
# 3   1 0.50  3678.794412 5177.749  678.8210 -0.0022297559 -1.318812e-05 -5.608541e-08
# 4   1 0.75  2231.301599 5678.265 1233.3871 -0.0041121287 -3.824790e-05 -2.508047e-07
# 5   1 1.00  1353.352829 5597.489 1787.0134 -0.0060546255 -7.847821e-05 -7.039447e-07
# 6   1 1.25   820.849983 5233.332 2295.7780 -0.0079139199 -1.335979e-04 -1.533960e-06
# 7   1 1.50   497.870681 4753.865 2741.1870 -0.0096246889 -2.025267e-04 -2.852515e-06
# 8   1 1.75   301.973832 4250.712 3118.6808 -0.0111649478 -2.838628e-04 -4.760265e-06
# 9   1 2.00   183.156387 3771.009 3430.9355 -0.0125357580 -3.761421e-04 -7.345522e-06
# 10  1 4.00     3.354626 1601.493 4374.6747 -0.0192607813 -1.370742e-03 -5.951920e-05





