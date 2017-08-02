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

$CMT GUT CENT PERI PROL TRANSIT1 TRANSIT2 TRANSIT3 CIRC

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
dxdt_PROL = ktr * (PROL + Circ0) * 
  ((1 - alpha * CENT/VC) * pow(Circ0/(CIRC + Circ0),gamma) - 1);
dxdt_TRANSIT1 = ktr * (PROL - TRANSIT1);
dxdt_TRANSIT2 = ktr * (TRANSIT1 - TRANSIT2);
dxdt_TRANSIT3 = ktr * (TRANSIT2 - TRANSIT3);
dxdt_CIRC = ktr * (TRANSIT3 - CIRC);
'

mod <- mread("acum", tempdir(), code)
e1 <- ev(amt = 80 * 1000, ii = 12, ss = 1)
mod %>% ev(e1) %>% mrgsim(end = 200) %>% plot  # plot data

time <- seq(from = 0, to = 11.25, by = 1.25)
xdata <- mod %>% ev(e1) %>% mrgsim(Req = "GUT, CENT, PERI, PROL,
                                         TRANSIT1, TRANSIT2, TRANSIT3, CIRC",
                                   end = -1, add = time,
                                   rescort = 3) %>% as.data.frame
xdata

# ID  time          GUT     CENT     PERI      PROL  TRANSIT2  TRANSIT3      CIRC
# 1   1  0.00 0.000000e+00     0.00     0.00  0.000000  0.000000  0.000000  0.000000
# 2   1  0.00 8.000000e+04 11996.63 55694.35 -3.636308 -3.653933 -3.653748 -3.653622
# 3   1  1.25 6.566800e+03 53123.67 70649.28 -3.650990 -3.653910 -3.653755 -3.653627
# 4   1  2.50 5.390358e+02 34202.00 80161.15 -3.662446 -3.653883 -3.653761 -3.653632
# 5   1  3.75 4.424675e+01 23849.69 80884.40 -3.665321 -3.653870 -3.653765 -3.653637
# 6   1  5.00 3.631995e+00 19166.83 78031.24 -3.664114 -3.653876 -3.653769 -3.653642
# 7   1  6.25 2.981323e-01 16799.55 74020.00 -3.660988 -3.653896 -3.653774 -3.653647
# 8   1  7.50 2.447219e-02 15333.26 69764.65 -3.656791 -3.653926 -3.653779 -3.653653
# 9   1  8.75 2.008801e-03 14233.96 65591.05 -3.651854 -3.653957 -3.653786 -3.653658
# 10  1 10.00 1.648918e-04 13303.26 61607.92 -3.646317 -3.653983 -3.653793 -3.653663
# 11  1 11.25 1.353552e-05 12466.56 57845.10 -3.640244 -3.653995 -3.653801 -3.653668

xdata <- mod %>% ev(e1) %>% mrgsim(Req = "TRANSIT1",
                                   end = -1, add = time,
                                   rescort = 3) %>% as.data.frame

xdata
# ID  time  TRANSIT1
# 1   1  0.00  0.000000
# 2   1  0.00 -3.653620
# 3   1  1.25 -3.653172
# 4   1  2.50 -3.653349
# 5   1  3.75 -3.653782
# 6   1  5.00 -3.654219
# 7   1  6.25 -3.654550
# 8   1  7.50 -3.654722
# 9   1  8.75 -3.654708
# 10  1 10.00 -3.654488
# 11  1 11.25 -3.654050

