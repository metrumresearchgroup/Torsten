## Simulate solutions for unit tests in Torsten
rm(list = ls())
gc()

# .libPaths("~/svn-StanPmetrics/script/lib")
library(mrgsolve)
library(dplyr)

## Simulate ME-2 plasma concentrations and ANC values
## for a simplified Friberg-Karlsson model
code <- '
$PARAM CL = 10, VC = 35, KA = 2.0, MTT = 125, 
Circ0 = 5, alpha = 3E-4, gamma = 0.17

$SET delta=0.1 // simulation grid

$CMT GUT CENT PROL TRANSIT CIRC

$MAIN

// Reparametrization
double k10 = CL / VC;
double ktr = 4/MTT;

$ODE 
dxdt_GUT = -KA * GUT;
dxdt_CENT = KA * GUT - k10 * CENT;
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
xdata <- mod %>% ev(e1) %>% mrgsim(Req = "GUT, CENT, PROL, TRANSIT, CIRC",
                                   end = -1, add = time,
                                   rescort = 3) %>% as.data.frame
xdata
# ID time          GUT     CENT          PROL       TRANSIT          CIRC
# 1   1 0.00 10000.000000    0.000  0.0000000000  0.000000e+00  0.000000e+00
# 2   1 0.25  6065.306623 3786.208 -0.0007126788 -1.985660e-06 -4.076891e-09
# 3   1 0.50  3678.794399 5821.649 -0.0023972982 -1.390885e-05 -5.850529e-08
# 4   1 0.75  2231.301565 6813.189 -0.0045843443 -4.140146e-05 -2.670477e-07
# 5   1 1.00  1353.352804 7188.323 -0.0069950553 -8.713363e-05 -7.646749e-07
# 6   1 1.25   820.849975 7205.188 -0.0094660436 -1.520317e-04 -1.698982e-06
# 7   1 1.50   497.870679 7019.273 -0.0119035120 -2.360139e-04 -3.219375e-06
# 8   1 1.75   301.973833 6723.888 -0.0142554906 -3.384324e-04 -5.470933e-06
# 9   1 2.00   183.156389 6374.696 -0.0164950219 -4.583390e-04 -8.591089e-06
# 10  1 4.00     3.354626 3716.663 -0.0300528403 -1.915294e-03 -7.813746e-05


e1 <- ev(amt = 12000, rate = 12000, addl = 14, ii = 12)
mod %>% ev(e1) %>% mrgsim(end = 500) %>% plot # plot data
mod %>% ev(e1) %>% mrgsim(end = 500) %>% plot # plot data

## save some data for unit tests (see amounts at t = 1 hour, with no noise)
time <- seq(from = 0.25, to = 2, by = 0.25)
time <- c(time, 4)
xdata <- mod %>% ev(e1) %>% mrgsim(Req = "GUT, CENT, PROL, TRANSIT, CIRC",
                                   end = -1, add = time,
                                   rescort = 3) %>% as.data.frame
xdata
# ID time        GUT      CENT          PROL       TRANSIT          CIRC
# 1   1 0.00    0.00000    0.0000  0.000000e+00  0.000000e+00  0.000000e+00
# 2   1 0.25 2360.81604  623.6384 -7.461806e-05 -1.531363e-07 -2.492751e-10
# 3   1 0.50 3792.72335 2098.1390 -5.238333e-04 -2.201387e-06 -7.280973e-09
# 4   1 0.75 4661.21904 4013.1415 -1.562825e-03 -1.006606e-05 -5.066933e-08
# 5   1 1.00 5187.98830 6124.9596 -3.296762e-03 -2.887519e-05 -1.964005e-07
# 6   1 1.25 3146.67405 7667.0022 -5.691112e-03 -6.411818e-05 -5.529458e-07
# 7   1 1.50 1908.55429 8329.8566 -8.448480e-03 -1.198063e-04 -1.267227e-06
# 8   1 1.75 1157.59668 8478.2378 -1.133512e-02 -1.976540e-04 -2.507387e-06
# 9   1 2.00  702.11787 8332.0619 -1.421568e-02 -2.979247e-04 -4.447571e-06
# 10  1 4.00   12.85974 5152.8451 -3.269696e-02 -1.786827e-03 -6.367482e-05
