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

## Steady State induced by repeated bolus doses.
e1 <- ev(amt = 1200, ii = 12, ss = 1, cmt = 1)
mod %>% ev(e1) %>% mrgsim(end = 500) %>% plot
time <- seq(from = 1, to = 10, by = 1)
xdata <- mod %>% ev(e1) %>% mrgsim(Req = "GUT, CENT, PROL, TRANSIT, CIRC",
                                   end = -1, add = time,
                                   rescort = 3) %>% as.data.frame
xdata
# ID time          GUT      CENT        PROL     TRANSIT        CIRC
# 1   1    0 1.200000e+03  46.92858 -0.08650993 -0.08766715 -0.08766290
# 2   1    1 1.624023e+02 897.86458 -0.08691954 -0.08763443 -0.08766247
# 3   1    2 2.197877e+01 791.46490 -0.08761369 -0.08762332 -0.08766135
# 4   1    3 2.974503e+00 610.56695 -0.08808584 -0.08763114 -0.08766024
# 5   1    4 4.025552e-01 460.96537 -0.08833222 -0.08764989 -0.08765960
# 6   1    5 5.447993e-02 346.69437 -0.08840099 -0.08767287 -0.08765965
# 7   1    6 7.373058e-03 260.57211 -0.08833521 -0.08769507 -0.08766043
# 8   1    7 9.978351e-04 195.81933 -0.08816816 -0.08771281 -0.08766182
# 9   1    8 1.350418e-04 147.15449 -0.08792498 -0.08772347 -0.08766361
# 10  1    9 1.827653e-05 110.58336 -0.08762457 -0.08772519 -0.08766555
# 11  1   10 2.473349e-06  83.10090 -0.08728113 -0.08771668 -0.08766732

## Steady State induced by repeated infusions
e1 <- ev(amt = 1200, rate = 150, ii = 12, ss = 1, cmt = 1)
mod %>% ev(e1) %>% mrgsim(end = 500) %>% plot
time <- seq(from = 1, to = 10, by = 1)
xdata <- mod %>% ev(e1) %>% mrgsim(Req = "GUT, CENT, PROL, TRANSIT, CIRC",
                                   end = -1, add = time,
                                   rescort = 3) %>% as.data.frame
xdata
# ID time        GUT     CENT        PROL     TRANSIT        CIRC
# 1   1    0  0.0251597 181.3172 -0.08776322 -0.08770075 -0.08766391
# 2   1    1 64.8532585 212.8358 -0.08754360 -0.08769911 -0.08766506
# 3   1    2 73.6267877 283.1219 -0.08740574 -0.08769178 -0.08766603
# 4   1    3 74.8141559 342.2470 -0.08735655 -0.08768178 -0.08766669
# 5   1    4 74.9748488 387.5317 -0.08737769 -0.08767172 -0.08766700
# 6   1    5 74.9965962 421.6776 -0.08745217 -0.08766351 -0.08766701
# 7   1    6 74.9995393 447.3531 -0.08756680 -0.08765858 -0.08766681
# 8   1    7 74.9999377 466.6498 -0.08771162 -0.08765792 -0.08766653
# 9   1    8 74.9999915 481.1511 -0.08787912 -0.08766221 -0.08766631
# 10  1    9 10.1501451 415.4866 -0.08802304 -0.08767157 -0.08766632
# 11  1   10  1.3736728 319.5250 -0.08804504 -0.08768332 -0.08766667
