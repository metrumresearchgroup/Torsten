context("MixedTwoCptModel_bdf")

library(mrgsolve)  ## tested with version 0.7.6.9029
library(testthat)

torstenDir <- Sys.getenv("TORSTEN_PATH")
try(if(torstenDir == "") stop("Missing environment variable 'TORSTEN_PATH'"))
source(file.path(torstenDir, "tests", "simulation_test_template.R"))

mrgcode <- '
$PARAM CL = 10, Q = 15, V2 = 35, V3 = 105, KA = 2.0, KOUT = 0.05, effect0=10, ec50=400

$CMT GUT CENT PERI RESP

$GLOBAL
#define CP (CENT/V2)

$TABLE
capture DV = RESP * exp(EPS(1));

$MAIN
double k10 = CL/V2;
double k12 = Q/V2;
double k21 = Q/V3;

$ODE
double Edrug = CP / (ec50 + CP);
double effect = RESP + effect0;
double kin0 = KOUT * effect0;
double kin = kin0 * (1.0 - Edrug);

dxdt_GUT = -KA*GUT;
dxdt_CENT = KA*GUT - (k10 + k12) * CENT + k21 * PERI;
dxdt_PERI = k12 * CENT - k21 * PERI;
dxdt_RESP = kin - KOUT * effect;

$CAPTURE CP
'

run_test_MixedTwoCptModel_bdf <- function(event, steptime) {
modelName <- "test_MixedTwoCptModel_bdf"
par <- default_stan_parameters()
return (simulation_test (event, steptime, modelName, mrgcode, par))
}
testmodel = "MixedTwoCptModel_bdf"
dose = "instantaneous bolus doses"
context(paste(testmodel, "with", dose))
test_that(dose, {
    e1 <- ev(amt = 80 * 1000)
    time <- c(0, 0.083, 0.167) ## Observation times
    res <- run_test_MixedTwoCptModel_bdf(e1, time)
    fit <- res$fit
    data <- res$data
    do_tests(fit, data)
})

dose = "infusion"
context(paste(testmodel, "with", dose))
test_that(dose, {
    e1 <- ev(amt = 80 * 1000, rate = 8000)
    time <- c(0, 0.083, 0.167) ## Observation times
    res <- run_test_MixedTwoCptModel_bdf(e1, time)
    fit <- res$fit
    data <- res$data
    do_tests(fit, data)
})

dose = "steady state with multiple bolus"
context(paste(testmodel, "with", dose))
test_that(dose, {
    e1 <- ev(amt = 80 * 1000, rate = 0, ss=1, ii=12, addl=2)
    time <- c(0, 12, 22) ## Observation times
    res <- run_test_MixedTwoCptModel_bdf(e1, time)
    fit <- res$fit
    data <- res$data
    do_tests(fit, data)
})

dose = "steady state with multiple infusion"
context(paste(testmodel, "with", dose))
test_that(dose, {
    e1 <- ev(amt = 10, rate=5, ss=1, ii=12, addl=2)
    time <- c(0, 12, 22) ## Observation times
    res <- run_test_MixedTwoCptModel_bdf(e1, time)
    fit <- res$fit
    data <- res$data
    do_tests(fit, data)
})
