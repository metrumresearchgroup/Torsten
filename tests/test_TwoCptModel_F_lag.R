context("TwoCptModel_F_lag")

library(mrgsolve)  ## tested with version 0.7.6.9029
library(testthat)
testmodel = "TwoCptModel_F_lag"

torstenDir <- Sys.getenv("TORSTEN_PATH")
try(if(torstenDir == "") stop("Missing environment variable 'TORSTEN_PATH'"))
source(file.path(torstenDir, "tests", "simulation_test_template.R"))

run_test_TwoCptModel_F_lag <- function(event, steptime) {
modelName <- "test_TwoCptModel_F_lag"

code <- '
$PARAM CL = 5, Q = 8, V2 = 20, V3 = 70, KA = 1.2

$PKMODEL cmt="GUT CENT PERI", depot = TRUE

$MAIN
F_CENT = 0.7;
ALAG_CENT = 0.25;

$GLOBAL
#define CP (CENT/V2)

$TABLE
capture DV = CP * exp(EPS(1));

$CAPTURE CP
'
par <- default_stan_parameters()
return (simulation_test(event, steptime, modelName, code, par))
}

dose = "IV doses with F and lag"
context(paste(testmodel, "with", dose))
test_that(dose, {
    e1 <- ev(amt = 80 * 1000, cmt = "CENT")
    time <- c(0, 1.0, 2.0) ## Observation times
    res <- run_test_TwoCptModel_F_lag(e1, time)
    fit <- res$fit
    data <- res$data
    do_tests(fit, data)
})
