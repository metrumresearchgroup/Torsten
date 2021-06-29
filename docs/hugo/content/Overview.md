+++
title = "Overview"
author = ["Yi Zhang"]
date = 2021-06-25T00:00:00-07:00
lastmod = 2021-06-28T19:35:09-07:00
draft = false
weight = 2001
+++

Torsten is a collection of Stan functions to facilitate analysis of
pharmacometric data using Stan. The current version
includes:

-   Specific linear compartment models:
    -   One compartment model with first order absorption.
    -   Two compartment model with elimination from and first order absorption into central compartment
-   General linear compartment model described by a system of first-order \underline{linear} Ordinary Differential Equations (ODEs).
-   General compartment model described by a system of first order ODEs.
-   Mix compartment model with PK forcing function described by a linear one or two compartment model.

The models and data format are based on
NONMEM\textregistered{}[^fn:1]/NMTRAN/PREDPP
conventions including:

-   Recursive calculation of model predictions
    -   This permits piecewise constant covariate values
-   Bolus or constant rate inputs into any compartment
-   Single dose and multiple dose events
-   Handles steady state dosing events
-   Implemented NMTRAN data items include: TIME, EVID, CMT, AMT, RATE, ADDL, II, SS

In general, all real variables may be passed as model parameters. A
few exceptions apply /to functions which use a numerical
integrator(i.e. the general and the mix compartment
models). The below listed cases present technical difficulties, which we expect to
overcome in Torsten's next release:

-   In the case of a multiple truncated infusion rate dosing regimen:
    -   The bioavailability (F) and the amount (AMT) must be fixed.

This library provides Stan language functions that calculate amounts
in each compartment, given an event schedule and an ODE system.

[^fn:1]: NONMEM\textregistered{} is licensed and distributed by ICON Development Solutions.
