<b> Torsten </b> (math) is a library of C++ functions that supports applications of Stan in Pharmacometrics. The current prototype provides:
* One and Two Compartment (with first-order absoprtion) Model functions, that solve ordinary differential equations (ODEs) analytically
* A Linear ODE-based Model function, that solves ODEs using a matrix exponential solution
* General Compartment Model functions, that solve ODEs numerically using one of two ODE integrators available on stan (rk45 and bdf)
* A flexible mechanism that handles the event schedule of clinical trials and solves differential equations within that schedule.

This prototype is still under development and has been uploaded to facilitate working with the community of Stan developers. We welcome users to try Torsten and will gladly assist them. Please use gitHub to report issues and bugs, and request new features. People who wish to contribute to this project can email charlesm@metrumrg.com or billg@metrumrg.com. The current prototype was written by Charles Margossian, Bill Gillespie, and Metrum Research Group, LLC. This project is and will remain open-source.

Version: 0.82

Licensing
---------
The Torsten library is licensed under the BSD 3-clause license. 

Updates
------
01/29/2017 (0.82)
* Split the parameter argument into three arguments: pMatrix (parameters for the ODEs -- note: for linOdeModel, pMatrix is replaced by the constant rate matrix K), biovar (parameters for the biovariability), and tlag (parameters for the lag time).
* Allow parameter arguments to be passed as 1D or 2D arrays
* Increase the numbers of unit tests
* Unit tests check automatic differentiation against finite differentiation.
* Fix minor bugs and report issues

09/27/2016 (0.81)
* Add linCptModel (linear compartmental model) function

09/21/2016 (0.80a)
* Add check_finite statements in pred_1 and pred_2 to reject metropolis proposal if initial conditions are not finite
