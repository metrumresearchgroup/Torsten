<b> Torsten </b> is a library of C++ functions that support applications of Stan in Pharmacometrics. The current prototype provides:
* One and Two Compartment Model functions, that compute solutions analytically
* A Linear Compartment Model function, that computes solutions semi-analytically using matrix exponentials
* General Compartment Model functions, that compute solutions numerically using ODE integrators
* A flexible mechanism that handles the event schedule of clinical trials

This prototype is still under development and has been uploaded to facilitate working with the community of Stan developers. The current prototype was written by Charles Margossian, Bill Gillespie, and Metrum Research Group, LLC.

Version: 0.81

Licensing
---------
The Torsten library is licensed under the BSD 3-clause license. 

Updates
------
09/27 (0.81)
* Add linCptModel (linear compartmental model) function

09/21 (0.80a)
* Add check_finite statements in pred_1 and pred_2 to reject metropolis proposal if initial conditions are not finite
