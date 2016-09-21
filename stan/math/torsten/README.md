<b> Torsten </b> is a library of C++ functions that support applications of Stan in Pharmacometrics. The current prototype provides:
* One and Two Compartment Model functions, that compute solutions analytically
* General Compartment Model functions, that compute solutions numerically
* A flexible mechanism that handles the event schedule of clinical trials

This prototype is still under development and has been uploaded to facilitate working with the community of Stan developers. The current prototype was written by Charles Margossian, Bill Gillespie, and Metrum Research Group, LLC.

Version: 0.81

Licensing
---------
The Torsten library is licensed under the BSD 3-clause license. 

Update
------
09/21
* Add check_finite statements in pred_1 and pred_2 to reject metropolis proposal if initial conditions is not finite
* Create a C++ unit test for PKModelOneCpt() for quick tests
