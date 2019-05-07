Torsten 0.84
------------

<b> Torsten </b> is a library of C++ functions that support applications of Stan in Pharmacometrics. The current version provides:
* Analytical solution for the one and two linear compartment model with a first-order absorption
* Matrix exponential solution for a general linear compartment model
* Numerical solution for general compartment model using built-in integrators:
  * Runge-Kutta 4th/5th (rk45) order for non-stiff ODE systems
  * Backward Differentiation (bdf) for stiff ODE systems
* Mixed solver for general compartment model with a forcing PK function:
  * Computes analytical solutions for forcing PK One and Two compartment model with a first-order absorption
  * Computes numerical solution for the forced compartment model using the rk45 or bdf method
  
__This prototype is still under development__ and has been uploaded to
facilitate working with the community of Stan developers. The current
version was written by Charles Margossian (@charlesm93), Yi Zhang
(@yizhang-cae), Bill Gillespie (@billgillespie), and Metrum Research
Group, LLC (@metrumresearchgroup. We have recieved extensive help from the Stan development team.

See the user manual (`torstenManual.pdf`) for more information and guidance on the examples. If you have any questions, please raise an issue on GitHub or send us an e-mail at billg@metrumrg.com. 

Licensing
---------
The Torsten library is open-source and licensed under the BSD 3-clause license. 


Install
-------
To install CmdStan with Torsten, run the shell script script `setupTorsten.sh`.
To install RStan with Torsten, install the R package TorstenHeaders
(https://github.com/metrumresearchgroup/TorstenHeaders) and run
install_torsten():
```r
devtools::install_github('metrumresearchgroup/TorstenHeaders')
library(torstenHeaders)
install_torsten()
```

We are working with Stan's development team to create a system to add and share Stan packages. In the mean time, users can download a forked version of Stan with Torsten from GitHub. The latest version of Torsten (v0.84) is compatible with Stan v2.17.1. Torsten is agnostic to which Stan interface you use.


Examples
---------
For each model, we provide the following files:
* *modelName*.stan
* *modelName*.data.R
* *modelName*.init.R
* *modelName*Simulation.R 

The simulation file can be used to create the data and initial estimate files. 

There are three Stan files for the  two compartment model: `TwoCptModel`, `LinTwoCptModel`, and `GenTwoCptModel`. We also provide more sophisticated PK/PD examples. `TwoCptModelPopulation` extends the two compartment model to the case where we have data from multiple patients with inter-individual variability. 

See the manual for more assistance.

Under the R directory, we provide scripts to run the Stan models.

C++ Code
--------
The C++ code for Torsten can be found on the following repos:

Math Library: https://github.com/metrumresearchgroup/math

Stan Grammar: https://github.com/metrumresearchgroup/stan
