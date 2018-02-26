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

See the user manual (`torstenManual.pdf`) for more
information and guidance on the examples. If you have any
questions, please raise an issue on GitHub or send us an
e-mail at billg@metrumrg.com. 

