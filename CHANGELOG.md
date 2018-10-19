# 0.85 <span class="timestamp-wrapper"><span class="timestamp">&lt;2018-10-20 Sat&gt;</span></span>


## Added

-   Dosing rate as parameter


## Changed

-   Update with Stan version 2.18.0.


# 0.84 <span class="timestamp-wrapper"><span class="timestamp">&lt;2018-02-24 Sat&gt;</span></span>


## Added

-   Piecewise linear interpolation function.
-   Univariate integral functions.


## Changed

-   Update with Stan version 2.17.1.
-   Minor revisions to User Manual.
-   Bugfixes.


# 0.83 <span class="timestamp-wrapper"><span class="timestamp">&lt;2017-08-02 Wed&gt;</span></span>


## Added

-   Work with TorstenHeaders
-   Each chain has a different initial estimate


## Changed

-   User manual
-   Fix misspecification in ODE system for TwoCpt example.
-   Other bugfixes


# 0.82 <span class="timestamp-wrapper"><span class="timestamp">&lt;2017-01-29 Sun&gt;</span></span>


## Added

-   Allow parameter arguments to be passed as 1D or 2D arrays
-   More unit tests
-   Unit tests check automatic differentiation against finite differentiation.


## Changed

-   Split the parameter argument into three arguments: pMatrix (parameters for the ODEs &#x2013; note: for linOdeModel, pMatrix is replaced by the constant rate matrix K), biovar (parameters for the biovariability), and tlag (parameters for the lag time).
-   bugfixes


# 0.81 <span class="timestamp-wrapper"><span class="timestamp">&lt;2016-09-27 Tue&gt;</span></span>


## Added

linCptModel (linear compartmental model) function


# 0.80a <span class="timestamp-wrapper"><span class="timestamp">&lt;2016-09-21 Wed&gt;</span></span>


## Added

check<sub>finite</sub> statements in pred<sub>1</sub> and pred<sub>2</sub> to reject metropolis proposal if initial conditions are not finite
