- [Build](#org42d552e)
  - [Edit/Add `cmdstan/make/local`](#org2dfef2b)
  - [Build in `cmdstan`](#org2c18206)
- [Run](#org6c3b73f)
- [Results](#org537e1cc)

Predator-Prey model based on <https://mc-stan.org/users/documentation/case-studies/lotka-volterra-predator-prey.html>


<a id="org42d552e"></a>

# Build


<a id="org2dfef2b"></a>

## Edit/Add `cmdstan/make/local`

```sh
TORSTEN_MPI = 1                                         # flag on torsten's MPI solvers
CXXFLAGS += -isystem /usr/local/include                 # path to MPI library's headers
```


<a id="org2c18206"></a>

## Build in `cmdstan`

```sh
make ../example-models/lotka_volterra_ode_group_model/lv_group
```


<a id="org6c3b73f"></a>

# Run

```sh
mpiexec -n 2 lv_group sample data file=lv_group.data.R
```


<a id="org537e1cc"></a>

# Results

Three binaries are built and used to run:

-   sequential run using Stan's `rk45` integrator.
-   sequential run using Torsten's `rk45` integrator.
-   MPI run using Torsten's `rk45` group integrator(with 2, 4, 8 processes, respectively).

The wall time of sequential runs and MPI runs(in seconds) with ODE group size 16.

| run                          | wall time(s) |
|---------------------------- |------------ |
| Sequential stan              | 2135         |
| sequential pmx               | 2570         |
| MPI pmx (n<sub>proc</sub>=1) | 2107         |
| MPI pmx (n<sub>proc</sub>=2) | 1187         |
| MPI pmx (n<sub>proc</sub>=4) | 715          |
| MPI pmx (n<sub>proc</sub>=8) | 467          |
