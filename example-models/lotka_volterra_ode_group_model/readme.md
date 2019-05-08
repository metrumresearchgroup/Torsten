- [Build](#org034b2b5)
  - [Edit/Add `cmdstan/make/local`](#org427b3cf)
  - [Build in `cmdstan`](#org33b2ed2)
- [Run](#org7eba283)
- [Results](#org2fe173c)


<a id="org034b2b5"></a>

# Build


<a id="org427b3cf"></a>

## Edit/Add `cmdstan/make/local`

```sh
TORSTEN_MPI = 1                                         # flag on torsten's MPI solvers
-include $(MATH)make/setup_torsten.mk
CXXFLAGS += $(CXXFLAGS_MPI) -isystem /usr/local/include # path to MPI library's headers
LDFLAGS += $(LDFLAGS_MPI)
CC=mpicxx                                               # name of mpi compilers
CXX=mpicxx
```


<a id="org33b2ed2"></a>

## Build in `cmdstan`

```sh
make ../example-models/lotka_volterra_ode_group_model/lv_group
```


<a id="org7eba283"></a>

# Run

```sh
mpiexec -n 2 lv_group sample data file=lv_group.data.R
```


<a id="org2fe173c"></a>

# Results

Three binaries are built and used to run:

-   sequential run using Stan's `rk45` integrator.
-   sequential run using Torsten's `rk45` integrator.
-   MPI run using Torsten's `rk45` group integrator(with 2, 4, 8 processes, respectively).

The wall time of sequential runs and MPI runs(in seconds) with ODE group size 4

| run                          | wall time(s) |
|---------------------------- |------------ |
| Sequential stan              | 2600         |
| sequential pmx               | 2679         |
| MPI pmx (n<sub>proc</sub>=2) | 1200         |
| MPI pmx (n<sub>proc</sub>=4) | 696          |
| MPI pmx (n<sub>proc</sub>=8) | 467          |
