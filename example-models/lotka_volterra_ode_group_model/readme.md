- [Build](#org49bd101)
  - [Edit/Add `cmdstan/make/local`](#orgf8a8d8c)
  - [Build in `cmdstan`](#org0d267d6)
- [Run](#org060aaa1)
- [Results](#org3a79b6d)


<a id="org49bd101"></a>

# Build


<a id="orgf8a8d8c"></a>

## Edit/Add `cmdstan/make/local`

```sh
TORSTEN_MPI = 1                                         # flag on torsten's MPI solvers
-include $(MATH)make/setup_torsten.mk
CXXFLAGS += $(CXXFLAGS_MPI) -isystem /usr/local/include # path to MPI library's headers
LDFLAGS += $(LDFLAGS_MPI)
CC=mpicxx                                               # name of mpi compilers
CXX=mpicxx
```


<a id="org0d267d6"></a>

## Build in `cmdstan`

```sh
make ../example-models/lotka_volterra_ode_group_model/lv_group
```


<a id="org060aaa1"></a>

# Run

```sh
mpiexec -n 2 lv_group sample data file=lv_group.data.R
```


<a id="org3a79b6d"></a>

# Results

Three binaries are built and used to run:

-   sequential run using Stan's `rk45` integrator.
-   sequential run using Torsten's `rk45` integrator.
-   MPI run using Torsten's `rk45` group integrator(with 2, 4, 8 processes, respectively).

The wall time of sequential runs and MPI runs(in seconds) with ODE group size 4

| run                          | wall time(s) |
|---------------------------- |------------ |
| Sequential stan              | 2135         |
| sequential pmx               | 2570         |
| MPI pmx (n<sub>proc</sub>=1) | 2107         |
| MPI pmx (n<sub>proc</sub>=2) | 1187         |
| MPI pmx (n<sub>proc</sub>=4) | 715          |
| MPI pmx (n<sub>proc</sub>=8) | 467          |
