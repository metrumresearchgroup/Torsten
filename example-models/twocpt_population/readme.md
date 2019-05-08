- [Description](#org34faa9f)
- [Build](#org95832f5)
  - [Edit/Add `cmdstan/make/local`](#orgf50e734)
  - [Build in `cmdstan`](#orgb1683f8)
- [Run](#orgcf5e02b)
- [Results](#org81b1f68)


<a id="org34faa9f"></a>

# Description

Parameter inference for the two-compartment solved using numerical integrator `pmx_integrate_ode_group_bdf`.


<a id="org95832f5"></a>

# Build


<a id="orgf50e734"></a>

## Edit/Add `cmdstan/make/local`

```sh
TORSTEN_MPI = 1                                         # flag on torsten's MPI solvers
-include $(MATH)make/setup_torsten.mk
CXXFLAGS += $(CXXFLAGS_MPI) -isystem /usr/local/include # path to MPI library's headers
LDFLAGS += $(LDFLAGS_MPI)
CC=mpicxx                                               # name of mpi compilers
CXX=mpicxx
```


<a id="orgb1683f8"></a>

## Build in `cmdstan`

```sh
make ../example-models/twocpt_population/twocpt_population
```


<a id="orgcf5e02b"></a>

# Run

```sh
mpiexec -n 2 twocpt_population sample data file=twocpt_population.data.R init=twocpt_population.init.R
```


<a id="org81b1f68"></a>

# Results

MPI run results with 1, 2, 4, 8 processes, respectively, for a population size 8.

| run                | wall time(s) |
|------------------ |------------ |
| n<sub>proc</sub>=1 | 10791        |
| n<sub>proc</sub>=2 | 5282         |
| n<sub>proc</sub>=4 | 2847         |
| n<sub>proc</sub>=8 | 2971         |
