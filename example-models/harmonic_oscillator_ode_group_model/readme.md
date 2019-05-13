- [Description](#org577d566)
- [Build](#orgefda1ad)
  - [Edit/Add `cmdstan/make/local`](#org4de764d)
  - [Build in `cmdstan`](#org6d50c58)
- [Run](#org774e77f)
- [Results](#org510429e)


<a id="org577d566"></a>

# Description

Parameter inference for the harmonic oscillator ODE using Torsten's ODE group integrator `pmx_integrate_ode_group_rk45`.


<a id="orgefda1ad"></a>

# Build


<a id="org4de764d"></a>

## Edit/Add `cmdstan/make/local`

```sh
TORSTEN_MPI = 1                                         # flag on torsten's MPI solvers
-include $(MATH)make/setup_torsten.mk
CXXFLAGS += $(CXXFLAGS_MPI) -isystem /usr/local/include # path to MPI library's headers
LDFLAGS += $(LDFLAGS_MPI)
CC=mpicxx                                               # name of mpi compilers
CXX=mpicxx
```


<a id="org6d50c58"></a>

## Build in `cmdstan`

```sh
make ../example-models/harmonic_oscillator_ode_group_model/sho_group
```


<a id="org774e77f"></a>

# Run

```sh
mpiexec -n 2 sho_group sample data file=sho_group.data.R
```


<a id="org510429e"></a>

# Results

Three binaries are built and used to run:

-   sequential run using Stan's `bdf` integrator.
-   sequential run using Torsten's `bdf` integrator.
-   MPI run using Torsten's `bdf` group integrator(with 1, 2, 4, 8, 16 processes, respectively).

The wall time of sequential runs and MPI runs(in seconds) with ODE group size 4

| run                          | wall time(s) |
|---------------------------- |------------ |
| Sequential stan              | 39           |
| sequential pmx               | 33           |
| MPI pmx (n<sub>proc</sub>=1) | 34           |
| MPI pmx (n<sub>proc</sub>=2) | 18           |
| MPI pmx (n<sub>proc</sub>=4) | 10           |

The wall time of sequential runs and MPI runs(in seconds) with ODE group size 20

| run                           | wall time(s) |
|----------------------------- |------------ |
| Sequential stan               | 891          |
| sequential pmx                | 319          |
| MPI pmx (n<sub>proc</sub>=1)  | 476          |
| MPI pmx (n<sub>proc</sub>=2)  | 251          |
| MPI pmx (n<sub>proc</sub>=4)  | 136          |
| MPI pmx (n<sub>proc</sub>=8)  | 90           |
| MPI pmx (n<sub>proc</sub>=16) | 97           |
