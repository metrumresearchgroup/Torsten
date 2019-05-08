- [Description](#org48a23c6)
- [Build](#orgf1aba9c)
  - [Edit/Add `cmdstan/make/local`](#org8081c39)
  - [Build in `cmdstan`](#org51e3b15)
- [Run](#org5967c70)
- [Results](#orgea6b51c)


<a id="org48a23c6"></a>

# Description

Parameter inference for the harmonic oscillator ODE using Torsten's ODE group integrator `pmx_integrate_ode_group_rk45`.


<a id="orgf1aba9c"></a>

# Build


<a id="org8081c39"></a>

## Edit/Add `cmdstan/make/local`

```sh
TORSTEN_MPI = 1                                         # flag on torsten's MPI solvers
-include $(MATH)make/setup_torsten.mk
CXXFLAGS += $(CXXFLAGS_MPI) -isystem /usr/local/include # path to MPI library's headers
LDFLAGS += $(LDFLAGS_MPI)
CC=mpicxx                                               # name of mpi compilers
CXX=mpicxx
```


<a id="org51e3b15"></a>

## Build in `cmdstan`

```sh
make ../example-models/harmonic_oscillator_ode_group_model/sho_group
```


<a id="org5967c70"></a>

# Run

```sh
mpiexec -n 2 sho_group sample data file=sho_group.data.R
```


<a id="orgea6b51c"></a>

# Results

Three chains are run using

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
| Sequential stan               | 262          |
| sequential pmx                | 276          |
| MPI pmx (n<sub>proc</sub>=1)  | 253          |
| MPI pmx (n<sub>proc</sub>=2)  | 119          |
| MPI pmx (n<sub>proc</sub>=4)  | 62           |
| MPI pmx (n<sub>proc</sub>=8)  | 74           |
| MPI pmx (n<sub>proc</sub>=16) | 58           |
