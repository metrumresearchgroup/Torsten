- [Build](#orgf9267f0)
  - [Edit/Add `cmdstan/make/local`](#org5d1ec6b)
  - [Build in `cmdstan`](#org70610a8)
- [Run](#org779a276)
- [Results](#orgaffa7a3)

TTPN model from Hands-on session 4 of the advanced Stan by Bill Gillespie


<a id="orgf9267f0"></a>

# Build


<a id="org5d1ec6b"></a>

## Edit/Add `cmdstan/make/local`

```sh
TORSTEN_MPI = 1                                         # flag on torsten's MPI solvers
CXXFLAGS += -isystem /usr/local/include                 # path to MPI library's headers
```


<a id="org70610a8"></a>

## Build in `cmdstan`

```sh
make ../example-models/ttpn2/ttpn2_group
```


<a id="org779a276"></a>

# Run

```sh
mpiexec -n 2 ttpn2_group sample num_warmup=500 num_samples=500 data file=ttpn2.data.R init=ttpn2.init.R
```


<a id="orgaffa7a3"></a>

# Results

Two binaries are built and used to run:

-   sequential run using Torsten's `rk45` solver `pmx_solve_rk45` by looping through the population.
-   MPI run using Torsten's group solver `pmx_solve_group_rk45` (with 1, 2, 4, 8, 16, 32, 64 processes, respectively). Note that with 64 processes the number of processes is greater than the number of subjects(60).

The wall time of sequential runs and MPI runs(in seconds) with ODE group size 4

| run                           | wall time(s) | Speedup |
|----------------------------- |------------ |------- |
| sequential pmx(v0.86)         | 1251         | 1.000   |
| MPI pmx (n<sub>proc</sub>=1)  | 1419         | 0.882   |
| MPI pmx (n<sub>proc</sub>=2)  | 862          | 1.451   |
| MPI pmx (n<sub>proc</sub>=4)  | 501          | 2.497   |
| MPI pmx (n<sub>proc</sub>=8)  | 411          | 3.044   |
| MPI pmx (n<sub>proc</sub>=16) | 289          | 4.329   |
| MPI pmx (n<sub>proc</sub>=32) | 263          | 4.757   |
| MPI pmx (n<sub>proc</sub>=64) | 230          | 5.439   |
