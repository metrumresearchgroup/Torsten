Since metworx doesn't give users access to SGE PE modification, there's no viable way to use SGE for core-binding run with AWS. The following method for MPI runs uses `mpiexec` option to directly control that.

For model `foo`, run

```bash
mpiexec -l -bind-to core -n 4 -f hostfile foo sample save_warmup=1 adapt num_cross_chains=4 cross_chain_ess=400 data file=chem.data.R init=init.R random seed=123
```
