- [Benchmark](#orgc770c59)



<a id="orgc770c59"></a>

# Benchmark

Dummy model `vec.stan` and `non_vec.stan` contains a minimalist implementation using `wiener_lpdf` for likelihood calculation. They are used to examine the effect of vectorization of `wiener_lpdf` function on performance.

First download the repo:

```bash
git clone -b wiener_demo git@github.com:metrumresearchgroup/Torsten.git
```

To build the models, do

```bash
cd Torsten/cmdstan
make examples/wiener_demo/vec
make examples/wiener_demo/non_vec
```

To run, do

```bash
./non_vec sample num_warmup=10 num_samples=0 data file=data.json init=init.json random seed=987 output file=output.non_vec.csv
./vec sample num_warmup=10 num_samples=0 data file=data.json init=init.json random seed=987 output file=output.vec.csv
```

The total time can be compared as

```bash
grep -A 2 "Elapsed" *.csv
output.non_vec.csv:#  Elapsed Time: 80.249 seconds (Warm-up)
output.non_vec.csv-#                0 seconds (Sampling)
output.non_vec.csv-#                80.249 seconds (Total)
--
output.vec.csv:#   Elapsed Time: 82.575 seconds (Warm-up)
output.vec.csv-#                 0 seconds (Sampling)
output.vec.csv-#                 82.575 seconds (Total)
```
