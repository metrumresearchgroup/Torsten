// This stan code file declares NLTRAN-like event data variables.
// to use this file, include it at the beginning of the "data" block, i.e.
//
// data {
// #include nltran_data.stan
// ...
// }
//
// and then build the model with include path specified:
//
// cmdstanr::cmdstan_model("model_file", include_paths="path_to_torsten/example-models/include")

              int<lower  = 1>  nId;		// number of subjects
              int<lower  = 1>  nt;		// number of events (rows in data set)
              int<lower  = 1>  nObs;		// number of PK observations
  array[nObs] int<lower  = 1>  iObs;		// event indices for PK observations
  array[nt]   real<lower = 0>  amt;		// AMT
  array[nt]   real<lower = 0>  time;		// TIME
  array[nt]   int<lower  = 1>  cmt;		// CMT
  array[nt]   int<lower  = 0>  evid;		// EVID
