check_sampling_csv_info_matches <- function (a, b) 
{
  if (a$model_name != b$model_name) {
    return("Supplied CSV files were not generated wtih the same model!")
  }
  if ((length(a$model_params) != length(b$model_params)) || 
      !(all(a$model_params == b$model_params) && all(a$sampler_diagnostics == 
                                                     b$sampler_diagnostics))) {
    return("Supplied CSV files have samples for different parameters!")
  }
  dont_match_list <- c("id", "inverse_metric", "step_size", 
                       "seed", "init")
  for (name in names(a)) {
    if (!(name %in% dont_match_list) && (is.null(b[[name]]) || 
                                         all(a[[name]] != b[[name]]))) {
      return("Supplied CSV files do not match in all sampling settings!")
    }
  }
  NULL
}

read_sample_csv <- function (output_files) 
{
  sampling_info <- NULL
  warmup_draws_array <- list()
  warmup_sampler_diagnostics_draws <- list()
  post_warmup_draws_array <- list()
  post_warmup_sampler_diagnostics_draws <- list()
  inverse_metric = list()
  step_size = list()
  for (output_file in output_files) {
    checkmate::assert_file_exists(output_file, access = "r", 
                                  extension = "csv")
    if (is.null(sampling_info)) {
      sampling_info <- cmdstanr:::read_sample_info_csv(output_file)
      inverse_metric[[sampling_info$id]] <- sampling_info$inverse_metric
      step_size[[sampling_info$id]] <- sampling_info$step_size
      id <- sampling_info$id
    }
    else {
      csv_file_info <- cmdstanr:::read_sample_info_csv(output_file)
      error <- check_sampling_csv_info_matches(sampling_info, 
                                               csv_file_info)
      if (!is.null(error)) {
        stop(error)
      }
      sampling_info$id <- c(sampling_info$id, csv_file_info$id)
      inverse_metric[[csv_file_info$id]] <- csv_file_info$inverse_metric
      step_size[[csv_file_info$id]] <- csv_file_info$step_size
      id <- csv_file_info$id
    }
    draws <- utils::read.csv(output_file, header = TRUE, 
                             comment.char = "#")
    if (sampling_info$save_warmup == 1) {
      warmup_draws_array[[id]] <- draws[1:sampling_info$num_warmup/sampling_info$thin, 
                                        sampling_info$model_params]
      warmup_sampler_diagnostics_draws[[id]] <- draws[1:sampling_info$num_warmup/sampling_info$thin, 
                                                      sampling_info$sampler_diagnostics]
      post_warmup_draws_array[[id]] <- draws[(sampling_info$num_warmup/sampling_info$thin + 
                                                1):sampling_info$num_iter, sampling_info$model_params]
      post_warmup_sampler_diagnostics_draws[[id]] <- draws[(sampling_info$num_warmup/sampling_info$thin + 
                                                              1):sampling_info$num_iter, sampling_info$sampler_diagnostics]
    }
    else {
      warmup_draws_array <- NULL
      warmup_sampler_diagnostics_draws <- NULL
      post_warmup_draws_array[[id]] <- draws[, sampling_info$model_params]
      post_warmup_sampler_diagnostics_draws[[id]] <- draws[, 
                                                           sampling_info$sampler_diagnostics]
    }
  }
  sampling_info$model_params <- cmdstanr:::repair_variable_names(sampling_info$model_params)
  num_chains <- length(sampling_info$id)
  if (!is.null(warmup_draws_array)) {
    if (!is.null(warmup_draws_array) && (length(warmup_draws_array) > 
                                         0)) {
      warmup_draws_array <- posterior::as_draws_array(array(unlist(do.call(rbind, 
                                                                           warmup_draws_array)), dim = c(sampling_info$num_warmup/sampling_info$thin, 
                                                                                                         num_chains, length(sampling_info$model_params)), 
                                                            dimnames = list(NULL, NULL, sampling_info$model_params)))
    }
    if (!is.null(warmup_sampler_diagnostics_draws) && (length(warmup_sampler_diagnostics_draws) > 
                                                       0)) {
      warmup_sampler_diagnostics_draws <- posterior::as_draws_array(array(unlist(do.call(rbind, 
                                                                                         warmup_sampler_diagnostics_draws)), dim = c(sampling_info$num_warmup/sampling_info$thin, 
                                                                                                                                     num_chains, length(sampling_info$sampler_diagnostics)), 
                                                                          dimnames = list(NULL, NULL, sampling_info$sampler_diagnostics)))
    }
  }
  if (!is.null(post_warmup_draws_array) && (length(post_warmup_draws_array) > 
                                            0)) {
    post_warmup_draws_array <- posterior::as_draws_array(array(unlist(do.call(rbind, 
                                                                              post_warmup_draws_array)), dim = c(sampling_info$num_samples/sampling_info$thin, 
                                                                                                                 num_chains, length(sampling_info$model_params)), 
                                                               dimnames = list(NULL, NULL, sampling_info$model_params)))
  }
  if (!is.null(post_warmup_sampler_diagnostics_draws) && (length(post_warmup_sampler_diagnostics_draws) > 
                                                          0)) {
    post_warmup_sampler_diagnostics_draws <- posterior::as_draws_array(array(unlist(do.call(rbind, 
                                                                                            post_warmup_sampler_diagnostics_draws)), dim = c(sampling_info$num_samples/sampling_info$thin, 
                                                                                                                                             num_chains, length(sampling_info$sampler_diagnostics)), 
                                                                             dimnames = list(NULL, NULL, sampling_info$sampler_diagnostics)))
  }
  sampling_info$inverse_metric <- NULL
  sampling_info$step_size <- NULL
  sampling_info$num_iter <- NULL
  list(sampling_info = sampling_info, inverse_metric = inverse_metric, 
       step_size = step_size, warmup_draws = warmup_draws_array, 
       post_warmup_draws = post_warmup_draws_array, warmup_sampler_diagnostics = warmup_sampler_diagnostics_draws, 
       post_warmup_sampler_diagnostics = post_warmup_sampler_diagnostics_draws)
}

