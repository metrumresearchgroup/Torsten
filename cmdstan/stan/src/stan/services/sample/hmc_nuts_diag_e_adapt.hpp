#ifndef STAN_SERVICES_SAMPLE_HMC_NUTS_DIAG_E_ADAPT_HPP
#define STAN_SERVICES_SAMPLE_HMC_NUTS_DIAG_E_ADAPT_HPP

#include <stan/math/prim.hpp>
#include <stan/callbacks/interrupt.hpp>
#include <stan/callbacks/logger.hpp>
#include <stan/callbacks/writer.hpp>
#include <stan/io/var_context.hpp>
#include <stan/mcmc/fixed_param_sampler.hpp>
#include <stan/services/error_codes.hpp>
#include <stan/mcmc/hmc/nuts/adapt_diag_e_nuts.hpp>
#include <stan/services/util/run_adaptive_sampler.hpp>
#include <stan/services/util/create_rng.hpp>
#include <stan/services/util/initialize.hpp>
#include <stan/services/util/inv_metric.hpp>
#include <vector>

namespace stan {
namespace services {
namespace sample {

/**
 * Runs HMC with NUTS with adaptation using diagonal Euclidean metric
 * with a pre-specified Euclidean metric.
 *
 * @tparam Model Model class
 * @param[in] model Input model to test (with data already instantiated)
 * @param[in] init var context for initialization
 * @param[in] init_inv_metric var context exposing an initial diagonal
              inverse Euclidean metric (must be positive definite)
 * @param[in] random_seed random seed for the random number generator
 * @param[in] chain chain id to advance the pseudo random number generator
 * @param[in] init_radius radius to initialize
 * @param[in] num_warmup Number of warmup samples
 * @param[in] num_samples Number of samples
 * @param[in] num_thin Number to thin the samples
 * @param[in] save_warmup Indicates whether to save the warmup iterations
 * @param[in] refresh Controls the output
 * @param[in] stepsize initial stepsize for discrete evolution
 * @param[in] stepsize_jitter uniform random jitter of stepsize
 * @param[in] max_depth Maximum tree depth
 * @param[in] delta adaptation target acceptance statistic
 * @param[in] gamma adaptation regularization scale
 * @param[in] kappa adaptation relaxation exponent
 * @param[in] t0 adaptation iteration offset
 * @param[in] init_buffer width of initial fast adaptation interval
 * @param[in] term_buffer width of final fast adaptation interval
 * @param[in] window initial width of slow adaptation interval
 * @param[in,out] interrupt Callback for interrupts
 * @param[in,out] logger Logger for messages
 * @param[in,out] init_writer Writer callback for unconstrained inits
 * @param[in,out] sample_writer Writer for draws
 * @param[in,out] diagnostic_writer Writer for diagnostic information
 * @return error_codes::OK if successful
 */
template <class Model>
int hmc_nuts_diag_e_adapt(
    Model& model, const stan::io::var_context& init,
    const stan::io::var_context& init_inv_metric, unsigned int random_seed,
    unsigned int chain, double init_radius,
    int num_cross_chains, int cross_chain_window, double cross_chain_rhat, int cross_chain_ess,
    int num_warmup, int num_samples,
    int num_thin, bool save_warmup, int refresh, double stepsize,
    double stepsize_jitter, int max_depth, double delta, double gamma,
    double kappa, double t0, unsigned int init_buffer, unsigned int term_buffer,
    unsigned int window, callbacks::interrupt& interrupt,
    callbacks::logger& logger, callbacks::writer& init_writer,
    callbacks::writer& sample_writer, callbacks::writer& diagnostic_writer) {
  boost::ecuyer1988 rng = util::create_rng(random_seed, chain);

  std::vector<int> disc_vector;
  std::vector<double> cont_vector = util::initialize(
      model, init, rng, init_radius, true, logger, init_writer);

  Eigen::VectorXd inv_metric;
  try {
    inv_metric = util::read_diag_inv_metric(init_inv_metric,
                                            model.num_params_r(), logger);
    util::validate_diag_inv_metric(inv_metric, logger);
  } catch (const std::domain_error& e) {
    return error_codes::CONFIG;
  }

  stan::mcmc::adapt_diag_e_nuts<Model, boost::ecuyer1988> sampler(model, rng);

  sampler.set_metric(inv_metric);
  sampler.set_nominal_stepsize(stepsize);
  sampler.set_stepsize_jitter(stepsize_jitter);
  sampler.set_max_depth(max_depth);

  sampler.get_stepsize_adaptation().set_mu(log(10 * stepsize));
  sampler.get_stepsize_adaptation().set_delta(delta);
  sampler.get_stepsize_adaptation().set_gamma(gamma);
  sampler.get_stepsize_adaptation().set_kappa(kappa);
  sampler.get_stepsize_adaptation().set_t0(t0);

  sampler.set_window_params(num_warmup, init_buffer, term_buffer, window,
                            logger);

  // cross chain adaptation
  sampler.set_cross_chain_adaptation_params(init_buffer, term_buffer, num_warmup,
                                            cross_chain_window, num_cross_chains,
                                            cross_chain_rhat, cross_chain_ess);
  mcmc::mpi_var_adaptation var_adapt(model.num_params_r(), num_warmup, cross_chain_window);
  sampler.set_cross_chain_metric_adaptation(&var_adapt);

  util::run_adaptive_sampler(
      sampler, model, cont_vector, num_warmup, num_samples, num_thin, refresh,
      save_warmup, rng, interrupt, logger, sample_writer, diagnostic_writer);

  return error_codes::OK;
}

/**
 * Runs HMC with NUTS with adaptation using diagonal Euclidean metric.
 *
 * @tparam Model Model class
 * @param[in] model Input model to test (with data already instantiated)
 * @param[in] init var context for initialization
 * @param[in] random_seed random seed for the random number generator
 * @param[in] chain chain id to advance the pseudo random number generator
 * @param[in] init_radius radius to initialize
 * @param[in] num_warmup Number of warmup samples
 * @param[in] num_samples Number of samples
 * @param[in] num_thin Number to thin the samples
 * @param[in] save_warmup Indicates whether to save the warmup iterations
 * @param[in] refresh Controls the output
 * @param[in] stepsize initial stepsize for discrete evolution
 * @param[in] stepsize_jitter uniform random jitter of stepsize
 * @param[in] max_depth Maximum tree depth
 * @param[in] delta adaptation target acceptance statistic
 * @param[in] gamma adaptation regularization scale
 * @param[in] kappa adaptation relaxation exponent
 * @param[in] t0 adaptation iteration offset
 * @param[in] init_buffer width of initial fast adaptation interval
 * @param[in] term_buffer width of final fast adaptation interval
 * @param[in] window initial width of slow adaptation interval
 * @param[in,out] interrupt Callback for interrupts
 * @param[in,out] logger Logger for messages
 * @param[in,out] init_writer Writer callback for unconstrained inits
 * @param[in,out] sample_writer Writer for draws
 * @param[in,out] diagnostic_writer Writer for diagnostic information
 * @return error_codes::OK if successful
 */
template <class Model>
int hmc_nuts_diag_e_adapt(
    Model& model, const stan::io::var_context& init, unsigned int random_seed,
    unsigned int chain, double init_radius,
    int num_cross_chains, int cross_chain_window, double cross_chain_rhat, int cross_chain_ess,
    int num_warmup, int num_samples,
    int num_thin, bool save_warmup, int refresh, double stepsize,
    double stepsize_jitter, int max_depth, double delta, double gamma,
    double kappa, double t0, unsigned int init_buffer, unsigned int term_buffer,
    unsigned int window, callbacks::interrupt& interrupt,
    callbacks::logger& logger, callbacks::writer& init_writer,
    callbacks::writer& sample_writer, callbacks::writer& diagnostic_writer) {
  stan::io::dump dmp
      = util::create_unit_e_diag_inv_metric(model.num_params_r());
  stan::io::var_context& unit_e_metric = dmp;

  return hmc_nuts_diag_e_adapt(
      model, init, unit_e_metric, random_seed, chain, init_radius,
      num_cross_chains, cross_chain_window, cross_chain_rhat, cross_chain_ess, num_warmup,
      num_samples, num_thin, save_warmup, refresh, stepsize, stepsize_jitter,
      max_depth, delta, gamma, kappa, t0, init_buffer, term_buffer, window,
      interrupt, logger, init_writer, sample_writer, diagnostic_writer);
}

}  // namespace sample
}  // namespace services
}  // namespace stan
#endif
