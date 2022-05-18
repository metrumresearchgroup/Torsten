#ifdef MPI_ADAPTED_WARMUP

#include <gtest/gtest.h>
#include <stan/mcmc/cross_chain/mpi_cross_chain_adapter.hpp>
#include <stan/analyze/mcmc/compute_potential_scale_reduction.hpp>
#include <stan/math/torsten/mpi.hpp>
#include <stan/callbacks/stream_logger.hpp>
#include <stan/mcmc/hmc/nuts/adapt_unit_e_nuts.hpp>
#include <boost/random/additive_combine.hpp>
#include <stan/callbacks/writer.hpp>
#include <stan/services/util/mcmc_writer.hpp>
#include <stan/services/util/generate_transitions.hpp>
#include <stan/io/dump.hpp>
#include <fstream>

#include <stan/services/sample/hmc_nuts_diag_e_adapt.hpp>
#include <stan/services/sample/hmc_nuts_dense_e_adapt.hpp>
#include <stan/io/empty_var_context.hpp>
#include <test/test-models/good/mcmc/hmc/common/gauss3D.hpp>
#include <test/unit/services/instrumented_callbacks.hpp>
#include <test/unit/services/check_adaptation.hpp>
#include <iostream>
#include <gtest/gtest.h>

TORSTEN_MPI_SESSION_INIT;

using Eigen::MatrixXd;
using Eigen::Matrix;
using std::vector;
using stan::math::mpi::Session;
using stan::math::mpi::Communicator;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using boost::accumulators::accumulator_set;
using boost::accumulators::stats;
using boost::accumulators::tag::mean;
using boost::accumulators::tag::variance;

class CrossChainAdapterTest : public testing::Test {
 public:
  CrossChainAdapterTest() :
    num_chains(3),
    max_num_windows(5),
    cross_chain_window_size(3),
    comm(Session::inter_chain_comm(num_chains)),
    num_iterations(cross_chain_window_size * max_num_windows),
    cross_chain_rhat(1.1),
    cross_chain_ess(100),
    num_warmup(100)
  {}

  std::stringstream model_log;
  stan::callbacks::logger logger;
  stan::callbacks::writer init, parameter, diagnostic;
  int num_chains;     // must run with mpiexec -n with n>=4
  int max_num_windows;
  int cross_chain_window_size;
  const Communicator& comm;
  int num_iterations;
  double cross_chain_rhat;
  double cross_chain_ess;
  int num_warmup;
};

TEST_F(CrossChainAdapterTest, stepsize) {
  const int n_par = 4;
  Eigen::VectorXd q = Eigen::VectorXd::Zero(n_par);

  stan::mcmc::mpi_cross_chain_adapter adapter;
  adapter.set_cross_chain_adaptation_params(75, 50, num_warmup,
                                            cross_chain_window_size, num_chains,
                                            cross_chain_rhat, cross_chain_ess);  
  stan::mcmc::mpi_var_adaptation var_adapt(n_par, num_warmup, cross_chain_window_size);
  adapter.set_cross_chain_metric_adaptation(&var_adapt);

  double chain_stepsize_base = 0.13;
  double chain_stepsize = chain_stepsize_base + comm.rank();
  double harmonic_mean_stepsize = 0.0;
  for (int rank = 0; rank < num_chains; ++rank) {
    harmonic_mean_stepsize += 1.0 / (chain_stepsize_base + rank);
  }
  harmonic_mean_stepsize = num_chains / harmonic_mean_stepsize;

  // double harmonic_mean_stepsize = 0.427066;
  adapter.add_cross_chain_sample(1.1, q);
  EXPECT_FLOAT_EQ(adapter.cross_chain_stepsize(chain_stepsize), chain_stepsize);
  adapter.add_cross_chain_sample(1.2, q);
  EXPECT_FLOAT_EQ(adapter.cross_chain_stepsize(chain_stepsize), chain_stepsize);
  adapter.add_cross_chain_sample(1.2, q);
  adapter.set_cross_chain_adapted(true);
  EXPECT_FLOAT_EQ(adapter.cross_chain_stepsize(chain_stepsize), harmonic_mean_stepsize);
}

TEST_F(CrossChainAdapterTest, gather) {
  const int n_par = 4;
  Eigen::VectorXd q = Eigen::VectorXd::Zero(n_par);
  if (comm.rank() >= 0) {       // only inter-chains
    q(comm.rank()) = comm.rank();    
  }

  stan::mcmc::mpi_cross_chain_adapter adapter;
  adapter.set_cross_chain_adaptation_params(75, 50, num_warmup,
                                            cross_chain_window_size, num_chains,
                                            cross_chain_rhat, cross_chain_ess);  
  stan::mcmc::mpi_var_adaptation var_adapt(n_par, num_warmup, cross_chain_window_size);
  adapter.set_cross_chain_metric_adaptation(&var_adapt);

  double lp;
  // window 1
  lp = 1.1 + comm.rank();
  adapter.add_cross_chain_sample(lp, q);
  lp = 3.8 + comm.rank();
  adapter.add_cross_chain_sample(lp, q);
  lp = 0.1 + comm.rank();
  adapter.add_cross_chain_sample(lp, q);
  Eigen::MatrixXd all_lp_draws;
  std::vector<double> all_chain_gather;
  adapter.cross_chain_gather(all_chain_gather, all_lp_draws);
  if (comm.rank() == 0) {
    EXPECT_EQ(all_chain_gather.size(), num_chains * 5);
    EXPECT_FLOAT_EQ(all_chain_gather[0], 5.0/3.0);
    EXPECT_FLOAT_EQ(all_chain_gather[1], 3.6633333333);
    EXPECT_FLOAT_EQ(all_chain_gather[2], 1.1);
    EXPECT_FLOAT_EQ(all_chain_gather[3], 3.8);
    EXPECT_FLOAT_EQ(all_chain_gather[4], 0.1);
    EXPECT_FLOAT_EQ(all_chain_gather[num_chains * 5 - 5], 5.0/3.0 + num_chains - 1);
    EXPECT_FLOAT_EQ(all_chain_gather[num_chains * 5 - 4], 3.6633333333);
    EXPECT_FLOAT_EQ(all_chain_gather[num_chains * 5 - 3], 1.1 + num_chains - 1);
    EXPECT_FLOAT_EQ(all_chain_gather[num_chains * 5 - 2], 3.8 + num_chains - 1);
    EXPECT_FLOAT_EQ(all_chain_gather[num_chains * 5 - 1], 0.1 + num_chains - 1);
  } else {
    EXPECT_EQ(all_chain_gather.size(), 0);
  }
  EXPECT_EQ(adapter.num_active_cross_chain_windows(), 1);

  // window 2
  lp = 2.1 + comm.rank();
  adapter.add_cross_chain_sample(lp, q);
  lp = 4.8 + comm.rank();
  adapter.add_cross_chain_sample(lp, q);

  // only gather at window end
  std::vector<double> empty_chain_gather;
  int n_gather = adapter.cross_chain_gather(empty_chain_gather, all_lp_draws);
  EXPECT_EQ(n_gather, 0);
  EXPECT_EQ(empty_chain_gather.size(), 0);

  lp = 1.1 + comm.rank();
  adapter.add_cross_chain_sample(lp, q);
  adapter.cross_chain_gather(all_chain_gather, all_lp_draws);
  if (comm.rank() == 0) {
    EXPECT_EQ(all_chain_gather.size(), num_chains * 7);
    EXPECT_FLOAT_EQ(all_chain_gather[0], 13.0/6.0);
    EXPECT_FLOAT_EQ(all_chain_gather[1], 3.23066666);
    EXPECT_FLOAT_EQ(all_chain_gather[2], 8.0/3.0);
    EXPECT_FLOAT_EQ(all_chain_gather[3], 3.6633333333);
    EXPECT_FLOAT_EQ(all_chain_gather[4], 2.1);
    EXPECT_FLOAT_EQ(all_chain_gather[5], 4.8);
    EXPECT_FLOAT_EQ(all_chain_gather[6], 1.1);
    EXPECT_FLOAT_EQ(all_chain_gather[num_chains * 7 - 7], 13.0/6.0 + num_chains - 1);
    EXPECT_FLOAT_EQ(all_chain_gather[num_chains * 7 - 6], 3.23066666);
    EXPECT_FLOAT_EQ(all_chain_gather[num_chains * 7 - 5], 8.0/3.0 + num_chains - 1);
    EXPECT_FLOAT_EQ(all_chain_gather[num_chains * 7 - 4], 3.6633333333);
    EXPECT_FLOAT_EQ(all_chain_gather[num_chains * 7 - 3], 2.1 + num_chains - 1);
    EXPECT_FLOAT_EQ(all_chain_gather[num_chains * 7 - 2], 4.8 + num_chains - 1);
    EXPECT_FLOAT_EQ(all_chain_gather[num_chains * 7 - 1], 1.1 + num_chains - 1);
  } else {
    EXPECT_EQ(all_chain_gather.size(), 0);
  }
  EXPECT_EQ(adapter.num_active_cross_chain_windows(), 2);
}

TEST_F(CrossChainAdapterTest, negative_adapted_window) {
  const int n_par = 4;
  Eigen::VectorXd q = Eigen::VectorXd::Zero(n_par);
  if (comm.rank() >= 0) {       // only inter-chains
    q(comm.rank()) = comm.rank();    
  }

  stan::mcmc::mpi_cross_chain_adapter adapter;
  adapter.set_cross_chain_adaptation_params(75, 50, num_warmup,
                                            cross_chain_window_size, num_chains,
                                            cross_chain_rhat, cross_chain_ess);  
  stan::mcmc::mpi_var_adaptation var_adapt(n_par, num_warmup, cross_chain_window_size);
  adapter.set_cross_chain_metric_adaptation(&var_adapt);

  std::vector<double> all_chain_gather;
  Eigen::MatrixXd all_lp_draws;

  double lp;
  // window 1
  lp = 1.1 + comm.rank();
  adapter.add_cross_chain_sample(lp, q);
  lp = 3.8 + comm.rank();
  adapter.add_cross_chain_sample(lp, q);
  lp = 0.1 + comm.rank();
  adapter.add_cross_chain_sample(lp, q);
  int n_gather = adapter.cross_chain_gather(all_chain_gather, all_lp_draws);

  // window 2
  lp = 2.1 + comm.rank();
  adapter.add_cross_chain_sample(lp, q);
  lp = 4.8 + comm.rank();
  adapter.add_cross_chain_sample(lp, q);
  lp = 1.1 + comm.rank();
  adapter.add_cross_chain_sample(lp, q);
  n_gather = adapter.cross_chain_gather(all_chain_gather, all_lp_draws);

  // standard rhat & ess
  std::vector<std::vector<double>> chain_lp(num_chains);
  for (int chain = 0; chain < num_chains; ++chain) {
    chain_lp[chain].push_back(1.1 + chain);
    chain_lp[chain].push_back(3.8 + chain);
    chain_lp[chain].push_back(0.1 + chain);
    chain_lp[chain].push_back(2.1 + chain);
    chain_lp[chain].push_back(4.8 + chain);
    chain_lp[chain].push_back(1.1 + chain);
  }
  std::vector<const double*> draws(num_chains);
  for (int chain = 0; chain < num_chains; ++chain) {
    draws[chain] = chain_lp[chain].data();
  }

  // compare window 1 of total 2 windows against standard results
  {
    double ess =
      stan::analyze::compute_effective_sample_size(draws, 2 * cross_chain_window_size);
    double rhat =
      stan::analyze::compute_potential_scale_reduction(draws, 2 * cross_chain_window_size);

    // compare window 2
    stan::callbacks::logger logger;
    if (stan::math::mpi::Session::is_in_inter_chain_comm(num_chains)) {
      int adapted_win = adapter.cross_chain_adapted_window(n_gather,
                                                           all_chain_gather,
                                                           all_lp_draws,
                                                           logger);
      if (comm.rank() == 0) {
        EXPECT_FLOAT_EQ(adapter.cross_chain_adapt_rhat()[0], rhat);
        EXPECT_FLOAT_EQ(adapter.cross_chain_adapt_ess()[0], ess);
      }
      EXPECT_EQ(adapted_win, -1);
    }
  }

  // window 3
  lp = 2.5 + comm.rank();
  adapter.add_cross_chain_sample(lp, q);
  lp = 4.5 + comm.rank();
  adapter.add_cross_chain_sample(lp, q);
  lp = 1.5 + comm.rank();
  adapter.add_cross_chain_sample(lp, q);
  n_gather = adapter.cross_chain_gather(all_chain_gather, all_lp_draws);

  for (int chain = 0; chain < num_chains; ++chain) {
    chain_lp[chain].push_back(2.5 + chain);
    chain_lp[chain].push_back(4.5 + chain);
    chain_lp[chain].push_back(1.5 + chain);
  }
  for (int chain = 0; chain < num_chains; ++chain) {
    draws[chain] = chain_lp[chain].data();
  }

  // compare window 1 of total 3 windows against standard results
  {
    double ess =
      stan::analyze::compute_effective_sample_size(draws, 3 * cross_chain_window_size);
    double rhat =
      stan::analyze::compute_potential_scale_reduction(draws, 3 * cross_chain_window_size);

    if (stan::math::mpi::Session::is_in_inter_chain_comm(num_chains)) {
      int adapted_win = adapter.cross_chain_adapted_window(n_gather,
                                                           all_chain_gather,
                                                           all_lp_draws,
                                                           logger);
      if (comm.rank() == 0) {
        EXPECT_FLOAT_EQ(adapter.cross_chain_adapt_rhat()[0], rhat);
        EXPECT_FLOAT_EQ(adapter.cross_chain_adapt_ess()[0], ess);
      }
      EXPECT_EQ(adapted_win, -2);
    }    
  }

  // compare window 2 of total 3 windows against standard results
  // for standard computation we need recollect draws
  {
    std::vector<std::vector<double>> chain_lp(num_chains);
    for (int chain = 0; chain < num_chains; ++chain) {
      chain_lp[chain].push_back(2.1 + chain);
      chain_lp[chain].push_back(4.8 + chain);
      chain_lp[chain].push_back(1.1 + chain);
      chain_lp[chain].push_back(2.5 + chain);
      chain_lp[chain].push_back(4.5 + chain);
      chain_lp[chain].push_back(1.5 + chain);
    }
    std::vector<const double*> draws(num_chains);
    for (int chain = 0; chain < num_chains; ++chain) {
      draws[chain] = chain_lp[chain].data();
    }

    double ess =
      stan::analyze::compute_effective_sample_size(draws, 2 * cross_chain_window_size);
    double rhat =
      stan::analyze::compute_potential_scale_reduction(draws, 2 * cross_chain_window_size);

    if (stan::math::mpi::Session::is_in_inter_chain_comm(num_chains)) {
      int adapted_win = adapter.cross_chain_adapted_window(n_gather,
                                                           all_chain_gather,
                                                           all_lp_draws,
                                                           logger);
      if (comm.rank() == 0) {
        EXPECT_FLOAT_EQ(adapter.cross_chain_adapt_rhat()[1], rhat);
        EXPECT_FLOAT_EQ(adapter.cross_chain_adapt_ess()[1], ess);
      }
      EXPECT_EQ(adapted_win, -2);
    }
  }
}

TEST_F(CrossChainAdapterTest, positive_adapted_window) {
  const int n_par = 4;
  Eigen::VectorXd q = Eigen::VectorXd::Zero(n_par);
  if (comm.rank() >= 0) {       // only inter-chains
    q(comm.rank()) = comm.rank();    
  }

  cross_chain_rhat = 1.12;
  cross_chain_ess = 15.0;
  stan::mcmc::mpi_cross_chain_adapter adapter;
  adapter.set_cross_chain_adaptation_params(75, 50, num_warmup,
                                            cross_chain_window_size, num_chains,
                                            cross_chain_rhat, cross_chain_ess);
  stan::mcmc::mpi_var_adaptation var_adapt(n_par, num_warmup, cross_chain_window_size);
  adapter.set_cross_chain_metric_adaptation(&var_adapt);

  std::vector<double> all_chain_gather;
  Eigen::MatrixXd all_lp_draws;
  stan::callbacks::logger logger;

  double lp;
  int n_gather = 0;
  // window 1
  lp = 1.1 + comm.rank();
  adapter.add_cross_chain_sample(lp, q);
  lp = 3.8 + comm.rank();
  adapter.add_cross_chain_sample(lp, q);
  lp = 0.1 + comm.rank();
  adapter.add_cross_chain_sample(lp, q);
  n_gather = adapter.cross_chain_gather(all_chain_gather, all_lp_draws);

  // window 2
  lp = 2.1 + comm.rank();
  adapter.add_cross_chain_sample(lp, q);
  lp = 4.8 + comm.rank();
  adapter.add_cross_chain_sample(lp, q);
  lp = 1.1 + comm.rank();
  adapter.add_cross_chain_sample(lp, q);
  n_gather = adapter.cross_chain_gather(all_chain_gather, all_lp_draws);

  if (stan::math::mpi::Session::is_in_inter_chain_comm(num_chains)) {
    int adapted_win = adapter.cross_chain_adapted_window(n_gather,
                                                         all_chain_gather,
                                                         all_lp_draws,
                                                         logger);
    EXPECT_EQ(adapted_win, 0);
  }

  // window 3
  lp = 2.5 + comm.rank();
  adapter.add_cross_chain_sample(lp, q);
  lp = 4.5 + comm.rank();
  adapter.add_cross_chain_sample(lp, q);
  lp = 1.5 + comm.rank();
  adapter.add_cross_chain_sample(lp, q);
  n_gather = adapter.cross_chain_gather(all_chain_gather, all_lp_draws);

  if (stan::math::mpi::Session::is_in_inter_chain_comm(num_chains)) {
    int adapted_win = adapter.cross_chain_adapted_window(n_gather,
                                                         all_chain_gather,
                                                         all_lp_draws,
                                                         logger);
    EXPECT_EQ(adapted_win, 1);
  }    
}

#endif
