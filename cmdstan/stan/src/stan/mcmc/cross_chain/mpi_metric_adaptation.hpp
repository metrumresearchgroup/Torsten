#ifndef STAN_MCMC_MPI_METRIC_ADAPTATION_HPP
#define STAN_MCMC_MPI_METRIC_ADAPTATION_HPP

#include <stan/math/prim/fun/Eigen.hpp>

#ifdef STAN_LANG_MPI
#include <stan/math/torsten/mpi/session.hpp>
#endif

namespace stan {

namespace mcmc {

  class mpi_metric_adaptation {
  public:
    virtual void add_sample(const Eigen::VectorXd& q, int curr_win_count) {};

    virtual void restart() {}

#ifdef MPI_ADAPTED_WARMUP
    virtual void learn_metric(Eigen::VectorXd& var, int win, int curr_win_count,
                              const stan::math::mpi::Communicator& comm)
    {}
      

    virtual void learn_metric(Eigen::MatrixXd& covar, int win, int curr_win_count,
                              const stan::math::mpi::Communicator& comm)
    {}
#endif
  };

}  // namespace mcmc

}  // namespace stan

#endif
