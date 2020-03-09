#ifndef STAN_MATH_TORSTEN_MPI_SESSION_HPP
#define STAN_MATH_TORSTEN_MPI_SESSION_HPP

#ifdef TORSTEN_MPI

#include <stan/math/torsten/mpi/communicator.hpp>

/*
 * Given a Torsten communicator comm that has m workers and solves a group of
 * size n, the total nb. of processes n_proc involved is
 *
 * n_proc = m * b
 *
 * with b the nb. of cores used in the braid for a single subject/ode.
 * The population/group is evenly distributed
 * among m workers, with each worker employs b
 * processes/cores to solve a single subject/ode.
 * This implies there is effectively no braid(in-time-parallel) when b = 1.
 */

#define NUM_TORSTEN_COMM 4
#define TORSTEN_COMM_PMX_PARM 0
#define TORSTEN_COMM_PMX_DATA 1
#define TORSTEN_COMM_ODE_PARM 2
#define TORSTEN_COMM_ODE_DATA 3

namespace torsten {
  namespace mpi {
    /*
     * MPI communicator wrapper for RAII. Note that no
     * MPI's predfined comm sich as @c MPI_COMM_WOLRD are allowed.
     */
    template<int N_comm>
    struct Session {
      static torsten::mpi::Envionment env;
      static std::vector<torsten::mpi::Communicator> comms;
    };

    template<int N_comm>
    torsten::mpi::Envionment Session<N_comm>::env;

    template<int N_comm>
    std::vector<Communicator> Session<N_comm>::comms(N_comm, Communicator(Session<N_comm>::env, MPI_COMM_WORLD));
  }
}

#endif

#endif
