#ifndef TORSTEN_MPI_SESSION_INIT_HPP
#define TORSTEN_MPI_SESSION_INIT_HPP

/**
 * Define static variables of <code>Session</code> \
 * Only to be included in tests & cmdstan top level #<code>command.hpp</code> \
 */
#if defined(STAN_LANG_MPI) || defined(TORSTEN_MPI)
#define TORSTEN_MPI_SESSION_INIT \
 stan::math::mpi::Envionment   stan::math::mpi::Session::env; \
                      MPI_Comm stan::math::mpi::Session::MPI_COMM_INTER_CHAIN(MPI_COMM_NULL); \
                      MPI_Comm stan::math::mpi::Session::MPI_COMM_INTRA_CHAIN(MPI_COMM_NULL); \
 stan::math::mpi::Communicator stan::math::mpi::Session::inter_chain(MPI_COMM_NULL); \
 stan::math::mpi::Communicator stan::math::mpi::Session::intra_chain(MPI_COMM_NULL); \
 int torsten::mpi::Session::num_chains = 1;
#else
#define TORSTEN_MPI_SESSION_INIT
#endif

#endif
