#ifndef STAN_MATH_TORSTEN_MPI_ENVIONMENT_HPP
#define STAN_MATH_TORSTEN_MPI_ENVIONMENT_HPP

#include <stan/math/prim/scal/err/check_greater_or_equal.hpp>
#ifdef TORSTEN_MPI
#include <boost/mpi.hpp>
#endif

namespace torsten {
  namespace mpi {
    struct Envionment {
      Envionment() {
        init();
      }
      ~Envionment() {
        finalize();
      }

      static void init() {
#ifdef TORSTEN_MPI
        int flag;
        MPI_Initialized(&flag);
        if(!flag) {
#ifdef _OPENMP
          static const char* caller("torsten::mpi::init");
          int provided;
          MPI_Init_thread(NULL, NULL, MPI_THREAD_FUNNELED, &provided);
          stan::math::check_greater_or_equal(caller, "level of threads support", provided, MPI_THREAD_FUNNELED); // NOLINT
#else
          MPI_Init(NULL, NULL);
#endif
        }
#endif
      }

      static void finalize() {
#ifdef TORSTEN_MPI
        int flag;
        MPI_Finalized(&flag);
        if(!flag) MPI_Finalize();
#endif
      }
    };
  }
}

#endif
