#ifndef STAN_MATH_TORSTEN_MPI_ENVIONMENT_HPP
#define STAN_MATH_TORSTEN_MPI_ENVIONMENT_HPP

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
        if(!flag) MPI_Init(NULL, NULL);
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
