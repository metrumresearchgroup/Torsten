#ifndef STAN_MATH_TORSTEN_MPI_MY_WORKER_HPP
#define STAN_MATH_TORSTEN_MPI_MY_WORKER_HPP

namespace torsten {
  namespace mpi {

    /*
     * Round-Robin style allocation of individules to MPI
     * nodes.
     *
     * @param i ID of the individule
     * @param np size of the population
     * @param nproc number of processors
     * @return true if @c i is current node's responsiblity
     * for computation, false otherwise.
     */
    inline int my_worker(int i, int np, int nproc) {
      int n = np / nproc;
      int r = np % nproc;
      if(i < r * (n + 1)) {
        return i / (n + 1);
      } else {
        return r + (i - r * (n + 1)) / n;
      }
    }
  }
}

#endif
