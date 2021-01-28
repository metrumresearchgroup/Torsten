#ifndef STAN_MATH_TORSTEN_MPI_CVODES_BRAID_HPP
#define STAN_MATH_TORSTEN_MPI_CVODES_BRAID_HPP

#include <braid.hpp>
#include <stan/math/torsten/dsolve/sundials_check.hpp>
#include <nvector/nvector_serial.h>
#include <cvodes/cvodes.h>
#include <vector>

namespace torsten {
  namespace mpi {
    struct CVBraidVec {
      N_Vector y;

      /*
       * constructor
       */
      CVBraidVec(int n) : y(N_VNew_Serial(n)) {
        N_VConst(0.0, y);
      }

      /*
       * constructor
       */
      CVBraidVec(const N_Vector& y_in) : CVBraidVec(NV_LENGTH_S(y_in))
      {
        N_VScale(1.0, y_in, y);
      }

      /*
       * copy constructor
       */
      CVBraidVec(const CVBraidVec& other) : CVBraidVec(other.y)
      {}

      ~CVBraidVec() {
        N_VDestroy_Serial(y);
      }
    };

    /*
     * braid_Vector == CVBraidVec*
     */
    struct CVBraidApp : BraidApp {
      void* mem;
      void* user_data;
      const N_Vector& y0;
      const int n;
      const int buf_size;
      CVBraidVec y0_v;
      CVBraidVec y_v;
      int rank;
      std::vector<MPI_Request> req;
      std::vector<std::vector<double> > result;

      /*
       * The constructor also initializes master to receive results
       * from the others.
       */
      CVBraidApp(void* mem_in, void* user_data_in,
                 const N_Vector& y0_in, MPI_Comm comm_in, double t0, double t1, int nt)
        : BraidApp(comm_in, t0, t1, nt),
          mem(mem_in), user_data(user_data_in),
          y0(y0_in), n(NV_LENGTH_S(y0_in)), buf_size(sizeof(double) * n),
          y0_v(y0), y_v(n), rank(-1), req(nt)
      {
        MPI_Comm_rank(comm_t, &rank);
        if (rank == 0) {
          result.resize(nt);
          for (int it = 0; it < nt; ++it) {
            result[it] = std::vector<double>(n, 0.0);
            MPI_Irecv(result[it].data(), n, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, comm_t, &req[it]);
          }
        }
      }

      virtual ~CVBraidApp() {}

      virtual int Step(braid_Vector u, braid_Vector ustop, braid_Vector fstop, BraidStepStatus &status)
      {
        double t0;             /* current time */
        double t1;              /* evolve to this time*/
        int braid_level;
        status.GetTstartTstop(&t0, &t1);
        status.GetLevel(&braid_level);
        CVBraidVec* vec = (CVBraidVec*) u;

        CHECK_SUNDIALS_CALL(CVodeReInit(mem, t0, vec -> y));
        CHECK_SUNDIALS_CALL(CVodeSStolerances(mem, 1.e-6, 1.e-6));
        CHECK_SUNDIALS_CALL(CVodeSetUserData(mem, user_data));
        CHECK_SUNDIALS_CALL(CVodeSetMaxNumSteps(mem, 10000));
        CHECK_SUNDIALS_CALL(CVodeSetMaxErrTestFails(mem, 20));
        CHECK_SUNDIALS_CALL(CVodeSetMaxConvFails(mem, 30));

        CHECK_SUNDIALS_CALL(CVode(mem, t1, vec -> y, &t0, CV_NORMAL));

        return 0;
      }

      virtual int Clone(braid_Vector u, braid_Vector *v_ptr)
      {
        CVBraidVec* vec = (CVBraidVec*) u;
        *v_ptr = (braid_Vector) new CVBraidVec(vec -> y);
        return 0;
      }

      virtual int Init(double t, braid_Vector *u_ptr)
      {
        if (t == tstart) {
          *u_ptr = (braid_Vector) new CVBraidVec(y0);
        } else {
          *u_ptr = (braid_Vector) new CVBraidVec(n);
        }
        return 0;
      }

      virtual int Free(braid_Vector u)
      {
        CVBraidVec* vec = (CVBraidVec*) u;
        N_VDestroy_Serial(vec -> y);

        return 0;    
      }

      virtual int Sum(double a, braid_Vector x, double b, braid_Vector y)
      {
        CVBraidVec* xp = (CVBraidVec*) x;
        CVBraidVec* yp = (CVBraidVec*) y;
        N_VLinearSum(a, xp -> y, b, yp -> y, yp -> y);
        return 0;
      }

      virtual int SpatialNorm(braid_Vector u, double *norm_ptr) {
        const CVBraidVec* const up = (CVBraidVec*) u;
        *norm_ptr = N_VWrmsNorm(up -> y, up -> y);
        return 0;     
      }

      virtual int BufSize(int *size_ptr, BraidBufferStatus &status)
      {
        *size_ptr = buf_size;
        return 0;
      }

      virtual int BufPack(braid_Vector u, void *buffer, BraidBufferStatus &status)
      {
        try {
          const CVBraidVec* const up = (CVBraidVec*) u;
          double *dbuf = (double *) buffer;
          N_Vector v = N_VMake_Serial(n, dbuf);
          N_VScale(1.0, up -> y, v);
        } catch (...) {
          throw;
        }

        return 0;
      }

      virtual int BufUnpack(void *buffer, braid_Vector *u_ptr, BraidBufferStatus &status)
      {
        try {
          *u_ptr = (braid_Vector) new CVBraidVec(n);
          CVBraidVec* const up = (CVBraidVec*) *u_ptr;
          double *dbuf = (double *) buffer;
          N_Vector v = N_VMake_Serial(n, dbuf);
          N_VScale(1.0, v, up -> y);
        } catch (...) {
          throw;
        }

        return 0;    
      }

      virtual int Coarsen(braid_Vector fu, braid_Vector *cu_ptr, BraidCoarsenRefStatus &status)
      {
        return 0;
      }

      virtual int Refine(braid_Vector cu, braid_Vector *fu_ptr, BraidCoarsenRefStatus &status)
      {
        return 0;    
      }

      /*
       * Assume access level = 1 
       */
      virtual int Access(braid_Vector u, BraidAccessStatus &astatus)
      {
        const CVBraidVec* const up = (CVBraidVec*) u;
        int it;
        astatus.GetTIndex(&it);
        if (it > 0) {
          MPI_Send(NV_DATA_S(up -> y), n, MPI_DOUBLE, 0, it - 1, comm_t);
        }
        return 0;
      }

      virtual int Residual(braid_Vector u, braid_Vector r, BraidStepStatus &pstatus)
      {
        return 0;
      }

      virtual std::vector<std::vector<double> > get_result() {
        if (rank == 0) {
          int nt = ntime;
          std::vector<std::vector<double> > res;
          res.resize(nt);
          for (int it = 0; it < nt; ++it) {
            MPI_Status status; 
            MPI_Wait(&req[it], &status);
            res[status.MPI_TAG] = result[it];
          }
          return res;
        } else {
          throw std::runtime_error("braid result collection is performed by non-zero rank");
        }
      }

    };

  }
}

#endif
