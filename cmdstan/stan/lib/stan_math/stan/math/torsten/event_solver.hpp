#ifndef STAN_MATH_TORSTEN_PKMODEL_REFACTOR_PRED_HPP
#define STAN_MATH_TORSTEN_PKMODEL_REFACTOR_PRED_HPP

#include <stan/math/torsten/dsolve/pk_vars.hpp>
#include <stan/math/torsten/events_manager.hpp>
#include <stan/math/torsten/mpi/session.hpp>
#include <stan/math/torsten/mpi/precomputed_gradients.hpp>
#include <Eigen/Dense>
#include <vector>

namespace torsten{
  /*
   * the solver wrapper is aware of @c T_model so it build model
   * accordingly.
   */
  template<typename T_model, typename... T_pred>
  struct EventSolver {

    /*
     * Data used to fill the results when computation throws exception.
     */
    static constexpr double invalid_res_d = std::numeric_limits<double>::quiet_NaN();

    /*
     * calculate the number of parameters according to the
     * model type.
     *
     * @tparam T_er type of events record
     * @param id subject id among a population
     * @param rec events record for the population
     */
    template<typename T_er>
    static int system_size(int id, const T_er& rec) {
      const int ncmt = rec.ncmt;
      const int npar = rec.parameter_size();
      const int nvar = T_model::nvars(ncmt, npar);
      const int nvar_ss = T_model::template nvars<typename T_er::T_amt, typename T_er::T_par_rate, typename T_er::T_par_ii>(npar); //NOLINT
      return rec.has_ss_dosing(id) ? pk_nsys(ncmt, nvar, nvar_ss) : pk_nsys(ncmt, nvar);
    }

    /**
     * Every Torsten function calls Pred.
     *
     * Predicts the amount in each compartment for each event,
     * given the event schedule and the parameters of the model.
     *
     * Proceeds in two steps. First, computes all the events that
     * are not included in the original data set but during which
     * amounts in the system get updated. Secondly, predicts
     * the amounts in each compartment sequentially by going
     * through the augmented schedule of events. The returned pred
     * Matrix only contains the amounts in the event originally
     * specified by the users.
     *
     * This function is valid for all models. What changes from one
     * model to the other are the Pred1 and PredSS functions, which
     * calculate the amount at an individual event.
     *
     * @tparam T_em type the @c EventsManager
     * @param[in] time times of events
     * @param[in] amt amount at each event
     * @param[in] rate rate at each event
     * @param[in] ii inter-dose interval at each event
     * @param[in] evid event identity:
     *                    (0) observation
     *                    (1) dosing
     *                    (2) other
     *                    (3) reset
     *                    (4) reset AND dosing
     * @param[in] cmt compartment number at each event (starts at 1)
     * @param[in] addl additional dosing at each event
     * @param[in] ss steady state approximation at each event
     * (0: no, 1: yes)
     * @param[in] pMatrix parameters at each event
     * @param[in] addParm additional parameters at each event
     * @parem[in] model basic info for ODE model and evolution operators
     * @param[in] SystemODE matrix describing linear ODE system that
     * defines compartment model. Used for matrix exponential solutions.
     * Included because it may get updated in modelParameters.
     * @return a matrix with predicted amount in each compartment
     * at each event.
     */
    template<typename T_events_record, typename... Ts>
    void pred(int id,
              const T_events_record& events_rec,
              Eigen::Matrix<typename EventsManager<T_events_record>::T_scalar, -1, -1>& res,
              const T_pred... pred_pars,
              const Ts... model_pars) {
      using Eigen::Matrix;
      using Eigen::Dynamic;
      using std::vector;
      using::stan::math::multiply;
      using refactor::PKRec;
      using EM = EventsManager<T_events_record>;

      using scalar = typename EM::T_scalar;

      int nCmt = EM::nCmt(events_rec);

      res.resize(nCmt, events_rec.num_event_times(id));
      PKRec<scalar> init(nCmt);
      init.setZero();

      try {
        EM em(id, events_rec);
        auto events = em.events();
        int ikeep = 0, iev = 0;
        while(ikeep < em.nKeep) {
          stepper(iev, init, em, pred_pars..., model_pars...);
          if(events.keep(iev)) {
            res.col(ikeep) = init;
            ikeep++;
          }
          iev++;
        }
      } catch (const std::exception& e) {
        throw;
      }
    }

    /*
     * Step through a range of events.
     */
    template<typename T_em, typename... Ts>
    void stepper(int i, refactor::PKRec<typename T_em::T_scalar>& init,
                        const T_em& em, const T_pred... pred_pars, const Ts... model_pars) {
      auto events = em.events();

      using scalar = typename T_em::T_scalar;
      typename T_em::T_time tprev = i == 0 ? events.time(0) : events.time(i-1);

      Eigen::Matrix<scalar, -1, 1> pred1;

      if (events.is_reset(i)) {
        init.setZero();
      } else if (events.is_ss_dosing(i)) {  // steady state event
        typename T_em::T_time model_time = events.time(i); // FIXME: time is not t0 but for adjust within SS solver
        std::vector<typename T_em::T_rate> dummy_rate;
        T_model pkmodel {model_time, init, dummy_rate, events.model_param(i), model_pars...};
        auto curr_amt = events.fractioned_amt(i);
        pred1 = stan::math::multiply(pkmodel.solve(curr_amt,
                                                   events.rate(i),
                                                   events.ii(i),
                                                   events.cmt(i),
                                                   pred_pars...),
                                     scalar(1.0));

        if (events.ss(i) == 2)
          init += pred1;  // steady state without reset
        else
          init = pred1;  // steady state with reset (ss = 1)
      } else if (events.time(i) > tprev) {           // non-steady dosing event
        typename T_em::T_time model_time = tprev;
        auto curr_rates = events.fractioned_rates(i);
        T_model pkmodel {model_time, init, curr_rates, events.model_param(i), model_pars...};
        pred1 = pkmodel.solve(events.time(i), pred_pars...);
        init = pred1;
      }

      if (events.is_bolus_dosing(i)) {
        init(0, events.cmt(i) - 1) += events.fractioned_amt(i);
      }
      tprev = events.time(i);
    }

    template<typename T_em, typename... Ts>
    void stepper_solve(int i, refactor::PKRec<typename T_em::T_scalar>& init,
                        refactor::PKRec<double>& sol_d,
                        const T_em& em, const T_pred... pred_pars, const Ts... model_pars) {
      using std::vector;
      using stan::math::var;

      auto events = em.events();

      typename T_em::T_time tprev = i == 0 ? events.time(0) : events.time(i-1);

      if (events.is_reset(i)) {
        init.setZero();
      } else if (events.is_ss_dosing(i)) {  // steady state event
        typename T_em::T_time model_time = events.time(i);
        std::vector<typename T_em::T_rate> dummy_rate;
        T_model pkmodel {model_time, init, dummy_rate, events.model_param(i), model_pars...};
        auto curr_amt = events.fractioned_amt(i);
        vector<var> v_i = pkmodel.vars(curr_amt, events.rate(i), events.ii(i));
        sol_d = pkmodel.solve_d(curr_amt, events.rate(i), events.ii(i), events.cmt(i), pred_pars...);
        if (events.ss(i) == 2)
          init += torsten::mpi::precomputed_gradients(sol_d, v_i);  // steady state without reset
        else
          init = torsten::mpi::precomputed_gradients(sol_d, v_i);  // steady state with reset (ss = 1)
      } else if (events.time(i) > tprev) {
          typename T_em::T_time model_time = tprev;
          auto curr_rates = events.fractioned_rates(i);
          T_model pkmodel {model_time, init, curr_rates, events.model_param(i), model_pars...};
          vector<var> v_i = pkmodel.vars(events.time(i));
          sol_d = pkmodel.solve_d(events.time(i), pred_pars...);
          init = torsten::mpi::precomputed_gradients(sol_d, v_i);
      }

      if (events.is_bolus_dosing(i)) {
        init(0, events.cmt(i) - 1) += events.fractioned_amt(i);
      }
      tprev = events.time(i);
    }

    template<typename T_em, typename... Ts>
    void stepper_sync(int i, refactor::PKRec<typename T_em::T_scalar>& init,
                      refactor::PKRec<double>& sol_d,
                      const T_em& em, const T_pred... pred_pars, const Ts... model_pars) {
      using std::vector;
      using stan::math::var;

      auto events = em.events();

      typename T_em::T_time tprev = i == 0 ? events.time(0) : events.time(i-1);

      if (events.is_reset(i)) {
        init.setZero();
      } else if (events.is_ss_dosing(i)) {  // steady state event
        typename T_em::T_time model_time = events.time(i);
        std::vector<typename T_em::T_rate> dummy_rate;
        T_model pkmodel {model_time, init, dummy_rate, events.model_param(i), model_pars...};
        auto curr_amt = events.fractioned_amt(i);
        vector<var> v_i = pkmodel.vars(curr_amt, events.rate(i), events.ii(i));
        int nsys = torsten::pk_nsys(em.ncmt, v_i.size());
        if (events.ss(i) == 2)
          init += torsten::mpi::precomputed_gradients(sol_d.segment(0, nsys), v_i);  // steady state without reset
        else
          init = torsten::mpi::precomputed_gradients(sol_d.segment(0, nsys), v_i);  // steady state with reset (ss = 1)
      } else if (events.time(i) > tprev) {
          typename T_em::T_time model_time = tprev;
          auto curr_rates = events.fractioned_rates(i);
          T_model pkmodel {model_time, init, curr_rates, events.model_param(i), model_pars...};
          vector<var> v_i = pkmodel.vars(events.time(i));
          int nsys = torsten::pk_nsys(em.ncmt, v_i.size());
          init = torsten::mpi::precomputed_gradients(sol_d.segment(0, nsys), v_i);
      }

      if (events.is_bolus_dosing(i)) {
        init(0, events.cmt(i) - 1) += events.fractioned_amt(i);
      }
    }

#ifdef TORSTEN_MPI

    /*
     * MPI solution when the population
     * information passed in as ragged arrays.
     *
     */
    template<typename T_events_record, typename... Ts,
             typename std::enable_if_t<stan::is_var<typename EventsManager<T_events_record>::T_scalar>::value >* = nullptr> //NOLINT
    void pred(const T_events_record& events_rec,
              Eigen::Matrix<typename EventsManager<T_events_record>::T_scalar, -1, -1>& res,
              const T_pred... pred_pars,
              const Ts... model_pars) {
      using Eigen::Matrix;
      using Eigen::MatrixXd;
      using Eigen::VectorXd;
      using Eigen::Dynamic;
      using std::vector;
      using::stan::math::var;
      using::stan::math::multiply;
      using refactor::PKRec;

      using ER = T_events_record;
      using EM = EventsManager<ER>;
      using scalar = typename ER::T_scalar;

      const int nCmt = EM::nCmt(events_rec);
      const int np = events_rec.num_subjects();
      bool is_invalid = false;
      std::ostringstream rank_fail_msg;

      MPI_Comm comm;
      comm = torsten::mpi::Session<NUM_TORSTEN_COMM>::comms[TORSTEN_COMM_PMX_PARM];
      int rank = torsten::mpi::Session<NUM_TORSTEN_COMM>::comms[TORSTEN_COMM_PMX_PARM].rank;
      int size = torsten::mpi::Session<NUM_TORSTEN_COMM>::comms[TORSTEN_COMM_PMX_PARM].size;

      std::vector<MPI_Request> req(np);
      vector<MatrixXd> res_d(np);
      
      res.resize(nCmt, events_rec.total_num_event_times);

      PKRec<scalar> init(nCmt);
      PKRec<double> pred1;
      for (int id = 0; id < np; ++id) {

        /* For every rank */

        const int nKeep = events_rec.num_event_times(id);

        int nev = EM::num_events(id, events_rec);
        res_d[id].resize(system_size(id, events_rec), nev);
        res_d[id].setConstant(0.0);

        int my_worker_id = torsten::mpi::my_worker(id, np, size);

        /* only solver rank */

        if (rank == my_worker_id) {
          if (is_invalid) {
            res_d[id].setConstant(invalid_res_d);
          } else {
            try {
              EM em(id, events_rec);
              auto events = em.events();
              assert(nev == events.num_state_times());
              assert(nKeep == em.nKeep);
              init.setZero();

              int ikeep = 0, iev = 0;
              while(ikeep < em.nKeep) {
                stepper_solve(iev, init, pred1, em, pred_pars..., model_pars...);
                res_d[id].col(iev).segment(0, pred1.size()) = pred1;
                if(events.keep(iev)) {
                  res.col(EM::begin(id, events_rec) + ikeep) = init;
                  ikeep++;
                }
                iev++;
              }
            } catch (const std::exception& e) {
              is_invalid = true;
              res_d[id].setConstant(invalid_res_d);
              rank_fail_msg << "Rank " << rank << " failed to solve id " << id << ": " << e.what();
            }
          }
        }
        MPI_Ibcast(res_d[id].data(), res_d[id].size(), MPI_DOUBLE, my_worker_id, comm, &req[id]);
      }

      // int finished = 0;
      // int index;
      // int flag = 0;
      // for(int id = 0; id < np; ++id) {
      //   MPI_Wait(&req[id], MPI_STATUS_IGNORE);

      //   if (is_invalid) continue;
      //   if (std::isnan(res_d[id](0))) {
      //     assert(rank != torsten::mpi::my_worker(id, np, size));
      //     is_invalid = true;
      //     rank_fail_msg << "Rank " << rank << " received invalid data for id " << id;
      //   } else {
      //     EM em(id, events_rec);
      //     auto events = em.events();
      //     PKRec<scalar> init(nCmt); init.setZero();
      //     PKRec<double> pred1 = VectorXd::Zero(res_d[id].rows());

      //     int ikeep = 0, iev = 0;
      //     while(ikeep < em.nKeep) {
      //       pred1 = res_d[id].col(iev);
      //       stepper_sync(iev, init, pred1, em, pred_pars..., model_pars...);
      //       if(events.keep(iev)) {
      //         res.col(EM::begin(id, events_rec) + ikeep) = init;
      //         ikeep++;
      //       }
      //       iev++;
      //     }
      //   }
      // }

      int finished = 0;
      int index = 0;
      int flag = 0;
      while (finished < np) {
        MPI_Test(&req[index], &flag, MPI_STATUS_IGNORE);
        if (flag) {
          int id = index;
          index++;
          finished++;
          if (is_invalid) continue;
          if (std::isnan(res_d[id](0))) {
            assert(rank != torsten::mpi::my_worker(id, np, size));
            is_invalid = true;
            rank_fail_msg << "Rank " << rank << " received invalid data for id " << id;
          } else {
            EM em(id, events_rec);
            auto events = em.events();
            PKRec<scalar> init(nCmt); init.setZero();
            PKRec<double> pred1 = VectorXd::Zero(res_d[id].rows());

            int ikeep = 0, iev = 0;
            while(ikeep < em.nKeep) {
              pred1 = res_d[id].col(iev);
              stepper_sync(iev, init, pred1, em, pred_pars..., model_pars...);
              if(events.keep(iev)) {
                res.col(EM::begin(id, events_rec) + ikeep) = init;
                ikeep++;
              }
              iev++;
            }
          }
        } 
      }

      MPI_Barrier(comm);

      if(is_invalid) {
        throw std::runtime_error(rank_fail_msg.str());
      }
    }

    /*
     * Data-only MPI solver that takes ragged arrays as input.
     */
    template<typename T_events_record, typename... Ts>
    void pred(const T_events_record& events_rec, Eigen::MatrixXd& res,
                     const T_pred... pred_pars,
                     const Ts... model_pars) {
      using Eigen::Matrix;
      using Eigen::MatrixXd;
      using Eigen::VectorXd;
      using Eigen::Dynamic;
      using std::vector;
      using::stan::math::var;
      using::stan::math::multiply;
      using refactor::PKRec;

      using ER = NONMENEventsRecord<double, double, double, double, std::vector<double>, double, double>;
      using EM = EventsManager<ER>;

      const int nCmt = EM::nCmt(events_rec);
      const int np = events_rec.num_subjects();
      bool is_invalid = false;
      std::ostringstream rank_fail_msg;

      MPI_Comm comm;
      comm = torsten::mpi::Session<NUM_TORSTEN_COMM>::comms[TORSTEN_COMM_PMX_DATA];
      int rank = torsten::mpi::Session<NUM_TORSTEN_COMM>::comms[TORSTEN_COMM_PMX_DATA].rank;
      int size = torsten::mpi::Session<NUM_TORSTEN_COMM>::comms[TORSTEN_COMM_PMX_DATA].size;

      std::vector<MPI_Request> req(np);

      res.resize(nCmt, events_rec.total_num_event_times);

      PKRec<double> init(nCmt);
      for (int id = 0; id < np; ++id) {

        /* For every rank */

        int nKeep = events_rec.num_event_times(id);
        int my_worker_id = torsten::mpi::my_worker(id, np, size);
        int begin_id = EM::begin(id, events_rec) * nCmt;
        int size_id = nKeep * nCmt;

        /* only solver rank */

        if (rank == my_worker_id) {
          try {
            EM em(id, events_rec);
            auto events = em.events();
            init.setZero();
            int ikeep = 0, iev = 0;
            while(ikeep < em.nKeep) {
              stepper(iev, init, em, pred_pars..., model_pars...);
              if(events.keep(iev)) {
                res.col(EM::begin(id, events_rec) + ikeep) = init;
                ikeep++;
              }
              iev++;
            }
          } catch (const std::exception& e) {
            is_invalid = true;
            res(begin_id) = invalid_res_d;
            rank_fail_msg << "Rank " << rank << " failed to solve id " << id << ": " << e.what();
          }
        }

        MPI_Ibcast(res.data() + begin_id, size_id, MPI_DOUBLE, my_worker_id, comm, &req[id]);
      }

      // make sure every rank throws in case any rank fails
      int finished = 0;
      int index = 0;
      int flag = 0;
      while (finished < np && size > 1) {
        MPI_Testany(np, req.data(), &index, &flag, MPI_STATUS_IGNORE);
        if (flag) {
          finished++;
          if(is_invalid) continue;
          int id = index;
          int begin_id = EM::begin(id, events_rec) * nCmt;
          if (std::isnan(res(begin_id))) {
            is_invalid = true;
            rank_fail_msg << "Rank " << rank << " received invalid data for id " << id;
          }
        }
      }

      MPI_Barrier(comm);

      if(is_invalid) {
        MPI_Barrier(comm);
        throw std::runtime_error(rank_fail_msg.str());
      }
    }
#else

    /*
     * For population input in the form of ragged arrays,
     * addional information of the size of each individual
     * is required to locate the data in a single array for population.
     */
    template<typename T_events_record, typename... Ts> //NOLINT
    void pred(const T_events_record& events_rec,
                     Eigen::Matrix<typename EventsManager<T_events_record>::T_scalar, -1, -1>& res,
                     const T_pred... pred_pars,
                     const Ts... model_pars) {
      using ER = T_events_record;
      using EM = EventsManager<ER>;

      const int nCmt = EM::nCmt(events_rec);
      const int np = events_rec.num_subjects();
      
      res.resize(nCmt, events_rec.total_num_event_times);

      static bool has_warning = false;
      if (!has_warning) {
        std::cout << "Torsten Population PK solver " << "running sequentially" << "\n";
        has_warning = true;
      }

      for (int id = 0; id < np; ++id) {
        const int nKeep = events_rec.num_event_times(id);
        Eigen::Matrix<typename EventsManager<T_events_record>::T_scalar, -1, -1> res_id(nCmt, nKeep);
        pred(id, events_rec, res_id, pred_pars..., model_pars...);
        for (int j = 0; j < nKeep; ++j) {
          res.col(EM::begin(id, events_rec) + j) = res_id.col(j);
        }
      }
    }
#endif
  };

  template<typename T_model, typename... T_pred>
  constexpr double EventSolver<T_model, T_pred...>::invalid_res_d;
}
#endif
