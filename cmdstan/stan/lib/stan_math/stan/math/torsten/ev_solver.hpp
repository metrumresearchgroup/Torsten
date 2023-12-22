#ifndef STAN_MATH_TORSTEN_EVENT_SOLVER_HPP
#define STAN_MATH_TORSTEN_EVENT_SOLVER_HPP

#include <stan/math/prim/fun/multiply.hpp>
#include <stan/math/rev/fun/multiply.hpp>
#include <stan/math/torsten/dsolve/pk_vars.hpp>
#include <stan/math/torsten/dsolve/pmx_ode_integrator.hpp>
#include <stan/math/torsten/ev_manager.hpp>
#include <stan/math/torsten/model_factory.hpp>
#include <stan/math/torsten/event.hpp>
#include <stan/math/torsten/mpi/session.hpp>
#include <stan/math/torsten/mpi/my_worker.hpp>
#include <stan/math/torsten/mpi/precomputed_gradients.hpp>
#include <Eigen/Dense>
#include <vector>

namespace torsten{

  template<typename T_model, typename T_em>
  struct EventSolver;

  /** 
   * The solver that solves events one by one according to @c evid and
   * dosing information. It is a wrapper that is aware of @c T_model so it build model
   * and select solver accordingly.
   * 
   * @tparam T_model type of model
   * @tparam T_event_record type of event record
   * @tparam T0 TIME type 
   * @tparam T4 THETA type
   * @theta_container container for theta type, could be
   *                  @c std::vector or @c eigen::vector.
   * 
   */
  template<typename T_model, typename T_event_record,
           typename T0, typename T4, template<typename...> class theta_container,
           typename... event_ctrl_type,
           typename... ode_data_type>
  struct EventSolver<T_model,
                     EventsManager<T_event_record, NonEventParameters<T0, T4, theta_container,
                                                                      std::tuple<event_ctrl_type...>,
                                                                      ode_data_type...>> > {
    using EM = EventsManager<T_event_record, NonEventParameters<T0, T4, theta_container,
                                                                std::tuple<event_ctrl_type...>,
                                                                ode_data_type...>>;
    using T_param = NonEventParameters<T0, T4, theta_container,
                                       std::tuple<event_ctrl_type...>,
                                       ode_data_type...>;
    /**
     * Data used to fill the results when computation throws exception.
     */
    static constexpr double invalid_res_d = std::numeric_limits<double>::quiet_NaN();

    /**
     * calculate the number of parameters according to the
     * model type.
     *
     * @tparam T_er type of events record
     * @param id subject id among a population
     * @param rec events record for the population
     */
    static int system_size(int id, const T_event_record& rec,
                           const std::vector<theta_container<T4>>& theta) {
      const int ncmt = rec.ncmt;
      const int npar = theta[0].size();
      const int nvar = pmx_model_nvars<T_model,
                                       typename EM::T_time,
                                       typename EM::T_scalar,
                                       typename EM::T_rate>::nvars(ncmt, npar);
      const int nvar_ss = pmx_model_nvars_ss<T_model,
                                             typename EM::T_amt,
                                             typename EM::T_par_rate,
                                             typename EM::T_par_ii>::nvars(npar);
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
     */
    template<typename integrator_type, typename... scalar_pars_type>
    void pred(int id,
              const T_event_record& events_rec,
              Eigen::Matrix<typename EM::T_scalar, -1, -1>& res,
              const integrator_type& integrator,
              const std::vector<theta_container<T4>>& theta,
              const std::vector<std::vector<event_ctrl_type>>&... event_ctrl,
              const std::vector<std::vector<ode_data_type>>&... ode_data,
              const scalar_pars_type... scalar_pars) {
      using Eigen::Matrix;
      using Eigen::Dynamic;
      using std::vector;
      using torsten::PKRec;

      int nCmt = EM::nCmt(events_rec);

      res.resize(nCmt, events_rec.num_event_times(id));
      PKRec<typename EM::T_scalar> init(nCmt);
      init.setZero();

      try {
        EM em(id, events_rec, theta, event_ctrl..., ode_data...);
        int ikeep = 0, iev = 0;

        while(ikeep < em.nKeep) {
          stepper(iev, init, em, integrator, scalar_pars...);
          if(em.events().keep(iev)) {
            res.col(ikeep) = init;
            ikeep++;
          }
          iev++;
        }
      } catch (const std::exception& e) {
        throw;
      }
    }

    /**
     * Step through a range of events.
     */
    template<typename integrator_type, typename... scalar_pars_type>
    void stepper(int i, PKRec<typename EM::T_scalar>& init, const EM& em,
                 const integrator_type& integrator,
                 const scalar_pars_type... scalar_pars) {
      auto& events = em.events();
      using scalar = typename EM::T_scalar;
      typename EM::T_time tprev = i == 0 ? events.time(0) : events.time(i-1);

      T_model pkmodel {model_factory<T_model, EM, scalar_pars_type...>::model(em, i,
                                                                              par_index_seq<ode_data_type...>{},
                                                                              scalar_pars...)};
      auto ev = em.event(i);
      ev(init, pkmodel, integrator);
    }

    /**
     * Solve an event from event manager and store the result in both
     * data format(solution & gradient) and <code>var</code> format.
     *
     * @param i id of the event to be solved
     * @param init <code>var</code> solution
     * @param sol_d data solution
     * @param em event manager
     * @param integrator numerical integrator for ODE
     * @param scalar_pars params needed to construst the PMX model,
     *        such as functor & dimension of ODE system
     */
    template<typename integrator_type, typename... scalar_pars_type>
    void stepper_solve(int i, torsten::PKRec<typename EM::T_scalar>& init,
                       torsten::PKRec<double>& sol_d,
                       const EM& em,
                       const integrator_type& integrator,
                       const scalar_pars_type... scalar_pars) {
      using std::vector;
      using stan::math::var;

      auto& events = em.events();

      typename EM::T_time tprev = i == 0 ? events.time(0) : events.time(i-1);

      T_model pkmodel {model_factory<T_model, EM, scalar_pars_type...>::model(em, i,
                                                                              par_index_seq<ode_data_type...>{},
                                                                              scalar_pars...)};
      auto ev = em.event(i);
      ev(sol_d, init, pkmodel, integrator, scalar_pars...);
    }

    template<typename integrator_type, typename... scalar_pars_type>
    /**
     * For MPI solutions, after each rank solves its corresponding
     * subjects, all ranks sync their results among each other.
     *
     * @param i id of the event to be synced
     * @param init <code>var</code> solution
     * @param sol_d data solution
     * @param em event manager
     * @param integrator numerical integrator for ODE
     * @param scalar_pars params needed to construst the PMX model,
     *        such as functor & dimension of ODE system
     */
    void stepper_sync(int i, torsten::PKRec<typename EM::T_scalar>& init,
                      torsten::PKRec<double>& sol_d,
                      const EM& em,
                      const integrator_type& integrator,
                      const scalar_pars_type... scalar_pars) {
      using std::vector;
      using stan::math::var;

      auto& events = em.events();

      typename EM::T_time tprev = i == 0 ? events.time(0) : events.time(i-1);

      if (events.is_reset(i)) {
        init.setZero();
      } else if (events.is_ss_dosing(i)) {  // steady state event
        T_model pkmodel {model_factory<T_model, EM, scalar_pars_type...>::model(em, i,
                                                                                par_index_seq<ode_data_type...>{},
                                                                                scalar_pars...)};
        auto curr_amt = em.fractioned_amt(i);
        vector<var> v_i(dsolve::pk_vars(curr_amt, events.rate(i), events.ii(i), pkmodel.par()));
        int nsys = torsten::pk_nsys(em.ncmt, v_i.size());
        if (events.ss(i) == 2)
          init += torsten::mpi::precomputed_gradients(sol_d.segment(0, nsys), v_i);  // steady state without reset
        else
          init = torsten::mpi::precomputed_gradients(sol_d.segment(0, nsys), v_i);  // steady state with reset (ss = 1)
      } else if (events.time(i) > tprev) {
        // auto curr_rates = stan::math::value_of(em.fractioned_rates(i));
        auto curr_rates = em.fractioned_rates(i);
          T_model pkmodel {model_factory<T_model, EM, scalar_pars_type...>::model(em, i,
                                                                                  par_index_seq<ode_data_type...>{},
                                                                                  scalar_pars...)};
          vector<var> v_i =
            pmx_model_vars<T_model>::vars(events.time(i), init, curr_rates, pkmodel.par());
          int nsys = torsten::pk_nsys(em.ncmt, v_i.size());
          init = torsten::mpi::precomputed_gradients(sol_d.segment(0, nsys), v_i);
      }

      if (events.is_bolus_dosing(i)) {
        init(events.cmt(i) - 1) += em.fractioned_amt(i);
      }
    }

#ifdef TORSTEN_MPI

    /**
     * MPI solution when the population
     * information passed in as ragged arrays. The overloading occurs
     * on <code>res</code> arg for results. Here the result is <code>var</code>.
     *
     * @param events_rec event record
     * @param res solution 
     * @param integrator ODE integrator
     * @param theta PMX parameters passed into ODE function
     * @param event_ctrl optional PMX parameters: bioavailability & tlag
     * @param ode_data optional ODE parameters: real & integer data
     * @param scalar_pars params needed to construst the PMX model,
     *        such as functor & dimension of ODE system
     */
    template<typename integrator_type, typename... scalar_pars_type>
    void pred(const T_event_record& events_rec,
              Eigen::Matrix<stan::math::var, -1, -1>& res,
              const integrator_type& integrator,
              const std::vector<theta_container<T4>>& theta,
              const std::vector<std::vector<event_ctrl_type>>&... event_ctrl,
              const std::vector<std::vector<ode_data_type>>&... ode_data,
              const scalar_pars_type... scalar_pars) {
      using Eigen::Matrix;
      using Eigen::MatrixXd;
      using Eigen::VectorXd;
      using Eigen::Dynamic;
      using std::vector;
      using::stan::math::var;
      using torsten::PKRec;

      using scalar = typename EM::T_scalar;

      const int nCmt = EM::nCmt(events_rec);
      const int np = events_rec.num_subjects();
      bool is_invalid = false;
      std::ostringstream rank_fail_msg;

      const MPI_Comm& comm = torsten::mpi::Session::pmx_parm_comm().comm();
      int rank = torsten::mpi::Session::pmx_parm_comm().rank();
      int size = torsten::mpi::Session::pmx_parm_comm().size();

      std::vector<MPI_Request> req(np);
      vector<MatrixXd> res_d(np);
      
      res.resize(nCmt, events_rec.total_num_event_times);

      PKRec<scalar> init(nCmt);
      PKRec<double> pred1;
      for (int id = 0; id < np; ++id) {

        /* For every rank */

        const int nKeep = events_rec.num_event_times(id);

        int nev = EM::num_events(id, events_rec, theta, event_ctrl..., ode_data...);
        res_d[id].resize(system_size(id, events_rec, theta), nev);
        res_d[id].setConstant(0.0);

        int my_worker_id = torsten::mpi::my_worker(id, np, size);

        /* only solver rank */

        if (rank == my_worker_id) {
          if (is_invalid) {
            res_d[id].setConstant(invalid_res_d);
          } else {
            try {
              EM em(id, events_rec, theta, event_ctrl..., ode_data...);
              auto& events = em.events();
              assert(nev == events.size());
              assert(nKeep == em.nKeep);
              init.setZero();

              int ikeep = 0, iev = 0;
              while(ikeep < em.nKeep) {
                stepper_solve(iev, init, pred1, em, integrator, scalar_pars...);
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
            // assert(rank != torsten::mpi::my_worker(id, np, size));
            is_invalid = true;
            rank_fail_msg << "Rank " << rank << " received invalid data for id " << id;
          } else {
            EM em(id, events_rec, theta, event_ctrl..., ode_data...);
            auto& events = em.events();
            PKRec<scalar> init(nCmt); init.setZero();
            PKRec<double> pred1 = VectorXd::Zero(res_d[id].rows());

            int ikeep = 0, iev = 0;
            while(ikeep < em.nKeep) {
              pred1 = res_d[id].col(iev);
              stepper_sync(iev, init, pred1, em, integrator, scalar_pars...);
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
    template<typename integrator_type, typename... scalar_pars_type>
    void pred(const T_event_record& events_rec, Eigen::MatrixXd& res,
              const integrator_type& integrator,
              const std::vector<theta_container<T4>>& theta,
              const std::vector<std::vector<event_ctrl_type>>&... event_ctrl,
              const std::vector<std::vector<ode_data_type>>&... ode_data,
              const scalar_pars_type... scalar_pars) {
      using Eigen::Matrix;
      using Eigen::MatrixXd;
      using Eigen::VectorXd;
      using Eigen::Dynamic;
      using std::vector;
      using::stan::math::var;
      using torsten::PKRec;

      using ER = NONMENEventsRecord<double, double, double, double>;
      using EM_d = EventsManager<ER, NonEventParameters<T0, T4, theta_container, std::tuple<event_ctrl_type...>, ode_data_type...>>;

      const int nCmt = EM_d::nCmt(events_rec);
      const int np = events_rec.num_subjects();
      bool is_invalid = false;
      std::ostringstream rank_fail_msg;

      const MPI_Comm& comm = torsten::mpi::Session::pmx_data_comm().comm();
      int rank = torsten::mpi::Session::pmx_data_comm().rank();
      int size = torsten::mpi::Session::pmx_data_comm().size();

      std::vector<MPI_Request> req(np);

      res.resize(nCmt, events_rec.total_num_event_times);

      PKRec<double> init(nCmt);
      for (int id = 0; id < np; ++id) {

        /* For every rank */
        int nKeep = events_rec.num_event_times(id);
        int my_worker_id = torsten::mpi::my_worker(id, np, size);
        int begin_id = EM_d::begin(id, events_rec) * nCmt;
        int size_id = nKeep * nCmt;

        /* only solver rank */
        if (rank == my_worker_id) {
          try {
            EM em(id, events_rec, theta, event_ctrl..., ode_data...);
            auto& events = em.events();
            init.setZero();
            int ikeep = 0, iev = 0;
            while(ikeep < em.nKeep) {
              stepper(iev, init, em, integrator, scalar_pars...);
              if(events.keep(iev)) {
                res.col(EM_d::begin(id, events_rec) + ikeep) = init;
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
          int begin_id = EM_d::begin(id, events_rec) * nCmt;
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
    template<typename integrator_type, typename... scalar_pars_type> //NOLINT
    void pred(const T_event_record& events_rec,
              Eigen::Matrix<typename EM::T_scalar, -1, -1>& res,
              const integrator_type& integrator,
              const std::vector<theta_container<T4>>& theta,
              const std::vector<std::vector<event_ctrl_type>>&... event_ctrl,
              const std::vector<std::vector<ode_data_type>>&... ode_data,
              const scalar_pars_type... scalar_pars) {
      using ER = T_event_record;

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
        Eigen::Matrix<typename EM::T_scalar, -1, -1> res_id(nCmt, nKeep);
        pred(id, events_rec, res_id, integrator, theta,
             event_ctrl..., ode_data..., scalar_pars...);
        for (int j = 0; j < nKeep; ++j) {
          res.col(EM::begin(id, events_rec) + j) = res_id.col(j);
        }
      }
    }
#endif
  };

  template<typename T_model, typename T_event_record,
           typename T0, typename T4, template<typename...> class theta_container,
           typename... event_ctrl_type,
           typename... ode_data_type>
  constexpr double EventSolver<T_model,
                               EventsManager<T_event_record, NonEventParameters<T0, T4, theta_container,
                                                                                std::tuple<event_ctrl_type...>,
                                                                                ode_data_type...>> >::invalid_res_d;
}
#endif
