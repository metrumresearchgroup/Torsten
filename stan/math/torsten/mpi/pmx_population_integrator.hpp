#ifndef STAN_MATH_TORSTEN_DSOLVE_PMX_POPULATION_INTEGRATOR_HPP
#define STAN_MATH_TORSTEN_DSOLVE_PMX_POPULATION_INTEGRATOR_HPP

#include <stan/math/prim/mat/fun/to_matrix.hpp>
#include <stan/math/torsten/mpi/precomputed_gradients.hpp>
#include <stan/math/torsten/mpi/session.hpp>
#include <stan/math/torsten/mpi/my_worker.hpp>
#include <stan/math/torsten/return_type.hpp>
#include <stan/math/torsten/is_var.hpp>
#include <algorithm>
#include <vector>

namespace torsten {
  namespace mpi {

    /*
     * For MPI integrator the workflow is to solve at
     * designated rank and pass the results to the rest, while
     * waiting for the results from the other ranks to be
     * passed over. The results are in the form of data so
     * they can be sent over through MPI before assembled
     * into @c var at each rank. The data-only scenario is treated
     * separately as it doesn't require @c var generation.
     * Note that when one rank fails to solve the ODE and
     * throws an exception we need make sure the rest ranks
     * also throw to avoid locking, for that purpose an
     * invalid value is filled into the data passed to the rest
     * ranks so they can detect the exception.
     */
    // template <typename F, typename ode_t, typename service_t>
    template <typename F, typename solver_t, template<typename...> class ode_t, typename... ode_pars_t>
    struct PMXPopulationIntegrator {

      solver_t& solver;
      static constexpr double invalid_res_d = std::numeric_limits<double>::quiet_NaN();

      PMXPopulationIntegrator(solver_t& solver0) : solver(solver0) {}

#ifdef TORSTEN_MPI
      /*
       * MPI solution when the parameters contain @c var,
       * assuming per-subject parameters.
       *
       * @tparam Tt time type
       * @tparam T_initial @c y0 type
       * @tparam T_param @c theta type
       * @param f functor for the RHS of ODE
       * @param y0 init condition
       * @param t0 starting time
       * @param len length of data for each subject in @c ts
       * @param ts time step for the entire group, concatenated.
       * @param theta parameter for each ODE among the group
       * @param x_r real data for each ODE among the group
       * @param x_i int data for each ODE among the group
       */ 
      template <typename Tt, typename T_initial, typename T_param>
      inline
      Eigen::Matrix<typename torsten::return_t<Tt, T_initial, T_param>::type, // NOLINT
                    Eigen::Dynamic, Eigen::Dynamic>
      operator()(const F& f,
                 const std::vector<std::vector<T_initial> >& y0,
                 double t0,
                 const std::vector<int>& len,
                 const std::vector<Tt>& ts,
                 const std::vector<std::vector<T_param> >& theta,
                 const std::vector<std::vector<double> >& x_r,
                 const std::vector<std::vector<int> >& x_i,
                 std::ostream* msgs) {
        using std::vector;
        using Eigen::Matrix;
        using Eigen::MatrixXd;
        using Eigen::Dynamic;
        const int m = theta[0].size();
        const int n = y0[0].size();
        const int np = theta.size(); // population size

        using Ode = ode_t<F, Tt, T_initial, T_param, ode_pars_t...>;
        torsten::dsolve::PMXOdeService<Ode> serv(n, m);
    
        MPI_Comm comm = torsten::mpi::Session<NUM_TORSTEN_COMM>::comms[TORSTEN_COMM_ODE_PARM].comm;
        int rank = torsten::mpi::Session<NUM_TORSTEN_COMM>::comms[TORSTEN_COMM_ODE_PARM].rank;
        int size = torsten::mpi::Session<NUM_TORSTEN_COMM>::comms[TORSTEN_COMM_ODE_PARM].size;

        using scalar_type = typename torsten::return_t<Tt, T_initial, T_param>::type;

        vector<MatrixXd> res_d(np);
        Matrix<scalar_type, -1, -1> res(n, ts.size());
        res.setZero();
        std::vector<double> g;
        std::vector<MPI_Request> req(np);

        bool is_invalid = false;
        std::ostringstream rank_fail_msg;

        // n_req will be updated to contain the number of
        // active requests. These are the requests that we
        // need to make sure to finish.
        int n_req = np;

        typename std::vector<Tt>::const_iterator iter = ts.begin();
        for (int i = 0; i < np; ++i) {
          std::vector<Tt> ts_i(iter, iter + len[i]);
          iter += len[i];
          int my_worker_id = torsten::mpi::my_worker(i, np, size);
          try {
            Ode ode{serv, f, t0, ts_i, y0[i], theta[i], x_r[i], x_i[i], msgs};
            res_d[i].resize(ode.fwd_system_size(), ode.ts().size());
            res_d[i].setConstant(0.0);

            // success in creating ODE, solve it
            if(rank == my_worker_id) {
              try {
                res_d[i] = solver. template integrate<Ode, false>(ode);
              } catch (const std::exception& e) {
                is_invalid = true;
                res_d[i].setConstant(invalid_res_d);
                rank_fail_msg << "Rank " << rank << " failed to solve ODEs for id " << i << ": " << e.what();
              }
            }
          } catch (const std::exception& e) {
            is_invalid = true;
            rank_fail_msg << "Rank " << rank << " failed to create ODEs for id " << i << ": " << e.what();
            n_req = i;
            break;
          }

          MPI_Ibcast(res_d[i].data(), res_d[i].size(), MPI_DOUBLE, my_worker_id, comm, &req[i]);
        }

        int begin_id = 0;
        iter = ts.begin();
        for (int i = 0; i < n_req; ++i) {
          std::vector<Tt> ts_i(iter, iter + len[i]);
          iter += len[i];
          MPI_Wait(&req[i], MPI_STATUS_IGNORE);
          if(is_invalid) continue;
          if (std::isnan(res_d[i](0))) {
            is_invalid = true;
            rank_fail_msg << "Rank " << rank << " received invalid data for id " << i;
          } else {
            Ode ode{serv, f, t0, ts_i, y0[i], theta[i], x_r[i], x_i[i], msgs};
            Matrix<scalar_type, -1, -1>::Map(&res(begin_id), n, len[i]) = torsten::precomputed_gradients(res_d[i], ode.vars());
          }
          begin_id += len[i] * n;
        }
        
        if(is_invalid) {
          throw std::runtime_error(rank_fail_msg.str());
        }

        return res;
      }

      /*
       * MPI solution when the parameters contain @c var,
       * assuming seperate group-wide and per-subject parameters.
       *
       * @tparam Tt time type
       * @tparam T_initial @c y0 type
       * @tparam T_param @c theta type
       * @param f functor for the RHS of ODE
       * @param y0 init condition
       * @param t0 starting time
       * @param len length of data for each subject in @c ts
       * @param ts time step for the entire group, concatenated.
       * @param group_theta ODE parameters shared among the group
       * @param theta parameter for each ODE among the group
       * @param x_r real data for each ODE among the group
       * @param x_i int data for each ODE among the group
       */ 
      template <typename Tt, typename T_initial, typename T_param>
      inline
      Eigen::Matrix<typename torsten::return_t<Tt, T_initial, T_param>::type, // NOLINT
                    Eigen::Dynamic, Eigen::Dynamic>
      operator()(const F& f,
                 const std::vector<std::vector<T_initial> >& y0,
                 double t0,
                 const std::vector<int>& len,
                 const std::vector<Tt>& ts,
                 const std::vector<T_param>& group_theta,
                 const std::vector<std::vector<T_param> >& theta,
                 const std::vector<std::vector<double> >& x_r,
                 const std::vector<std::vector<int> >& x_i,
                 std::ostream* msgs) {
        using std::vector;
        using Eigen::Matrix;
        using Eigen::MatrixXd;
        using Eigen::Dynamic;
        using Ode = ode_t<F, Tt, T_initial, T_param, ode_pars_t...>;
        const int m = theta[0].size() + group_theta.size();
        const int n = y0[0].size();
        const int np = len.size(); // population size

        torsten::dsolve::PMXOdeService<Ode> serv(n, m);
    
        MPI_Comm comm = torsten::mpi::Session<NUM_TORSTEN_COMM>::comms[TORSTEN_COMM_ODE_PARM].comm;
        int rank = torsten::mpi::Session<NUM_TORSTEN_COMM>::comms[TORSTEN_COMM_ODE_PARM].rank;
        int size = torsten::mpi::Session<NUM_TORSTEN_COMM>::comms[TORSTEN_COMM_ODE_PARM].size;

        using scalar_type = typename torsten::return_t<Tt, T_initial, T_param>::type;

        vector<MatrixXd> res_d(np);
        Matrix<scalar_type, -1, -1> res(n, ts.size());
        res.setZero();
        std::vector<double> g;
        std::vector<MPI_Request> req(np);

        bool is_invalid = false;
        std::ostringstream rank_fail_msg;

        // n_req will be updated to contain the number of
        // active requests. These are the requests that we
        // need to make sure to finish.
        int n_req = np;

        vector<T_param> theta_i(m);
        std::copy(group_theta.begin(), group_theta.end(), theta_i.begin());

        typename std::vector<Tt>::const_iterator iter = ts.begin();
        for (int i = 0; i < np; ++i) {
          std::vector<Tt> ts_i(iter, iter + len[i]);
          iter += len[i];
          int my_worker_id = torsten::mpi::my_worker(i, np, size);
          try {
            std::copy(theta[i].begin(), theta[i].end(), theta_i.begin() + group_theta.size());
            Ode ode{serv, f, t0, ts_i, y0[i], theta_i, x_r[i], x_i[i], msgs};
            res_d[i].resize(ode.fwd_system_size(), ode.ts().size());
            res_d[i].setConstant(0.0);

            // success in creating ODE, solve it
            if(rank == my_worker_id) {
              try {
                res_d[i] = solver. template integrate<Ode, false>(ode);
              } catch (const std::exception& e) {
                is_invalid = true;
                res_d[i].setConstant(invalid_res_d);
                rank_fail_msg << "Rank " << rank << " failed to solve ODEs for id " << i << ": " << e.what();
              }
            }
          } catch (const std::exception& e) {
            is_invalid = true;
            rank_fail_msg << "Rank " << rank << " failed to create ODEs for id " << i << ": " << e.what();
            n_req = i;
            break;
          }

          MPI_Ibcast(res_d[i].data(), res_d[i].size(), MPI_DOUBLE, my_worker_id, comm, &req[i]);
        }

        int begin_id = 0;
        iter = ts.begin();
        for (int i = 0; i < n_req; ++i) {
          std::vector<Tt> ts_i(iter, iter + len[i]);
          iter += len[i];
          MPI_Wait(&req[i], MPI_STATUS_IGNORE);
          if(is_invalid) continue;
          if (std::isnan(res_d[i](0))) {
            is_invalid = true;
            rank_fail_msg << "Rank " << rank << " received invalid data for id " << i;
          } else {
            std::copy(theta[i].begin(), theta[i].end(), theta_i.begin() + group_theta.size());
            Ode ode{serv, f, t0, ts_i, y0[i], theta_i, x_r[i], x_i[i], msgs};
            Matrix<scalar_type, -1, -1>::Map(&res(begin_id), n, len[i]) = torsten::precomputed_gradients(res_d[i], ode.vars());
          }
          begin_id += len[i] * n;
        }
        
        if(is_invalid) {
          throw std::runtime_error(rank_fail_msg.str());
        }

        return res;
      }

      /*
       * Data-only MPI solution when paramters are
       * per-subject based.
       *
       * @param f functor for the RHS of ODE
       * @param y0 init condition
       * @param t0 starting time
       * @param len length of data for each subject in @c ts
       * @param ts time step for the entire group, concatenated.
       * @param theta parameter for each ODE among the group
       * @param x_r real data for each ODE among the group
       * @param x_i int data for each ODE among the group
       */ 
      inline Eigen::MatrixXd operator()(const F& f,
                                        const std::vector<std::vector<double> >& y0,    // NOLINT
                                        double t0,
                                        const std::vector<int>& len,
                                        const std::vector<double>& ts,
                                        const std::vector<std::vector<double> >& theta, // NOLINT
                                        const std::vector<std::vector<double> >& x_r,   // NOLINT
                                        const std::vector<std::vector<int> >& x_i,      // NOLINT
                                        std::ostream* msgs) {
        using stan::math::var;
        using std::vector;
        using Eigen::Matrix;
        using Eigen::MatrixXd;
        using Eigen::Dynamic;
        using Ode = ode_t<F, double, double, double, ode_pars_t...>;
        const int m = theta[0].size();
        const int n = y0[0].size();
        const int np = theta.size(); // population size

        torsten::dsolve::PMXOdeService<Ode> serv(n, m);
    
        MPI_Comm comm = torsten::mpi::Session<NUM_TORSTEN_COMM>::comms[TORSTEN_COMM_ODE_DATA].comm;
        int rank = torsten::mpi::Session<NUM_TORSTEN_COMM>::comms[TORSTEN_COMM_ODE_DATA].rank;
        int size = torsten::mpi::Session<NUM_TORSTEN_COMM>::comms[TORSTEN_COMM_ODE_DATA].size;

        Eigen::MatrixXd res = Eigen::MatrixXd::Zero(n, ts.size());
        std::vector<MPI_Request> req(np);
        std::vector<int> begin_ids(np);

        bool is_invalid = false;
        int n_req = np;
        std::ostringstream rank_fail_msg;

        std::vector<double>::const_iterator iter = ts.begin();
        int begin_id = 0;
        for (int i = 0; i < np; ++i) {
          std::vector<double> ts_i(iter, iter + len[i]);
          iter += len[i];
          int size_id = len[i] * n;
          int my_worker_id = torsten::mpi::my_worker(i, np, size);
          try {
            Ode ode{serv, f, t0, ts_i, y0[i], theta[i], x_r[i], x_i[i], msgs};

            // success in creating ODE, solve it
            if(rank == my_worker_id) {
              try {
                MatrixXd::Map(&res(begin_id), n, len[i]) = solver.template integrate<Ode, false>(ode);
              } catch (const std::exception& e) {
                is_invalid = true;
                MatrixXd::Map(&res(begin_id), n, len[i]).setConstant(invalid_res_d);
                rank_fail_msg << "Rank " << rank << " failed to solve ODEs for id " << i << ": " << e.what();
              }
            }
          } catch (const std::exception& e) {
            is_invalid = true;
            rank_fail_msg << "Rank " << rank << " failed to create ODEs for id " << i << ": " << e.what();
            n_req = i;
            break;
          }

          MPI_Ibcast(&res(begin_id), size_id, MPI_DOUBLE, my_worker_id, comm, &req[i]);
          begin_ids[i] = begin_id;
          begin_id += size_id;
        }

        int finished = 0;
        int index;
        while(finished != n_req) {
          MPI_Waitany(n_req, req.data(), &index, MPI_STATUS_IGNORE);
          finished++;
          if(is_invalid) continue;
          int i = index;
          if (std::isnan(res(begin_ids[i]))) {
            is_invalid = true;
            rank_fail_msg << "Rank " << rank << " received invalid data for id " << i;
          }
        }

        if(is_invalid) {
          throw std::runtime_error(rank_fail_msg.str());
        }

        return res;
      }

      /*
       * Data-only MPI solution when paramters contain
       * group and per-subject components.
       *
       * @param f functor for the RHS of ODE
       * @param y0 init condition
       * @param t0 starting time
       * @param len length of data for each subject in @c ts
       * @param ts time step for the entire group, concatenated.
       * @param group_theta parameters shared among the group
       * @param theta parameter for each ODE among the group
       * @param x_r real data for each ODE among the group
       * @param x_i int data for each ODE among the group
       */ 
      inline Eigen::MatrixXd operator()(const F& f,
                                        const std::vector<std::vector<double> >& y0,    // NOLINT
                                        double t0,
                                        const std::vector<int>& len,
                                        const std::vector<double>& ts,
                                        const std::vector<double>& group_theta,
                                        const std::vector<std::vector<double> >& theta, // NOLINT
                                        const std::vector<std::vector<double> >& x_r,   // NOLINT
                                        const std::vector<std::vector<int> >& x_i,      // NOLINT
                                        std::ostream* msgs) {
        using std::vector;
        vector<vector<double> > theta_(theta.size(), vector<double>(group_theta.size() + theta[0].size()) ); // NOLINT
        std::transform(theta.begin(), theta.end(), theta_.begin(),
                       [&](const vector<double>& v1) {
                         vector<double> v(group_theta.size() + v1.size());
                         std::copy(group_theta.begin(), group_theta.end(), v.begin());
                         std::copy(v1.begin(), v1.end(), v.begin() + group_theta.size());
                         return v;
                       });
        return (*this)(f, y0, t0, len, ts, theta_, x_r, x_i, msgs);
      }

#else
      /*
       * In sequential run we simply solve the ODEs one-by-one.
       */
      template <typename Tt, typename T_initial, typename T_param>
      inline Eigen::Matrix<typename torsten::return_t<Tt, T_initial, T_param>::type, // NOLINT
                           Eigen::Dynamic, Eigen::Dynamic>
      operator()(const F& f,
                 const std::vector<std::vector<T_initial> >& y0,
                 double t0,
                 const std::vector<int>& len,
                 const std::vector<Tt>& ts,
                 const std::vector<std::vector<T_param> >& theta,
                 const std::vector<std::vector<double> >& x_r,
                 const std::vector<std::vector<int> >& x_i,
                 std::ostream* msgs) {
        using std::vector;
        using Eigen::Matrix;
        using Eigen::MatrixXd;
        using Eigen::Dynamic;
        using Ode = ode_t<F, Tt, T_initial, T_param, ode_pars_t...>;
        const int m = theta[0].size();
        const int n = y0[0].size();
        const int np = theta.size(); // population size

        static bool has_warning = false;
        if (!has_warning) {
          std::cout << "Warning: Torsten group ODE integrator " << "running sequentially" << "\n";
          has_warning = true;
        }

        static torsten::dsolve::PMXOdeService<Ode> serv(n, m);

        using scalar_type = typename torsten::return_t<Tt, T_initial, T_param>::type;
        Matrix<scalar_type, Dynamic, Dynamic> res(n, ts.size());

        typename std::vector<Tt>::const_iterator iter = ts.begin();
        int begin_id = 0;
        for (int i = 0; i < np; ++i) {
          std::vector<Tt> ts_i(iter, iter + len[i]);
          iter += len[i];
          Ode ode{serv, f, t0, ts_i, y0[i], theta[i], x_r[i], x_i[i], msgs};
          Matrix<scalar_type, -1, -1>::Map(&res(begin_id), n, len[i])
            = stan::math::to_matrix(solver.template integrate(ode)).transpose();
          begin_id += n * len[i];
        }

        assert(begin_id == res.size());

        return res;
      }

      /*
       * In sequential run we simply solve the ODEs one-by-one.
       */
      template <typename Tt, typename T_initial, typename T_param>
      inline Eigen::Matrix<typename torsten::return_t<Tt, T_initial, T_param>::type, // NOLINT
                           Eigen::Dynamic, Eigen::Dynamic>
      operator()(const F& f,
                 const std::vector<std::vector<T_initial> >& y0,
                 double t0,
                 const std::vector<int>& len,
                 const std::vector<Tt>& ts,
                 const std::vector<T_param>& group_theta,
                 const std::vector<std::vector<T_param> >& theta,
                 const std::vector<std::vector<double> >& x_r,
                 const std::vector<std::vector<int> >& x_i,
                 std::ostream* msgs) {
        using std::vector;
        using Eigen::Matrix;
        using Eigen::MatrixXd;
        using Eigen::Dynamic;
        using Ode = ode_t<F, Tt, T_initial, T_param, ode_pars_t...>;
        const int m = theta[0].size() + group_theta.size();
        const int n = y0[0].size();
        const int np = theta.size(); // population size

        static bool has_warning = false;
        if (!has_warning) {
          std::cout << "Torsten Population ODE solver " << "running sequentially" << "\n";
          has_warning = true;
        }

        static torsten::dsolve::PMXOdeService<Ode> serv(n, m);

        using scalar_type = typename torsten::return_t<Tt, T_initial, T_param>::type;
        Matrix<scalar_type, Dynamic, Dynamic> res(n, ts.size());

        vector<T_param> theta_i(m);
        std::copy(group_theta.begin(), group_theta.end(), theta_i.begin());

        vector<vector<scalar_type> > res_i;
        typename std::vector<Tt>::const_iterator iter = ts.begin();
        int begin_id = 0;
        for (int i = 0; i < np; ++i) {
          std::copy(theta[i].begin(), theta[i].end(), theta_i.begin() + group_theta.size());
          std::vector<Tt> ts_i(iter, iter + len[i]);
          iter += len[i];
          Ode ode{serv, f, t0, ts_i, y0[i], theta_i, x_r[i], x_i[i], msgs};
          Matrix<scalar_type, -1, -1>::Map(&res(begin_id), n, len[i])
            = stan::math::to_matrix(solver.template integrate(ode)).transpose();
          begin_id += n * len[i];
        }

        return res;
      }
#endif
    };

    template <typename F, typename solver_t, template<typename...> class ode_t, typename... ode_pars_t>
    constexpr double PMXPopulationIntegrator<F, solver_t, ode_t, ode_pars_t...>::invalid_res_d;

  }
}

#endif
