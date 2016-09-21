#ifndef STAN_MATH_PMETRICS_PKMODEL_PMETRICS_SOLVER_HPP
#define STAN_MATH_PMETRICS_PKMODEL_PMETRICS_SOLVER_HPP

#include <iostream>
#include <Eigen/Dense>
#include <stan/math/rev/mat/functor/integrate_ode_bdf.hpp>
#include <stan/math/prim/arr/functor/integrate_ode_rk45.hpp>

// FIX ME: using statements should appear within
// the scope of functions.
using std::vector;
using namespace Eigen;
using std::string;

/**
 *  Construct functors that run the ODE integrator. Specify integrator 
 *  type, the base ODE system, and the tuning parameters (relative tolerance,
 *  absolute tolerance, and maximum number of steps). 
 *	 
 */
struct pmetrics_solver_structure {

private:

	double rel_tol, abs_tol;
	long int max_num_steps;
	std::string solver_type; 
     
public:
        
    pmetrics_solver_structure() {
        rel_tol = 1e-10;
        abs_tol = 1e-10;
        max_num_steps = 1e8;
        solver_type = "bdf";
    }
    
    pmetrics_solver_structure (
      double p_rel_tol, double p_abs_tol, long int p_max_num_steps, 
      std::string p_solver_type) {
    	rel_tol = p_rel_tol;
    	abs_tol = p_abs_tol;
    	max_num_steps = p_max_num_steps;
    	solver_type = p_solver_type;
    }
    
	// CONSTRUCTOR FOR OPERATOR
	template<typename F, typename T1, typename T2>
	vector<vector<typename stan::return_type<T1, T2>::type> > // EXPLAIN ME! What is stan::return_type<T1, T2> ?
	operator()(const F& f,
			   const vector<T1> y0,
			   const double t0,
			   const vector<double>& ts,
			   const vector<T2>& theta,
			   const vector<double>& x,
			   const vector<int>& x_int) {	
		  
		  if(solver_type == "bdf") {
		    return stan::math::integrate_ode_bdf(f, y0, t0, ts, theta, x, x_int, 
                                             0, // std::ostream * msgs
								             rel_tol, abs_tol, max_num_steps);
		  }
		  else {  // if(solver_type == "rk45")
		    return stan::math::integrate_ode_rk45(f, y0, t0, ts, theta, x, x_int, 
                                              0, // std::ostream * msgs
                                              rel_tol, abs_tol, max_num_steps);
		  }
	}
};

#endif
