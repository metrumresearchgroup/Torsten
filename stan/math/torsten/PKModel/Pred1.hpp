// version 0.8
// DEV - I lay the foundation for models that are based on a
// system of linear ODE expressed through a matrix, but this is 
// still work in progress 

#ifndef STAN_MATH_TORSTEN_PKMODEL_PRED1_HPP
#define STAN_MATH_TORSTEN_PKMODEL_PRED1_HPP

#include <iostream>
#include <Eigen/Dense>
#include <stan/math/prim/scal/err/check_finite.hpp>
#include <stan/math/torsten/PKModel/Pred/Pred1_oneCpt.hpp>
#include <stan/math/torsten/PKModel/Pred/Pred1_twoCpt.hpp>
#include <stan/math/torsten/PKModel/Pred/Pred1_general_solver.hpp>

using std::vector;
using namespace Eigen;
using std::string;

typedef Matrix<double, Dynamic, Dynamic> DMatrix;

/**
 *   The Functor of Pred1, which predicts amounts in each compartment 
 *   at one event. Defines a class of Pred1 functions (functors). The 
 *   key components is the constructors used to create the pred operator.
 *   Calls the different versions of Pred1 stored in the Pred directory.
 *
 *   Built-in Model types:
 *       1 - One Compartment Model with first-order absorption
 *       2 - Two Compartment Model with first-order absorption
 *		 3 - General Compartment Model using numerical ODE solver
 *		 4 - EXPERIMENTAL: PKPD model using semi-analytical solver 
 *		   
 *	 @tparam T_time type of scalar for time
 *	 @tparam T_rate type of scalar for rate
 *	 @tparam T_parameters type of scalar for model parameters
 *	 @tparam F type of ODE system function
 *	 @param[in] dt time between current and previous event
 *	 @param[in] parameter model parameters at current event
 *	 @param[in] init amount in each compartment at previous event
 *	 @param[in] rate rate in each compartment
 *	 @param[in] f functor for base ordinary differential equation that defines 
 *              compartment model.
 *   @return an eigen vector that contains predicted amount in each compartment 
 *           at the current event.      
 */
struct Pred1_structure {

private:
    string modeltype;
    DMatrix K;
    
public:
    
    Pred1_structure() {  //default constructor
        Matrix<double,1,1> init_Matrix;
        init_Matrix(0,0) = 0;
        modeltype = "default";
        K = init_Matrix; // Matrix = [0]
    }
    
    Pred1_structure(string p_modeltype) {  //constructor for standard models
        Matrix<double,1,1> init_Matrix;
        init_Matrix(0,0) = 0;
        
        modeltype = p_modeltype;
        K = init_Matrix; // Matrix = [0]
    }

	template<typename T_time, typename T_rate, typename T_parameters, typename F>
	Matrix<typename promote_args< T_time, T_rate, T_parameters>::type, 1, Dynamic> 
    operator()( const T_time& dt,
    		    const ModelParameters<T_time, T_parameters>& parameter,
    		    const Matrix<typename promote_args<T_time, T_rate,
    		      T_parameters>::type, 1, Dynamic>& init, 
    		    const vector<T_rate>& rate,
    		    const F& f) {
    	
    	stan::math::check_finite("Pred1", "initial values", init);

    	typedef typename promote_args< T_time, T_rate, T_parameters>::type scalar; 
    	
        if(modeltype == "OneCptModel")
          return Pred1_one(dt, parameter, init, rate, f);
        else if(modeltype == "TwoCptModel")
          return Pred1_two(dt, parameter, init, rate, f);
        else if(modeltype == "GeneralCptModel")
          return Pred1_general_solver(dt, parameter, init, rate, f); 
        else {
        	Matrix<scalar, 1, Dynamic> default_pred =
        	  Matrix<scalar, 1, Dynamic>::Zero(1);
         	return default_pred;
        }
    }
};

#endif
