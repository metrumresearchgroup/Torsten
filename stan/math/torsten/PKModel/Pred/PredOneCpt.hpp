// version 0.8
// EXPERIMENTAL: will likely not be included in final library 


#ifndef PKMODEL_PRED_PREDONECPT_HPP
#define PKMODEL_PRED_PREDONECPT_HPP

#include <iostream>
#include "PolyExp.hpp"

using std::vector;
using boost::math::tools::promote_args;
using Eigen::Matrix;
using Eigen::Dynamic;
typedef vector<double> dvector;
typedef vector<int> ivector;

template<typename T_time, typename T_rate, typename T_parameters, typename F>
Matrix<typename promote_args< T_time, T_rate, T_parameters>::type, 1, Dynamic> 
Pred1_one(const T_time& dt,
		      const vector<T_parameters>& parameters, 
		      const vector<typename promote_args<T_time, T_rate, T_parameters>::type>& init,
		      const vector<T_rate>& rate)
{
	typedef typename promote_args<T_time, T_rate, T_parameters>::type scalar;
	
	T_parameters CL, V2, ka, k10;
	vector<scalar> a(2,0); 
	vector<T_parameters> alpha(2,0); 
	Matrix<scalar, 1, Dynamic> pred = Matrix<scalar, 1, Dynamic>::Zero(2); //initialize pred to a row-vector of 
											                                                  	//length 2, with elements equal to 0. 
	CL = parameters[0];
	V2 = parameters[1];
	ka = parameters[2];
	
	assert(((CL>0)&&(V2>0))&&(ka>0));
	k10 = CL/V2;
	alpha[0] = k10;
	alpha[1] = ka;
	
	if((init[0]!=0)||(rate[0]!=0))
	{
	
		pred(0,0)=init[0]*exp(-ka*dt) + rate[0]*(1-exp(-ka*dt))/ka;
		a[0] = ka/(ka-alpha[0]);
		a[1] = -a[0];
		pred(0,1) += PolyExp(dt, init[0],0,0,0,false,a,alpha,2) +
				   PolyExp(dt,0,rate[0],dt,0,false,a,alpha,2);
	}
			
	if((init[1]!=0)||(rate[1]!=0))
	{
		a[0]=1;
		pred(0,1) += PolyExp(dt,init[1],0,0,0,false,a,alpha,1)
				   + PolyExp(dt,0,rate[1],dt,0,false,a,alpha,1);
		
	}	
		
	return pred;
}
		
#endif