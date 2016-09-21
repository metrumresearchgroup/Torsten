// Define a few simple functions. I guess I didn't deem them worthy of 
// having their own file, but this might not be the best approach
// FIX ME: deprecate the print functions. 
// FIX ME: each function should have its own header file
// FIX ME: using statements should be inside the scope of functions

#ifndef STAN_MATH_TORSTEN_PKMODEL_FUNCTIONS_HPP
#define STAN_MATH_TORSTEN_PKMODEL_FUNCTIONS_HPP

#include <vector>
#include <iostream>
#include <string>
#include <Eigen/Dense>
#include <stan/math/torsten/PKModel/SearchReal.hpp>

using std::vector;
using std::string;
using stan::math::var;
using Eigen::Dynamic;
using Eigen::Matrix;

/**
 * Define some simple and handy functions.
 *
 * print(): print input and skips a line if second argument is TRUE.
 * Function is overloaded for AutoDiff variables
 *
 * printpMatrix(): print the variable type of elements in Matrix,
 * and the elements contained inside the matrix.
 *
 * line(): print an empty line
 *
 * Marker(): prints 'marker', followed by the number of
 * the marker. Useful for tracking lines at which errors occur.
 *
 * find_time(): determines whether a value is present inside
 * a vector. Used to determine if times in the Events Schedule
 * are present in the ModelParameterHistory object.
 *
 * void_function(): the void function does nothing. It can be used
 * as a dummy function to call Pred for models that use analytical 
 * solutions and do not require a system of ODEs to be solved. 
 */

////////// FOR DEVELOPMENT PURPOSES ////////////// 

template <typename printtype>
inline void print(printtype T, bool skipline){	
	if(skipline==true) std::cout<<T<<std::endl;
	else std::cout<<T<<" ";
	}

inline void print(stan::math::var x, bool skipline){
	if(skipline==true) std::cout<< x.val() <<std::endl;
	else std::cout<<x.val()<<" ";
	}


void printpMatrix(vector <Matrix <var, Dynamic, 1> > pMatrix)
	{	
		print("PMatrix contains var.",true);
		int i, j;
		for(i=0;i<pMatrix.size();i++){
		for(j=0; j<pMatrix[0].rows(); j++){
			print(pMatrix[i](j,0).val(),false);
			}
			print(" ",true);
		}
	}

void printpMatrix(vector < Matrix <double, Dynamic, 1> > pMatrix)
	{	
		print("PMatrix contains double.",true);
		int i, j;
		for(i=0;i<pMatrix.size();i++){
		for(j=0; j<pMatrix[0].rows(); j++){
			print(pMatrix[i](j,1),false);
			}
			print(" ",true);
		}
	}

void line() {print(" ",true);}

int marker_count = 0; //define global variable
void Marker()
{
	print("MARKER",false);
	print(marker_count,true);
	marker_count++;
}


////////// REQUIRED FOR TORSTEN //////////////

template<typename T>
bool find_time(vector<T> v, T time)
{
	int k, size;
	bool found=false;
	size = v.size();
	k = SearchReal(v, size, time); //find the index of the largest
                                   //element in v smaller than time.
	if((k <= size) && (time == v[k-1])) found = true; 
	return found;
}


// construct dummy ODE 
template <typename T0, typename T1, typename T2, typename T3>
inline
std::vector<typename boost::math::tools::promote_args<T0, T1, T2, T3>::type>
f(const T0& t,
  const vector<T1>& x,
  const vector<T2>& parms,
  const vector<T3>& rate,
  const vector<int>& dummy,
  std::ostream* pstream__)
{
		typedef typename boost::math::tools::promote_args<T0, T1, T2, T3>::type scalar;
		vector<scalar> returned_vector = vector<scalar>(0,scalar(0));
		return returned_vector;
}		

struct dummy_ode
{
	template <typename T0, typename T1, typename T2, typename T3>
	inline
	vector<typename boost::math::tools::promote_args<T0, T1, T2, T3>::type>
	operator()(const T0& t,
			   const vector<T1>& x,
			   const vector<T2>& parms,
			   const vector<T3>& rate,
			   const vector<int>& dummy, 
			   std::ostream* pstream__) const
			   {
			   		return f(t, x, parms, rate, dummy, pstream__);
			   }
};

#endif