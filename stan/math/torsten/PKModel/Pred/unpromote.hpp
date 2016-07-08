// version 0.8

#ifndef PKMODEL_CONVERT_HPP
#define PKMODEL_CONVERT_HPP

#include <iostream>
using std::vector;
using stan::math::var;

// DEV: may need to overload functions to work on fvar.
//      also probably doing something a little bit wrong if I am using
//      this function. 


/**
 * Functions that converts an autodiff variable into a double. 
 * The variable will either be a stan::math::var
 * or a double (in which case it will not be modified).
 *
 */

double unpromote(var x) {return x.val();}
double unpromote(double x) {return x;}

#endif