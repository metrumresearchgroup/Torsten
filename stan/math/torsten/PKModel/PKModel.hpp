#ifndef STAN_MATH_TORSTEN_PKMODEL_HPP
#define STAN_MATH_TORSTEN_PKMODEL_HPP

#include <stan/math/prim/scal/err/check_positive_finite.hpp>
#include <stan/math/torsten/PKModel/pmetricsCheck.hpp>
#include <stan/math/torsten/PKModel/functions.hpp>
#include <stan/math/torsten/PKModel/SearchReal.hpp>
#include <stan/math/torsten/PKModel/ExtractVector.hpp>
#include <stan/math/torsten/PKModel/PKModel_class.hpp>
#include <stan/math/torsten/PKModel/Event.hpp>
#include <stan/math/torsten/PKModel/Rate.hpp>
#include <stan/math/torsten/PKModel/ModelParameters.hpp>
#include <stan/math/torsten/PKModel/pmetrics_solver.hpp>
#include <stan/math/torsten/PKModel/Pred/PolyExp.hpp>

// Initial functor for ODE solver
pmetrics_solver_structure pmetrics_solver; 

#include <stan/math/torsten/PKModel/Pred1.hpp>
#include <stan/math/torsten/PKModel/PredSS.hpp>

// Initial functor for evolution operators
Pred1_structure Pred1("default");
PredSS_structure PredSS("default");

#include <stan/math/torsten/PKModel/Pred.hpp>

extern int marker_count; // For testing purposes

#endif
