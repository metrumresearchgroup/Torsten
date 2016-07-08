// version 0.8
// Headerfiles for the Torsten Library
// Order matters.

#ifndef PKMODEL_HPP
#define PKMODEL_HPP

#include <stan/model/model_header.hpp>

#include "functions.hpp"
#include "SearchReal.hpp"
#include "ExtractVector.hpp"

#include "PKModel_class.hpp"
#include "Event.hpp"
#include "Rate.hpp"
#include "ModelParameters.hpp"

#include "pmetrics_solver.hpp"
pmetrics_solver_structure pmetrics_solver; // Define functors for ODE solvers
										                       // as global variables

#include "Pred/PolyExp.hpp"
#include "Pred1.hpp"
#include "PredSS.hpp"

// Define Pred1 and PredSS as global variables
Pred1_structure Pred1("default");
PredSS_structure PredSS("default");

#include "Pred.hpp"

// For testing purposes, declare global variable marker_counter
extern int marker_count; 

#endif
