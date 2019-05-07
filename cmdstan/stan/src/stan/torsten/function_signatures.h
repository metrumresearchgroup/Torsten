// included from constructor for function_signatures() in src/stan/lang/ast.hpp

/****************************************
 * TORSTEN: function signatures
 ****************************************/

/* TIME, AMT,, RATE, II, EVID, CMT, ADDL, SS */
const int num_pmx_args = 11;
std::vector<function_arg_type> pmx_data_arg_types(num_pmx_args);
for (int i = 0; i < 4; i++) pmx_data_arg_types[i] = function_arg_type(vector_types[1]);
for (int i = 0; i < 4; i++) pmx_data_arg_types[i + 4] = function_arg_type(int_vector_types[1]);

/* THETA, BIOVAR, TLAG */
#define TORSTEN_PMX_FUNC_ARG_TYPES_TABLE \
  TORSTEN_PMX_ARGS(1U, 1U, 1U)           \
  TORSTEN_PMX_ARGS(2U, 1U, 1U)           \
  TORSTEN_PMX_ARGS(1U, 2U, 1U)           \
  TORSTEN_PMX_ARGS(2U, 2U, 1U)           \
  TORSTEN_PMX_ARGS(1U, 1U, 2U)           \
  TORSTEN_PMX_ARGS(2U, 1U, 2U)           \
  TORSTEN_PMX_ARGS(1U, 2U, 2U)           \
  TORSTEN_PMX_ARGS(2U, 2U, 2U)

#define TORSTEN_PMX_ARGS(A, B, C) pmx_data_arg_types[8] = function_arg_type(expr_type(double_type(), A)); \
  pmx_data_arg_types[9]  = function_arg_type(expr_type(double_type(), B)); \
  pmx_data_arg_types[10] = function_arg_type(expr_type(double_type(), C)); \
  add("PKModelOneCpt", expr_type(matrix_type()), pmx_data_arg_types); \
  add("PKModelTwoCpt", expr_type(matrix_type()), pmx_data_arg_types); \
  add("pmx_solve_onecpt", expr_type(matrix_type()), pmx_data_arg_types); \
  add("pmx_solve_twocpt", expr_type(matrix_type()), pmx_data_arg_types);
    TORSTEN_PMX_FUNC_ARG_TYPES_TABLE
#undef TORSTEN_PMX_ARGS

#undef TORSTEN_PMX_FUNC_ARG_TYPES_TABLE

#define TORSTEN_PMX_FUNC_ARG_TYPES_TABLE                                    \
  TORSTEN_PMX_ARGS(function_arg_type(expr_type(matrix_type(), 1U)), 1U, 1U) \
  TORSTEN_PMX_ARGS(function_arg_type(expr_type(matrix_type()))    , 1U, 1U) \
  TORSTEN_PMX_ARGS(function_arg_type(expr_type(matrix_type(), 1U)), 2U, 1U) \
  TORSTEN_PMX_ARGS(function_arg_type(expr_type(matrix_type()))    , 2U, 1U) \
  TORSTEN_PMX_ARGS(function_arg_type(expr_type(matrix_type(), 1U)), 1U, 2U) \
  TORSTEN_PMX_ARGS(function_arg_type(expr_type(matrix_type()))    , 1U, 2U) \
  TORSTEN_PMX_ARGS(function_arg_type(expr_type(matrix_type(), 1U)), 2U, 2U) \
  TORSTEN_PMX_ARGS(function_arg_type(expr_type(matrix_type()))    , 2U, 2U)

#define TORSTEN_PMX_ARGS(A, B, C) pmx_data_arg_types[8] = A; \
  pmx_data_arg_types[9]  = function_arg_type(expr_type(double_type(), B)); \
  pmx_data_arg_types[10] = function_arg_type(expr_type(double_type(), C)); \
  add("linOdeModel", expr_type(matrix_type()), pmx_data_arg_types); \
  add("pmx_solve_linode", expr_type(matrix_type()), pmx_data_arg_types);
    TORSTEN_PMX_FUNC_ARG_TYPES_TABLE
#undef TORSTEN_PMX_ARGS

#undef TORSTEN_PMX_FUNC_ARG_TYPES_TABLE

add("linear_interpolation", expr_type(double_type()), expr_type(double_type()), vector_types[1], vector_types[1]);
add("linear_interpolation", vector_types[1], vector_types[1], vector_types[1], vector_types[1]);
