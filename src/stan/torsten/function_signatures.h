// included from constructor for function_signatures() in src/stan/lang/ast.hpp

/****************************************
 * TORSTEN: function signatures
 ****************************************/

/* arg_types.push_back(function_arg_type(arg_type)) */

std::vector<function_arg_type> data_arg_types;
for (int i = 0; i < 4; i++)
  data_arg_types.push_back(function_arg_type(vector_types[1]));
for (int i = 0; i < 4; i++)
  data_arg_types.push_back(function_arg_type(int_vector_types[1]));

std::vector<function_arg_type> arg_types_222 = data_arg_types;
for (int i = 0; i < 3; i++)
  arg_types_222.push_back(function_arg_type(expr_type(double_type(), 2U)));

std::vector<function_arg_type> arg_types_122 = data_arg_types;
arg_types_122.push_back(function_arg_type(expr_type(double_type(), 1U)));
arg_types_122.push_back(function_arg_type(expr_type(double_type(), 2U)));
arg_types_122.push_back(function_arg_type(expr_type(double_type(), 2U)));

std::vector<function_arg_type> arg_types_112 = data_arg_types;
arg_types_112.push_back(function_arg_type(expr_type(double_type(), 1U)));
arg_types_112.push_back(function_arg_type(expr_type(double_type(), 1U)));
arg_types_112.push_back(function_arg_type(expr_type(double_type(), 2U)));

std::vector<function_arg_type> arg_types_111 = data_arg_types;
for (int i = 0; i < 3; i++)
  arg_types_111.push_back(function_arg_type(expr_type(double_type(), 1U)));

std::vector<function_arg_type> arg_types_121 = data_arg_types;
arg_types_121.push_back(function_arg_type(expr_type(double_type(), 1U)));
arg_types_121.push_back(function_arg_type(expr_type(double_type(), 2U)));
arg_types_121.push_back(function_arg_type(expr_type(double_type(), 1U)));

std::vector<function_arg_type> arg_types_212 = data_arg_types;
arg_types_212.push_back(function_arg_type(expr_type(double_type(), 2U)));
arg_types_212.push_back(function_arg_type(expr_type(double_type(), 1U)));
arg_types_212.push_back(function_arg_type(expr_type(double_type(), 2U)));

std::vector<function_arg_type> arg_types_211 = data_arg_types;
arg_types_211.push_back(function_arg_type(expr_type(double_type(), 2U)));
arg_types_211.push_back(function_arg_type(expr_type(double_type(), 1U)));
arg_types_211.push_back(function_arg_type(expr_type(double_type(), 1U)));

std::vector<function_arg_type> arg_types_221 = data_arg_types;
arg_types_221.push_back(function_arg_type(expr_type(double_type(), 2U)));
arg_types_221.push_back(function_arg_type(expr_type(double_type(), 2U)));
arg_types_221.push_back(function_arg_type(expr_type(double_type(), 1U)));

add("PKModelOneCpt", expr_type(matrix_type()), arg_types_222);
add("PKModelOneCpt", expr_type(matrix_type()), arg_types_122);
add("PKModelOneCpt", expr_type(matrix_type()), arg_types_112);
add("PKModelOneCpt", expr_type(matrix_type()), arg_types_111);
add("PKModelOneCpt", expr_type(matrix_type()), arg_types_121);
add("PKModelOneCpt", expr_type(matrix_type()), arg_types_212);
add("PKModelOneCpt", expr_type(matrix_type()), arg_types_211);
add("PKModelOneCpt", expr_type(matrix_type()), arg_types_221);

add("PKModelTwoCpt", expr_type(matrix_type()), arg_types_222);
add("PKModelTwoCpt", expr_type(matrix_type()), arg_types_122);
add("PKModelTwoCpt", expr_type(matrix_type()), arg_types_112);
add("PKModelTwoCpt", expr_type(matrix_type()), arg_types_111);
add("PKModelTwoCpt", expr_type(matrix_type()), arg_types_121);
add("PKModelTwoCpt", expr_type(matrix_type()), arg_types_212);
add("PKModelTwoCpt", expr_type(matrix_type()), arg_types_211);
add("PKModelTwoCpt", expr_type(matrix_type()), arg_types_221);

arg_types_222[8] = function_arg_type(expr_type(matrix_type(), 1U));
arg_types_122[8] = function_arg_type(expr_type(matrix_type()));
arg_types_112[8] = function_arg_type(expr_type(matrix_type()));
arg_types_111[8] = function_arg_type(expr_type(matrix_type()));
arg_types_121[8] = function_arg_type(expr_type(matrix_type()));
arg_types_212[8] = function_arg_type(expr_type(matrix_type(), 1U));
arg_types_211[8] = function_arg_type(expr_type(matrix_type(), 1U));
arg_types_221[8] = function_arg_type(expr_type(matrix_type(), 1U));

add("linOdeModel", expr_type(matrix_type()), arg_types_222);
add("linOdeModel", expr_type(matrix_type()), arg_types_122);
add("linOdeModel", expr_type(matrix_type()), arg_types_112);
add("linOdeModel", expr_type(matrix_type()), arg_types_111);
add("linOdeModel", expr_type(matrix_type()), arg_types_121);
add("linOdeModel", expr_type(matrix_type()), arg_types_212);
add("linOdeModel", expr_type(matrix_type()), arg_types_211);
add("linOdeModel", expr_type(matrix_type()), arg_types_221);

add("linear_interpolation", expr_type(double_type()), expr_type(double_type()), vector_types[1], vector_types[1]);
add("linear_interpolation", vector_types[1], vector_types[1], vector_types[1], vector_types[1]);
