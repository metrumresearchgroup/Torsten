#ifndef STAN_MATH_TORSTEN_PKMODEL_FUNCTORS_DUMMY_ODE_HPP
#define STAN_MATH_TORSTEN_PKMODEL_FUNCTORS_DUMMY_ODE_HPP

#include<vector>

/**
 * A structure for the functor that is passed in the ODE 
 * integrator. Can be used to specify a dummy ODE functor,
 * which allows functions not using an ODE integrator to 
 * call pred().
 */
struct dummy_ode {
  template <typename T0, typename T1, typename T2, typename T3>
  inline
  std::vector<typename boost::math::tools::promote_args<T0, T1, T2, T3>::type>
  rate_dbl(const T0& t,
           const std::vector<T1>& x,
           const std::vector<T2>& parms,
           const std::vector<T3>& rate,
           const std::vector<int>& dummy,
           std::ostream* pstream__) const {
              typedef typename boost::math::tools::
                promote_args<T0, T1, T2, T3>::type scalar;
    std::vector<scalar> returned_vector = std::vector<scalar>(0, scalar(0));
    return returned_vector;
  }

  template <typename T0, typename T1, typename T2, typename T3>
  inline
  std::vector<typename boost::math::tools::promote_args<T0, T1, T2, T3>::type>
  rate_var(const T0& t,
           const std::vector<T1>& x,
           const std::vector<T2>& parms,
           const std::vector<T3>& rate,
           const std::vector<int>& dummy,
           std::ostream* pstream__) const {
              typedef typename boost::math::tools::
                promote_args<T0, T1, T2, T3>::type scalar;
    std::vector<scalar> returned_vector = std::vector<scalar>(0, scalar(0));
    return returned_vector;
  }

};

#endif
