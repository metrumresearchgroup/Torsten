#ifndef TORSTEN_ODE_TAGGED_FUNCTOR_HPP
#define TORSTEN_ODE_TAGGED_FUNCTOR_HPP

struct torsten_ode_tagged_functor {
  const int tag;
  template<typename F>
  torsten_ode_tagged_functor(const F& f0) :
    tag(torsten_ode_tagged_functor::get_tag(f0))
  {}

  template<typename F>
  torsten_ode_tagged_functor() :
    tag(torsten_ode_tagged_functor::get_tag<F>())
  {}

  template<typename... Ts>
  auto operator()(const Ts... xs) {
#define ODE_FUNC(i, f) if (tag == i) return f()(xs...);
    TORSTEN_ODE_FUNC_TABLE
#undef ODE_FUNC
  }

  template<typename F>
  static constexpr int get_tag(const F& f0) {
#define ODE_FUNC(i, f) if (typeid(F) == typeid(f)) return i;
    TORSTEN_ODE_FUNC_TABLE
#undef ODE_FUNC
      return -1;
  }

  template<typename F>
  static constexpr int get_tag() {
#define ODE_FUNC(i, f) if (typeid(F) == typeid(f)) return i;
    TORSTEN_ODE_FUNC_TABLE
#undef ODE_FUNC
      return -1;
  }
};

#endif
