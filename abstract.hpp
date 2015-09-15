#pragma once
#include "macros.hpp"
#include "bases.hpp"


// Abstract classes with common functionality (grid evaluations)
// CRTP
template <class Subtype,
         template<int,int> class Basis,
         int N,
         int D >
struct LocalRemainderAbstract
{
  using Self = LocalRemainderAbstract<Subtype,Basis,N,D>;
  IMPORT_TYPES_FROM(Basis,N,D);
  
  void calculate_local_quadratic() {
    that.calculate_local_quadratic_implementation();
  }
  
  void calculate_local_remainder() {
    that.calculate_local_quadratic_implementation();
  }
  
   potential_evaluation_type evaluate_local_remainder_at( RVector<D> g,
      RVector<D> position ) {
    return that.evaluate_local_remainder_at_implementation( g, position );
  }
  
  template < template <typename...> class grid_in = std::vector,
           template <typename...> class grid_out = grid_in >
  grid_out<potential_evaluation_type> evaluate_local_remainder(
    const grid_in<RVector<D> > &args,
    RVector<D> position ) {
    return evaluate_function_in_grid < RVector<D>,
           potential_evaluation_type,
           grid_in,
           grid_out,
           function_t > (
             std::bind( &Self::evaluate_local_remainder_at,
                        this,
                        std::placeholders::_1,
                        position ),
             args );
  }
};


template <class Subtype,
         template<int,int> class Basis,
         int N,
         int D >
struct EvaluationAbstract
{
  using Self = EvaluationAbstract<Subtype,Basis,N,D>;
  IMPORT_TYPES_FROM(Basis,N,D);
  
  potential_evaluation_type evaluate_at( const RVector<D> &arg ) {
    return that.evaluate_at_implementation( arg );
  }
  
  template < template <typename...> class grid_in = std::vector,
           template <typename...> class grid_out = grid_in >
  grid_out<potential_evaluation_type> evaluate(
    const grid_in<RVector<D> > &args ) {
    return evaluate_function_in_grid < RVector<D>,
           potential_evaluation_type,
           grid_in,
           grid_out,
           function_t > (
             std::bind( &Self::evaluate_at, this, std::placeholders::_1 ), args );
  }
  
  jacobian_evaluation_type evaluate_jacobian_at( const RVector<D> &arg ) {
    return that.evaluate_jacobian_at_implementation( arg );
  }
  
  template < template <typename...> class grid_in = std::vector,
           template <typename...> class grid_out = grid_in >
  grid_out<jacobian_evaluation_type> evaluate_jacobian(
    const grid_in<RVector<D> > &args ) {
    return evaluate_function_in_grid < RVector<D>,
           jacobian_evaluation_type,
           grid_in,
           grid_out,
           function_t > (
             std::bind( &Self::evaluate_jacobian_at, this, std::placeholders::_1 ),
             args );
  }
  
  hessian_evaluation_type evaluate_hessian_at( const RVector<D> &arg ) {
    return that.evaluate_hessian_at_implementation( arg );
  }
  
  template < template <typename...> class grid_in = std::vector,
           template <typename...> class grid_out = grid_in >
  grid_out<hessian_evaluation_type> evaluate_hessian(
    const grid_in<RVector<D> > &args ) {
    return evaluate_function_in_grid < RVector<D>,
           hessian_evaluation_type,
           grid_in,
           grid_out,
           function_t > (
             std::bind( &Self::evaluate_hessian_at, this, std::placeholders::_1 ),
             args );
  }
  
  template <template <typename...> class Tuple = std::tuple>
  Tuple<> taylor_at( RVector<D> g ) {
    return Tuple<>(
             evaluate_at( g ), evaluate_jacobian_at( g ), evaluate_hessian_at( g ) );
  }
  
  template < template <typename...> class Tuple = std::tuple,
           template <typename...> class grid_in = std::vector,
           template <typename...> class grid_out = grid_in >
  grid_out<Tuple<> > taylor( const grid_in<RVector<D> > &args ) {
    return evaluate_function_in_grid < RVector<D>,
           Tuple<>,
           grid_in,
           grid_out,
           function_t > (
             std::bind( &Self::taylor_at, this, std::placeholders::_1 ), args );
  }
};

