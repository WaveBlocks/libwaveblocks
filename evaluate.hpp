#pragma once
#include "macros.hpp"
#include "bases.hpp"

template < template <template <int, int> class, int, int> class S,
         template <int, int> class B,
         int N,
         int D >
struct EvaluationAbstract
{
  GET_TYPES( EvaluationAbstract, S, B, N, D );
  
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
             std::bind( &self_type::evaluate_at, this, std::placeholders::_1 ), args );
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
             std::bind(
               &self_type::evaluate_jacobian_at, this, std::placeholders::_1 ),
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
             std::bind( &self_type::evaluate_hessian_at, this, std::placeholders::_1 ),
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
             std::bind( &self_type::taylor_at, this, std::placeholders::_1 ), args );
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
             std::bind( &self_type::evaluate_local_remainder_at,
                        this,
                        std::placeholders::_1,
                        position ),
             args );
  }
};

template < template <template <int, int> class, int, int> class S,
         template <int, int> class B,
         int N,
         int D >
class EvaluationImplementation;

template <template <template <int, int> class, int, int> class S, int N, int D>
struct EvaluationImplementation<S, CanonicalBasis, N, D> : EvaluationAbstract <
    S,
    CanonicalBasis,
    N,
    D > {
    GET_TYPES( EvaluationImplementation, S, CanonicalBasis, N, D );
    
    friend Subtype;
    
  protected:
    potential_type potential;
    jacobian_type jacobian;
    hessian_type hessian;
    
  public:
    EvaluationImplementation( potential_type potential,
                              jacobian_type jacobian,
                              hessian_type hessian )
      : potential( potential ), jacobian( jacobian ), hessian( hessian ) {}
      
  public:
    potential_evaluation_type evaluate_at_implementation( const RVector<D> &arg ) {
      return evaluate_function_matrix < N,
             GMatrix,
             RVector<D>,
             potential_return_type,
             function_t > ( this->potential, arg );
    }
    
    jacobian_evaluation_type evaluate_jacobian_at_implementation(
      const RVector<D> &arg ) {
      return evaluate_function_matrix < N,
             GMatrix,
             RVector<D>,
             jacobian_return_type,
             function_t > ( this->jacobian, arg );
    }
    
    hessian_evaluation_type evaluate_hessian_at_implementation(
      const RVector<D> &arg ) {
      return evaluate_function_matrix < N,
             GMatrix,
             RVector<D>,
             hessian_return_type,
             function_t > ( this->hessian, arg );
    }
    
    potential_evaluation_type evaluate_local_remainder_at_implementation(
      RVector<D> g,
      RVector<D> position ) {
      return evaluate_function_matrix < N,
             GMatrix,
             RVector<D>,
             potential_return_type,
             function_t > ( that.local_remainder( position ),
                            g );
    }
};

template <template <template <int, int> class, int, int> class S, int N, int D>
struct EvaluationImplementation < S,
    EigenBasis,
    N,
    D > : EvaluationAbstract<S, EigenBasis, N, D> {
    GET_TYPES( EvaluationImplementation, S, EigenBasis, N, D );
    
    friend Subtype;
    
  protected:
    potential_type potential;
    jacobian_type jacobian;
    hessian_type hessian;
    
  public:
    EvaluationImplementation( potential_type potential,
                              jacobian_type jacobian,
                              hessian_type hessian )
      : potential( potential ), jacobian( jacobian ), hessian( hessian ) {}
      
  public:
    potential_evaluation_type evaluate_at_implementation( const RVector<D> &arg ) {
      return evaluate_function_vector < N,
             GVector,
             RVector<D>,
             potential_return_type,
             function_t > ( this->potential, arg );
    }
    
    jacobian_evaluation_type evaluate_jacobian_at_implementation(
      const RVector<D> &arg ) {
      return evaluate_function_vector < N,
             GVector,
             RVector<D>,
             jacobian_return_type,
             function_t > ( this->jacobian, arg );
    }
    
    hessian_evaluation_type evaluate_hessian_at_implementation(
      const RVector<D> &arg ) {
      return evaluate_function_vector < N,
             GVector,
             RVector<D>,
             hessian_return_type,
             function_t > ( this->hessian, arg );
    }
    
    potential_evaluation_type evaluate_local_remainder_at_implementation(
      RVector<D> g,
      RVector<D> position ) {
      return evaluate_function_vector < N,
             GVector,
             RVector<D>,
             potential_return_type,
             function_t > ( this->local_remainder( position ),
                            g );
    }
};
