#pragma once
#include "macros.hpp"
#include "Bases.hpp"

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
  Tuple<> taylor_at(const RVector<D>& g ) {
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


template < 
         template<int,int> class  Basis,
         int N,
         int D >
class EvaluationImplementation;

template <int N, int D>
struct EvaluationImplementation<CanonicalBasis, N, D> : EvaluationAbstract < EvaluationImplementation<CanonicalBasis, N, D>,CanonicalBasis,N,D> {
  
    IMPORT_TYPES_FROM(CanonicalBasis,N,D);
    
    
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
             function_t > ( potential, arg );
    }
    
    jacobian_evaluation_type evaluate_jacobian_at_implementation(
      const RVector<D> &arg ) {
      return evaluate_function_matrix < N,
             GMatrix,
             RVector<D>,
             jacobian_return_type,
             function_t > ( jacobian, arg );
    }
    
    hessian_evaluation_type evaluate_hessian_at_implementation(
      const RVector<D> &arg ) {
      return evaluate_function_matrix < N,
             GMatrix,
             RVector<D>,
             hessian_return_type,
             function_t > ( hessian, arg );
    }
    

};

template <int N, int D>
struct EvaluationImplementation <
    EigenBasis,
    N,
    D > : EvaluationAbstract< EvaluationImplementation<EigenBasis,N,D>, EigenBasis, N, D> {
      
    IMPORT_TYPES_FROM(EigenBasis,N,D);
    
    
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
             function_t > ( potential, arg );
    }
    
    jacobian_evaluation_type evaluate_jacobian_at_implementation(
      const RVector<D> &arg ) {
      return evaluate_function_vector < N,
             GVector,
             RVector<D>,
             jacobian_return_type,
             function_t > ( jacobian, arg );
    }
    
    hessian_evaluation_type evaluate_hessian_at_implementation(
      const RVector<D> &arg ) {
      return evaluate_function_vector < N,
             GVector,
             RVector<D>,
             hessian_return_type,
             function_t > ( hessian, arg );
    }
};
