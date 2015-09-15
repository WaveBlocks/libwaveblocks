#pragma once
#include "macros.hpp"
#include "bases.hpp"
#include "abstract.hpp"

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
