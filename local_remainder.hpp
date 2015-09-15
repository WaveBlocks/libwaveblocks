#pragma once
#include "macros.hpp"
#include "abstract.hpp"
template <class EvalImpl,
         template<int,int> class Basis,
         int N,
         int D >
class LocalRemainderImplementation : public LocalRemainderAbstract<LocalRemainderImplementation<EvalImpl, Basis, N, D>,Basis,N,D>, public EvalImpl
{
  public:
    IMPORT_TYPES_FROM(Basis,N,D);
    
  protected:
    local_remainder_type local_remainder;
    local_quadratic_type local_quadratic;
    
  public:
    LocalRemainderImplementation(potential_type potential,
                              jacobian_type jacobian,
                              hessian_type hessian )
      : EvalImpl(potential,jacobian,hessian)
       {
      calculate_local_quadratic_implementation();
      calculate_local_remainder_implementation();
    }
    
    void calculate_local_quadratic_implementation() {
      local_quadratic = [ = ]( RVector<D> q ) {
        potential_type result_matrix;
        for ( int l = 0; l < N; ++l ) {
          for ( int m = 0; m < N; ++m ) {
            result_matrix( l, m ) = [ = ]( RVector<D> x ) {
              auto V = EvalImpl::potential( l, m )( q );
              auto J = EvalImpl::jacobian( l, m )( q );
              auto H = EvalImpl::hessian( l, m )( q );
              
              auto result = V;
              
              for ( int i = 0; i < D; ++i ) {
                auto xmqi = x[i] - q[i];
                result += J[i] * ( xmqi );
                
                for ( int j = 0; j < D; ++j ) {
                  result += 0.5 * xmqi * H( i, j ) * ( x[j] - q[j] );
                }
              }
              
              return result;
            };
          }
        }
        
        return result_matrix;
      };
    }
    
    void calculate_local_remainder_implementation() {
      local_remainder = [ = ]( RVector<D> q ) {
        potential_type result;
        auto local_quadratic_q = local_quadratic( q );
        
        for ( int i = 0; i < N; ++i ) {
          for ( int j = 0; j < N; ++j ) {
            result( i, j ) = [ = ]( RVector<D> x ) {
              return EvalImpl::potential( i, j )( x ) - local_quadratic_q( i, j )( x );
            };
          }
        }
        
        return result;
      };
    }
    
};
