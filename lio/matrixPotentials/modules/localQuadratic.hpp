#pragma once
#include "macros.hpp"
#include "types.hpp"
#include "utilities/evaluations.hpp"

namespace waveblocks
{
  namespace matrixPotentials
  {
    namespace modules
    {
      namespace localQuadratic
      {
         /**
         * \brief Abstract class for local quadratic evaluation
         * 
         * A matrix potential inheriting an implementation of this module
         * can evaluate the local quadratic approximation of its' potential
         * elementwise
         * 
         * This makes use of the CRTPattern
         * 
         * \tparam Subtype The type extending this interface (used for static polymorphism)
         * \tparam Basis
         * Which basis (bases::Eigen or bases::Canonical) the potential is given in
         * \tparam N
         * Number of levels (dimension of square matrix when evaluated)
         * \tparam D
         * Dimension of argument space
         */
        template <class Subtype, template <int, int> class Basis, int N, int D>
        struct Abstract {
            using Self = Abstract<Subtype, Basis, N, D>;
            IMPORT_TYPES_FROM( Basis, N, D );
            
          protected:
            local_quadratic_type local_quadratic;
            
          public:
            void calculate_local_quadratic() const {
              static_cast<const Subtype*>(this)->calculate_local_quadratic_implementation();
            }
            
            potential_evaluation_type evaluate_local_quadratic_at(
              const CVector<D> &arg,
              const CVector<D> &position ) const {
              return static_cast<const Subtype*>(this)->local_quadratic( position )( arg );
            }
            
            template < template <typename...> class grid_in = std::vector,
                     template <typename...> class grid_out = grid_in >
            grid_out<potential_evaluation_type> evaluate_local_remainder(
              const grid_in<CVector<D> > &args,
              CVector<D> position ) const {
              return utilities::evaluate_function_in_grid < CVector<D>,
                     potential_evaluation_type,
                     grid_in,
                     grid_out,
                     function_t > (
                       std::bind( &Self::evaluate_local_quadratic_at,
                                  this,
                                  std::placeholders::_1,
                                  position ),
                       args );
            }
        };
        
        template <class EvalImpl, template <int, int> class Basis, int N, int D>
        class Standard : public Abstract<Standard<EvalImpl, Basis, N, D>, Basis, N, D>,
          public EvalImpl
        {
          public:
            IMPORT_TYPES_FROM( Basis, N, D );
            using Super = Abstract<Standard<EvalImpl, Basis, N, D>, Basis, N, D>;
            
          public:
            Standard( potential_type potential,
                      jacobian_type jacobian,
                      hessian_type hessian )
              : EvalImpl( potential, jacobian, hessian ) {
              calculate_local_quadratic_implementation();
            }
            
            void calculate_local_quadratic_implementation() const {
              Super::local_quadratic = [ this ]( CVector<D> q ) {
                potential_type result_matrix;
                
                for ( int l = 0; l < N; ++l ) {
                  for ( int m = 0; m < N; ++m )  {
                    result_matrix( l, m ) = [ this, l, m, q]( CVector<D> x ) {
                  // Takes care of
                  // http://stackoverflow.com/questions/19850648/error-when-calling-base-member-function-from-within-a-lambda
#if defined(__clang__)
                  auto V = EvalImpl::potential( q );
                  auto J = EvalImpl::jacobian( q );
                  auto H = EvalImpl::hessian( q );
#elif defined(__GNUC__) || defined(__GNUG__)
                  auto V = this->potential( q );
                  auto J = this->jacobian( q );
                  auto H = this->hessian( q );
#else
                  auto V = EvalImpl::potential( q );
                  auto J = EvalImpl::jacobian( q );
                  auto H = EvalImpl::hessian( q );
#endif                      
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
        };
        
        template <class EvalImpl, template <int, int> class Basis, int D>
        class Standard<EvalImpl, Basis, 1, D> : public Abstract <
          Standard<EvalImpl, Basis, 1, D>,
          Basis,
          1,
          D > ,
        public EvalImpl
        {
          public:
            IMPORT_TYPES_FROM( Basis, 1, D );
            using Super = Abstract<Standard<EvalImpl, Basis, 1, D>, Basis, 1, D>;
            
          public:
            Standard( potential_type potential,
                      jacobian_type jacobian,
                      hessian_type hessian )
              : EvalImpl( potential, jacobian, hessian ) {
              calculate_local_quadratic_implementation();
            }
            
            void calculate_local_quadratic_implementation() {
              Super::local_quadratic = [ this ]( CVector<D> q ) {
                return [ this ,q ]( CVector<D> x ) {
                
                  // Takes care of
                  // http://stackoverflow.com/questions/19850648/error-when-calling-base-member-function-from-within-a-lambda
#if defined(__clang__)
                  auto V = EvalImpl::potential( q );
                  auto J = EvalImpl::jacobian( q );
                  auto H = EvalImpl::hessian( q );
#elif defined(__GNUC__) || defined(__GNUG__)
                  auto V = this->potential( q );
                  auto J = this->jacobian( q );
                  auto H = this->hessian( q );
#else
                  auto V = EvalImpl::potential( q );
                  auto J = EvalImpl::jacobian( q );
                  auto H = EvalImpl::hessian( q );
#endif
                  
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
              };
            }
        };
      }
      template <class EvalImpl, template <int, int> class Basis, int N, int D>
      using LocalQuadratic = localQuadratic::Standard<EvalImpl, Basis, N, D>;
    }
  }
}
