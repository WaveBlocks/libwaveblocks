#pragma once
#include "macros.hpp"
#include "matrixPotentials/bases.hpp"
#include "matrixPotentials/modules/localQuadratic.hpp"
#include "types.hpp"
#include "utilities/evaluations.hpp"

namespace waveblocks
{
  namespace matrixPotentials
  {
    namespace modules
    {
      namespace localRemainder
      {
        namespace helper
        {
        
          template <int N>
          struct DiagonalDifference {
            struct Inhomogenous {
              static CMatrix<N, N> apply( RMatrix<N, N> V, const CVector<N> &u ) {
                CMatrix<N, N> C = V.template cast<complex_t>();
                
                for ( int i = 0; i < N; ++i ) {
                  C( i, i ) -= u( i );
                }
                
                return C;
              }
            };
            struct Homogenous {
              static CMatrix<N, N> apply( RMatrix<N, N> V, const complex_t &u ) {
                CMatrix<N, N> C = V.template cast<complex_t>();
                
                for ( int i = 0; i < N; ++i ) {
                  C( i, i ) -= u;
                }
                
                return C;
              }
            };
          };
          
          template <>
          struct DiagonalDifference<1> {
            struct Homogenous {
              static complex_t apply( real_t V, const complex_t &u ) {
                return V - u;
              }
            };
            using Inhomogenous = Homogenous;
          };
        }
        
        template <class Subtype, int N, int D>
        struct Abstract {
          using Self = Abstract<Subtype, N, D>;
          CMatrix<N, N> evaluate_local_remainder_at( const RVector<D> &arg,
              const RVector<D> &q ) const {
            return static_cast<const Subtype*>(this)->evaluate_local_remainder_at_implementation( arg, q );
          }
          
          template < template <typename...> class grid_in = std::vector,
                   template <typename...> class grid_out = grid_in >
          grid_out<CMatrix<N, N> > evaluate_local_remainder(
            const grid_in<RVector<D> > &args,
            const RVector<D> &q ) const {
            return utilities::evaluate_function_in_grid < RVector<D>,
                   CMatrix<N, N>,
                   grid_in,
                   grid_out,
                   function_t > (
                     std::bind(
                       &Self::evaluate_local_remainder_at, this, std::placeholders::_1, q ),
                     args );
          }
          
        }; 
        
        template <class EvalImpl, class LeadingEvalImpl, int N, int D>
        class Homogenous
          : public Abstract<Homogenous<EvalImpl, LeadingEvalImpl, N, D>, N, D>,
            EvalImpl
        {
            IMPORT_TYPES_FROM( bases::Canonical, N, D );
                      
            LocalQuadratic<LeadingEvalImpl, bases::Eigen, 1, D> quadratic;
            
          public:
                      using leading_level_type = bases::Eigen<1, D>;

            Homogenous( potential_type pot,
                        jacobian_type jac,
                        hessian_type hess,
                        typename leading_level_type::potential_type lead_pot,
                        typename leading_level_type::jacobian_type lead_jac,
                        typename leading_level_type::hessian_type lead_hess )
              : EvalImpl( pot, jac, hess ), quadratic( lead_pot, lead_jac, lead_hess ) {}
              
            CMatrix<N, N> evaluate_local_remainder_at_implementation(
              const RVector<D> &arg,
              const RVector<D> &q ) const {
              auto u = quadratic.evaluate_local_quadratic_at( arg, q );
              return helper::DiagonalDifference<N>::Homogenous::apply(
                       EvalImpl::evaluate_at( arg ), u );
            }
            
          template <template <typename...> class Tuple = std::tuple>
          Tuple<typename leading_level_type::potential_evaluation_type,
          typename leading_level_type::jacobian_evaluation_type,
          typename leading_level_type::hessian_evaluation_type> taylor_leading_at( const RVector<D> &g ) const {
            return quadratic.taylor_at(g);
          }
        };
        
        template <class EvalImpl, class LeadingEvalImpl, int N, int D>
        class Inhomogenous
          : public Abstract<Inhomogenous<EvalImpl, LeadingEvalImpl, N, D>, N, D>,
            EvalImpl
        {
            IMPORT_TYPES_FROM( bases::Canonical, N, D );
            
            
            LocalQuadratic<LeadingEvalImpl, bases::Eigen, N, D> quadratic;
            
          public:
                      using leading_level_type = bases::Eigen<N, D>;

            Inhomogenous( potential_type pot,
                          jacobian_type jac,
                          hessian_type hess,
                          typename leading_level_type::potential_type lead_pot,
                          typename leading_level_type::jacobian_type lead_jac,
                          typename leading_level_type::hessian_type lead_hess )
              : EvalImpl( pot, jac, hess ), quadratic( lead_pot, lead_jac, lead_hess ) {}
              
            CMatrix<N, N> evaluate_local_remainder_at_implementation(
              const RVector<D> &arg,
              const RVector<D> &q ) const {
              auto u = quadratic.evaluate_local_quadratic_at( arg, q );
              return helper::DiagonalDifference<N>::Inhomogenous::apply(
                       EvalImpl::evaluate_at( arg ), u );
            }
            
            template <template <typename...> class Tuple = std::tuple>
              Tuple<typename leading_level_type::potential_evaluation_type,
          typename leading_level_type::jacobian_evaluation_type,
          typename leading_level_type::hessian_evaluation_type> taylor_leading_at( const RVector<D> &g ) const {
            return quadratic.taylor_at(g);
          }
        };
      }
    }
  }
}