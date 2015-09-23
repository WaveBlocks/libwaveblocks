#pragma once
#include "macros.hpp"
#include "matrixPotentials/bases.hpp"
#include "matrixPotentials/modules/localQuadratic.hpp"
#include "matrixPotentials/modules/leadingLevelOwner.hpp"
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
         
        template <class Subtype, int N, int D>
        struct Abstract {
          using Self = Abstract<Subtype, N, D>;
          IMPORT_TYPES_FROM( bases::Canonical, N, D );


          local_quadratic_evaluation_type evaluate_local_remainder_at( const CVector<D> &arg,
              const CVector<D> &q ) const {
            return static_cast<const Subtype*>(this)->evaluate_local_remainder_at_implementation( arg, q );
          }
          
          template < template <typename...> class grid_in = std::vector,
                   template <typename...> class grid_out = grid_in >
          grid_out<local_quadratic_evaluation_type > evaluate_local_remainder(
            const grid_in<CVector<D> > &args,
            const CVector<D> &q ) const {
            return utilities::evaluate_function_in_grid < CVector<D>,
                   local_quadratic_evaluation_type,
                   grid_in,
                   grid_out,
                   function_t > (
                     std::bind(
                       &Self::evaluate_local_remainder_at, this, std::placeholders::_1, q ),
                     args );
          }
          
        }; 
        
        template <class DiagDifference, class EvalImpl, class LocQuadraticImpl, int N, int D>
        class General : public Abstract<General<DiagDifference, EvalImpl, LocQuadraticImpl, N, D>, N, D>, public EvalImpl, public LeadingLevelOwner<LocQuadraticImpl>
        {
            IMPORT_TYPES_FROM( bases::Canonical, N, D );
                      
          public:

            General( potential_type pot,
                        typename LocQuadraticImpl::potential_type lead_pot,
                        typename LocQuadraticImpl::jacobian_type lead_jac,
                        typename LocQuadraticImpl::hessian_type lead_hess )
              : EvalImpl(pot), LeadingLevelOwner<LocQuadraticImpl>( lead_pot, lead_jac, lead_hess ) {}

              
            local_quadratic_evaluation_type evaluate_local_remainder_at_implementation(
              const CVector<D> &arg,
              const CVector<D> &q ) const {
              auto u = LeadingLevelOwner<LocQuadraticImpl>::get_leading_level().evaluate_local_quadratic_at( arg, q );
              return DiagDifference::apply(
                       EvalImpl::evaluate_at( arg ), u );
            }
        };

        
        namespace helper
        {       
          template <class LocQuadImpl, int N, int D>
          struct DiagonalDifference {
            IMPORT_TYPES_FROM( bases::Canonical, N, D );

            struct Inhomogenous {
              static local_quadratic_evaluation_type apply(const potential_evaluation_type& V, const typename LocQuadImpl::local_quadratic_evaluation_type &u ) {
                local_quadratic_evaluation_type C = V.template cast<local_quadratic_return_type>();
                
                for ( int i = 0; i < N; ++i ) {
                  C( i, i ) -= u( i );
                }
                
                return C;
              }
            };
            struct Homogenous {
              static local_quadratic_evaluation_type apply(const potential_evaluation_type& V, const typename LocQuadImpl::local_quadratic_evaluation_type &u ) {
                local_quadratic_evaluation_type C = V.template cast<local_quadratic_return_type>();
                
                for ( int i = 0; i < N; ++i ) {
                  C( i, i ) -= u;
                }
                
                return C;
              }
            };
          };
          
          template <class LocQuadImpl, int D>
          struct DiagonalDifference<LocQuadImpl,1,D> {
            IMPORT_TYPES_FROM( bases::Canonical, 1, D );

            struct Homogenous {
              static local_quadratic_evaluation_type apply(const potential_evaluation_type& V, const typename LocQuadImpl::local_quadratic_evaluation_type &u ) {
                return V - u;
              }
            };
            using Inhomogenous = Homogenous;
          };
        }
        

       

        template<class EvalImpl, class LocQuadraticImpl, int N, int D>
        using Homogenous = General<typename helper::DiagonalDifference<LocQuadraticImpl,N,D>::Homogenous,EvalImpl,LocQuadraticImpl,N,D>;

        template<class EvalImpl, class LocQuadraticImpl, int N, int D>
        using Inhomogenous = General<typename helper::DiagonalDifference<LocQuadraticImpl,N,D>::Inhomogenous,EvalImpl,LocQuadraticImpl,N,D>;
        
      }
    
      template<int N, int D>
      using Homogenous = localRemainder::Homogenous<Evaluation<bases::Canonical,N,D>,LocalQuadratic<bases::Eigen,1,D>,N,D>;

      template<int N, int D>
      using Inhomogenous = localRemainder::Inhomogenous<Evaluation<bases::Canonical,N,D>,LocalQuadratic<bases::Eigen,N,D>,N,D>;
      
    }
  }
}
