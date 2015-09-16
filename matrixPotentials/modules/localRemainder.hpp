#pragma once
#include "../macros.hpp"
#include "../bases.hpp"
#include "localQuadratic.hpp"
#include "../types.hpp"
#include "../utilities/evaluations.hpp"

namespace matrixPotentials {
  namespace modules {
    namespace localRemainder{
      namespace helper {
        
        template<int N>
        struct DiagonalDifference {
          struct Inhomogenous {
            static CMatrix<N,N> apply(RMatrix<N,N> V, const CVector<N>& u) {
              CMatrix<N,N> C = V.template cast<complex_t>();
              for (int i = 0; i < N; ++i) {
                    C(i,i) -= u(i);
              }
              return C;
            }
          };
          struct Homogenous {
            static CMatrix<N,N> apply(RMatrix<N,N> V, const complex_t& u) {
              CMatrix<N,N> C = V.template cast<complex_t>();
              for (int i = 0; i < N; ++i) {
                    C(i,i) -= u;
              }
              return C;
            }
          };
        };
        
        template<>
        struct DiagonalDifference<1> {
          struct Homogenous {
            static complex_t apply(real_t V, const complex_t& u) {
           return V - u;
        }
          };
          using Inhomogenous = Homogenous;
        };
      }
      
      template<class Subtype,
               int N,
               int D >
      struct Abstract {
        CMatrix<N,N> evaluate_local_remainder_at(const RVector<D>& arg, const RVector<D>& q) {
          return that.evaluate_local_remainder_at_implementation(arg,q);
        }
      };


     

      template <class EvalImpl, class LeadingEvalImpl, int N, int D>
      class Homogenous : public Abstract<Homogenous<EvalImpl, LeadingEvalImpl, N,D>,N,D>, EvalImpl{
        IMPORT_TYPES_FROM(bases::Canonical,N,D);

        using leading_level_type = bases::Eigen<1,D>;
        
        LocalQuadratic<LeadingEvalImpl,bases::Eigen,1,D> quadratic;
        
      public:
        Homogenous(
          potential_type pot,
          jacobian_type jac,
          hessian_type hess,
          typename leading_level_type::potential_type lead_pot,
          typename leading_level_type::jacobian_type lead_jac,
          typename leading_level_type::hessian_type lead_hess
        ) : EvalImpl(pot,jac,hess), quadratic(lead_pot,lead_jac,lead_hess) {} 
        
          CMatrix<N,N> evaluate_local_remainder_at_implementation(const RVector<D>& arg, const RVector<D>& q) {
            auto u = quadratic.evaluate_local_quadratic_at(arg,q);
            return helper::DiagonalDifference<N>::Homogenous::apply(EvalImpl::evaluate_at(arg),u);
          }  
      };

      template <class EvalImpl, class LeadingEvalImpl, int N, int D>
      class Inhomogenous : public Abstract<Inhomogenous<EvalImpl, LeadingEvalImpl, N,D>,N,D>, EvalImpl{
        IMPORT_TYPES_FROM(bases::Canonical,N,D);
        
        using leading_level_type = bases::Eigen<N,D>;

        LocalQuadratic<LeadingEvalImpl,bases::Eigen,N,D> quadratic;
      public:
        Inhomogenous(
          potential_type pot,
          jacobian_type jac,
          hessian_type hess,
          typename leading_level_type::potential_type lead_pot,
          typename leading_level_type::jacobian_type lead_jac,
          typename leading_level_type::hessian_type lead_hess
        ) : EvalImpl(pot,jac,hess), quadratic(lead_pot,lead_jac,lead_hess) {} 
        
        
          CMatrix<N,N> evaluate_local_remainder_at_implementation(const RVector<D>& arg, const RVector<D>& q) {
            auto u = quadratic.evaluate_local_quadratic_at(arg,q);
            return helper::DiagonalDifference<N>::Inhomogenous::apply(EvalImpl::evaluate_at(arg),u);
          }  
      };
    }
  }
}
