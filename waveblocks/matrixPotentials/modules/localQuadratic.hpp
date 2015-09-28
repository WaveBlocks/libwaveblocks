#pragma once
#include "macros.hpp"
#include "types.hpp"
#include "utilities/evaluations.hpp"
#include "matrixPotentials/modules/taylor.hpp"

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
            
            
          public:
            potential_evaluation_type evaluate_local_quadratic_at(
              const argument_type &arg,
              const argument_type &position ) const {
              return static_cast<const Subtype*>(this)->evaluate_local_quadratic_at_implementation( arg,position );
            }
            
            template < template <typename...> class grid_in = std::vector,
                     template <typename...> class grid_out = grid_in >
            grid_out<potential_evaluation_type> evaluate_local_remainder(
              const grid_in<argument_type > &args,
              argument_type position ) const {
              return utilities::evaluate_function_in_grid < argument_type,
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
        
        template <class TaylorImpl, template <int, int> class Basis, int N, int D>
        class Standard : public Abstract<Standard<TaylorImpl, Basis, N, D>, Basis, N, D>,
        public TaylorImpl
        {
        public:
          IMPORT_TYPES_FROM( Basis, N, D );

        public:
          Standard( potential_type potential,
                    jacobian_type jacobian,
                    hessian_type hessian )
            : TaylorImpl( potential, jacobian, hessian ) {
          }
          
          potential_evaluation_type evaluate_local_quadratic_at_implementation(
            const argument_type &x,
            const argument_type &q ) const {

            potential_evaluation_type result_matrix;
            
            auto V_mat = TaylorImpl::evaluate_at(q );
            auto J_mat = TaylorImpl::evaluate_jacobian_at(q );
            auto H_mat = TaylorImpl::evaluate_hessian_at(q );

            for ( int l = 0; l < N; ++l ) {
              for ( int m = 0; m < N; ++m )  {
                const auto& V = V_mat(l,m);
                const auto& J = J_mat(l,m);
                const auto& H = H_mat(l,m);

                auto result = V;

                for ( int i = 0; i < D; ++i ) {
                  auto xmqi = x[i] - q[i];
                  result += J[i] * ( xmqi );

                  for ( int j = 0; j < D; ++j ) {
                    result += 0.5 * xmqi * H( i, j ) * ( x[j] - q[j] );
                  }
                }

                result_matrix(l,m) = result;
              }
            }
              
            return result_matrix;
          }
        };

        template <class TaylorImpl, template <int, int> class Basis, int N>
        class Standard<TaylorImpl, Basis, N, 1> : public Abstract<Standard<TaylorImpl, Basis, N, 1>, Basis, N, 1>,
        public TaylorImpl
        {
        public:
          IMPORT_TYPES_FROM( Basis, N, 1);

        public:
          Standard( potential_type potential,
                    jacobian_type jacobian,
                    hessian_type hessian )
            : TaylorImpl( potential, jacobian, hessian ) {
          }
          
          potential_evaluation_type evaluate_local_quadratic_at_implementation(
            const argument_type &x,
            const argument_type &q ) const {
            auto xmq = x - q;

            potential_evaluation_type result_matrix;
            
            auto V_mat = TaylorImpl::evaluate_at(q );
            auto J_mat = TaylorImpl::evaluate_jacobian_at(q );
            auto H_mat = TaylorImpl::evaluate_hessian_at(q );

            for ( int l = 0; l < N; ++l ) {
              for ( int m = 0; m < N; ++m )  {
                const auto& V = V_mat(l,m);
                const auto& J = J_mat(l,m);
                const auto& H = H_mat(l,m);


                result_matrix(l,m) = V + J*xmq + 0.5*xmq*H*xmq;
              }
            }
              
            return result_matrix;
          }
        };

        
        
        template <class TaylorImpl, template <int, int> class Basis, int D>
        class Standard<TaylorImpl, Basis, 1, D> : public Abstract <
          Standard<TaylorImpl, Basis, 1, D>,
          Basis,
          1,
          D > ,
        public TaylorImpl
        {
          public:
            IMPORT_TYPES_FROM( Basis, 1, D );

            Standard( potential_type potential,
                      jacobian_type jacobian,
                      hessian_type hessian )
              : TaylorImpl( potential, jacobian, hessian ) {
            }
            
          potential_evaluation_type evaluate_local_quadratic_at_implementation(
            const argument_type &x,
            const argument_type &q ) const {
              auto V = TaylorImpl::evaluate_at( q );
              auto J = TaylorImpl::evaluate_jacobian_at(q );
              auto H = TaylorImpl::evaluate_hessian_at(q );
           
              auto result = V;
              
              for ( int i = 0; i < D; ++i ) {
                auto xmqi = x[i] - q[i];
                result += J[i] * ( xmqi );
                
                for ( int j = 0; j < D; ++j ) {
                  result += 0.5 * xmqi * H( i, j ) * ( x[j] - q[j] );
                }
              }
              
              return result;
          }
        };

        template <class TaylorImpl, template <int, int> class Basis>
        class Standard<TaylorImpl, Basis, 1, 1> : public Abstract<Standard<TaylorImpl, Basis, 1, 1>, Basis, 1, 1>,
        public TaylorImpl
        {
        public:
          IMPORT_TYPES_FROM( Basis, 1, 1);

        public:
          Standard( potential_type potential,
                    jacobian_type jacobian,
                    hessian_type hessian )
            : TaylorImpl( potential, jacobian, hessian ) {
          }
          
          potential_evaluation_type evaluate_local_quadratic_at_implementation(
            const argument_type &x,
            const argument_type &q ) const {
            auto xmq = x - q;

            potential_evaluation_type result_matrix;
            
            auto V = TaylorImpl::evaluate_at(q );
            auto J = TaylorImpl::evaluate_jacobian_at(q );
            auto H = TaylorImpl::evaluate_hessian_at(q );

            return V + J*xmq + 0.5*xmq*H*xmq;
          }

        };
      }
      template <template <int, int> class Basis, int N, int D>
      using LocalQuadratic = localQuadratic::Standard<Taylor<Basis,N,D>, Basis, N, D>;
    }
  }
}
