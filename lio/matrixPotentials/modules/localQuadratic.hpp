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
        template <class Subtype, class Basis>
        struct Abstract {
            using Self = Abstract<Subtype, Basis>;
            IMPORT_TYPES_FROM( Basis)
            
            
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
        
        template <class TaylorImpl, class Basis>
        class Standard : public Abstract<Standard<TaylorImpl, Basis>, Basis>,
        public TaylorImpl
        {
        public:
          IMPORT_TYPES_FROM( Basis)

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

            for ( int l = 0; l < Basis::number_of_levels; ++l ) {
              for ( int m = 0; m < Basis::number_of_columns; ++m )  {
                const auto& V = V_mat(l,m);
                const auto& J = J_mat(l,m);
                const auto& H = H_mat(l,m);

                auto result = V;

                for ( int i = 0; i < Basis::argument_dimension; ++i ) {
                  auto xmqi = x[i] - q[i];
                  result += J[i] * ( xmqi );

                  for ( int j = 0; j < Basis::argument_dimension; ++j ) {
                    result += 0.5 * xmqi * H( i, j ) * ( x[j] - q[j] );
                  }
                }

                result_matrix(l,m) = result;
              }
            }
              
            return result_matrix;
          }
        };

        template <class TaylorImpl, template <int, int, int> class B, int N, int C>
        class Standard<TaylorImpl, B<N, 1, C>> : public Abstract<Standard<TaylorImpl, B<N,1, C>>, B<N,1, C>>,
        public TaylorImpl
        {
        public:
          using Basis = B<N,1, C>;
          IMPORT_TYPES_FROM( Basis)

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

            for ( int l = 0; l < Basis::number_of_levels; ++l ) {
              for ( int m = 0; m < Basis::number_of_columns; ++m )  {
                const auto& V = V_mat(l,m);
                const auto& J = J_mat(l,m);
                const auto& H = H_mat(l,m);


                result_matrix(l,m) = V + J*xmq + 0.5*xmq*H*xmq;
              }
            }
              
            return result_matrix;
          }
        };

        
        
        template <class TaylorImpl, template <int, int, int> class B, int D, int C>
        class Standard<TaylorImpl, B<1, D,C>> : public Abstract <
          Standard<TaylorImpl, B<1,D,C>>, B<1,D,C> >,        public TaylorImpl
        {
          public:
            using Basis = B<1,D,C>;
            IMPORT_TYPES_FROM( Basis)

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

        template <class TaylorImpl, template <int, int, int> class B, int C>
        class Standard<TaylorImpl, B<1,1,C>> : public Abstract<Standard<TaylorImpl, B<1, 1,C>>, B<1, 1,C>>,
        public TaylorImpl
        {
        public:
          using Basis = B<1,1,C>;
          IMPORT_TYPES_FROM( Basis)

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
      template <class Basis>
      using LocalQuadratic = localQuadratic::Standard<Taylor<Basis>, Basis>;
    }
  }
}
