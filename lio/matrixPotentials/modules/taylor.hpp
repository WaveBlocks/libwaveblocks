#pragma once
#include "macros.hpp"
#include "matrixPotentials/bases.hpp"
#include "types.hpp"
#include "utilities/evaluations.hpp"
#include "matrixPotentials/modules/evaluation.hpp"
#include "matrixPotentials/modules/jacobian.hpp"
#include "matrixPotentials/modules/hessian.hpp"

#include <tuple>

namespace waveblocks
{
  namespace matrixPotentials
  {
    namespace modules
    {
      namespace taylor
      {
        /**
       * \brief Abstract class for potential evaluation
       * 
       * A matrix potential inheriting an implementation of this module
       * can evaluate its potential, jacobian and hessian in one or multiple points
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
          
          
          template <template <typename...> class Tuple = std::tuple>
          Tuple<potential_evaluation_type, jacobian_evaluation_type, hessian_evaluation_type> taylor_at( const CVector<D> &g ) const {
            static_cast<const Subtype*>(this)->taylor_at_implementation(g);
          }
          
          template < template <typename...> class Tuple = std::tuple,
                   template <typename...> class grid_in = std::vector,
                   template <typename...> class grid_out = grid_in >
          grid_out< Tuple<potential_evaluation_type,jacobian_evaluation_type,hessian_evaluation_type>>taylor( const grid_in<CVector<D> > &args ) const {
            return utilities::evaluate_function_in_grid < CVector<D>,
                   Tuple<>,
                   grid_in,
                   grid_out,
                   function_t > (
                     std::bind( &Self::taylor_at, this, std::placeholders::_1 ), args );
          }
        };
        
        /**
         * \brief Helper class for easier template specialization
         * 
         * This wraps concrete implementations of the Abstract base class
         * 
         * \tparam N
         * Number of levels (dimension of square matrix when evaluated)
         * \tparam D
         * Dimension of argument space
         */
        template <class EvalImpl, class JacImpl, class HessImpl, template<int,int> class Basis, int N, int D>
        struct Standard : public Abstract<Standard<EvalImpl, JacImpl, HessImpl, Basis, N, D>, Basis, N, D>, public EvalImpl, public JacImpl, public HessImpl {
                IMPORT_TYPES_FROM( Basis, N, D );
                
              public:
                Standard(potential_type potential,
                      jacobian_type jacobian,
                      hessian_type hessian) : EvalImpl(potential), JacImpl(jacobian), HessImpl(hessian) {}
                  
              public:
                template <template <typename...> class Tuple = std::tuple>
                Tuple<potential_evaluation_type, jacobian_evaluation_type, hessian_evaluation_type> taylor_at_implementation( const CVector<D> &g ) const {
                             return Tuple<potential_evaluation_type,jacobian_evaluation_type,hessian_evaluation_type>(
                     EvalImpl::evaluate_at( g ), JacImpl::evaluate_jacobian_at( g ), HessImpl::evaluate_hessian_at( g ) );
                };
            };
        };
        
      
      template <template <int, int> class Basis, int N, int D>
      using Taylor = taylor::Standard<Evaluation<Basis,N,D>, Jacobian<Basis,N,D>, Hessian<Basis,N,D>, Basis, N, D>;
    }
  }
}
