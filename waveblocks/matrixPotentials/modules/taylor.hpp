#pragma once

#include <tuple>

#include "../../types.hpp"
#include "../../utilities/evaluations.hpp"
#include "../bases.hpp"
#include "evaluation.hpp"
#include "jacobian.hpp"
#include "hessian.hpp"


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
       */
        template <class Subtype, class Basis>
        struct Abstract {
          using Self = Abstract<Subtype, Basis>;
          IMPORT_TYPES_FROM( Basis)


          template <template <typename...> class Tuple = std::tuple>
          Tuple<potential_evaluation_type, jacobian_evaluation_type, hessian_evaluation_type> taylor_at( const argument_type &g ) const {
            return static_cast<const Subtype*>(this)->taylor_at_implementation(g);
          }

          template < template <typename...> class Tuple = std::tuple,
                   template <typename...> class grid_in = std::vector,
                   template <typename...> class grid_out = grid_in >
          grid_out< Tuple<potential_evaluation_type,jacobian_evaluation_type,hessian_evaluation_type>>taylor( const grid_in<argument_type > &args ) const {
            return utilities::evaluate_function_in_grid < argument_type,
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
         * \tparam EvalImpl Implemntation of the evaluation module.
         * \tparam JacImpl Implementation of the jacobian module.
         * \tparam HessImpl Implementation of the hessian module.
         * \tparam Basis
         * Which basis (bases::Eigen or bases::Canonical) the potential is given in
         */
        template <class EvalImpl, class JacImpl, class HessImpl, class Basis>
        struct Standard : public Abstract<Standard<EvalImpl, JacImpl, HessImpl, Basis>, Basis>, public EvalImpl, public JacImpl, public HessImpl {
                IMPORT_TYPES_FROM( Basis)

              public:
                Standard(potential_type potential,
                      jacobian_type jacobian,
                      hessian_type hessian) : EvalImpl(potential), JacImpl(jacobian), HessImpl(hessian) {}

              public:
                template <template <typename...> class Tuple = std::tuple>
                Tuple<potential_evaluation_type, jacobian_evaluation_type, hessian_evaluation_type> taylor_at_implementation( const argument_type &g ) const {
                             return Tuple<potential_evaluation_type,jacobian_evaluation_type,hessian_evaluation_type>(
                     EvalImpl::evaluate_at( g ), JacImpl::evaluate_jacobian_at( g ), HessImpl::evaluate_hessian_at( g ) );
                }
            };
        }


      template <class Basis>
      using Taylor = taylor::Standard<Evaluation<Basis>, Jacobian<Basis>, Hessian<Basis>, Basis>;
    }
  }
}
