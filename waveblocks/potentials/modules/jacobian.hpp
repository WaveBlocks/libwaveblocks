#pragma once

#include "../../types.hpp"
#include "../../utilities/evaluations.hpp"

#include "../bases.hpp"


namespace waveblocks
{
  namespace potentials
  {
    namespace modules
    {
      namespace jacobian
      {
        /**
       * \brief Abstract class for potential evaluation
       *
       * A matrix potential inheriting an implementation of this module
       * can evaluate its jacobian in one or multiple points
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


          jacobian_evaluation_type evaluate_jacobian_at( const argument_type &arg ) const {
            return static_cast<const Subtype*>(this)->evaluate_jacobian_at_implementation( arg );
          }

          template < template <typename...> class grid_in = std::vector,
                   template <typename...> class grid_out = grid_in >
          grid_out<jacobian_evaluation_type> evaluate_jacobian(
            const grid_in<argument_type > &args ) const {
            return utilities::evaluate_function_in_grid < argument_type,
                   jacobian_evaluation_type,
                   grid_in,
                   grid_out,
                   function_t > (
                     std::bind( &Self::evaluate_jacobian_at, this, std::placeholders::_1 ),
                     args );
          }
      };

        /**
         * \brief Helper class for easier template specialization
         *
         * This wraps concrete implementations of the Abstract base class
         *
         * \tparam Basis
         * Which basis (bases::Eigen or bases::Canonical) the potential is given in
         */
        template <class Basis>
        struct Standard : Abstract<Standard<Basis>, Basis> {

                IMPORT_TYPES_FROM(Basis )

              private:
                jacobian_type jacobian;

              public:
                Standard(jacobian_type jacobian)
                  : jacobian( jacobian ) {}

              public:
                jacobian_evaluation_type evaluate_jacobian_at_implementation(
                  const argument_type &arg ) const {
                  return utilities::FunctionMatrixEvaluator < Basis::number_of_levels,
                        Basis::number_of_columns,
                         GMatrix,
                         argument_type,
                         jacobian_return_type,
                         function_t >::apply( jacobian, arg );
                }
            };

      }

      template <class Basis>
      using Jacobian = jacobian::Standard<Basis>;
    }
  }
}
