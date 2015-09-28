#pragma once
#include "macros.hpp"
#include "matrixPotentials/bases.hpp"
#include "types.hpp"
#include "utilities/evaluations.hpp"

namespace waveblocks
{
  namespace matrixPotentials
  {
    namespace modules
    {
      namespace evaluation
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
        template <class Subtype, class Basis>
        struct Abstract {
          using Self = Abstract<Subtype, Basis>;
          IMPORT_TYPES_FROM( Basis);
          
          potential_evaluation_type evaluate_at( const argument_type &arg ) const {

            return static_cast<const Subtype*>(this)->evaluate_at_implementation( arg );
          }
          
          template < template <typename...> class grid_in = std::vector,
                   template <typename...> class grid_out = grid_in >
          grid_out<potential_evaluation_type> evaluate(
            const grid_in<argument_type > &args ) const {
            return utilities::evaluate_function_in_grid < argument_type,
                   potential_evaluation_type,
                   grid_in,
                   grid_out,
                   function_t > (
                     std::bind( &Self::evaluate_at, this, std::placeholders::_1 ), args );
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
        template <class Basis>
        struct Standard : Abstract<Standard<Basis>, Basis> {
          IMPORT_TYPES_FROM( Basis);

          private:
          potential_type potential;

          public:
          Standard( potential_type potential)
            : potential( potential ){}
            
          public:
          potential_evaluation_type evaluate_at_implementation(
            const argument_type &arg ) const {
            return utilities::FunctionMatrixEvaluator < Basis::number_of_levels, Basis::number_of_columns,
                   GMatrix,
                   argument_type,
                   potential_return_type,
                   function_t >::apply( potential, arg );
          }
        };
      }
      
      template <class Basis>
      using Evaluation = evaluation::Standard<Basis>;
    }
  }
}
