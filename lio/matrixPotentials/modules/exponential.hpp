#pragma once
#include "macros.hpp"
#include "matrixPotentials/bases.hpp"
#include "utilities/evaluations.hpp"
#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

namespace lio
{
  namespace matrixPotentials
  {
    namespace modules
    {
      namespace exponential
      {
           /**
         * \brief Abstract class for exponential of potential evaluation
         * 
         * A matrix potential inheriting an implementation of this module
         * can evaluate the exponential of its potential
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
        class Abstract
        {
            using Self = Abstract<Subtype, Basis, N, D>;
            IMPORT_TYPES_FROM( Basis, N, D );
            
            potential_evaluation_type evaluate_exponential_at( const RVector<D> &arg,
                const real_t &factor = 1 ) {
              that.evaluate_exponential_at_implementation( arg, factor );
            }
            
            template < template <typename...> class grid_in = std::vector,
                     template <typename...> class grid_out = grid_in >
            grid_out<potential_evaluation_type> evaluate_exponential(
              grid_in<RVector<D> > args,
              real_t factor = 1 ) {
              return utilities::evaluate_function_in_grid < RVector<D>,
                     potential_evaluation_type,
                     grid_in,
                     grid_out,
                     function_t > (
                       std::bind( &Self::evaluate_exponential_at,
                                  this,
                                  std::placeholders::_1,
                                  factor ),
                       args );
            }
        };
        
         /**
         * \brief Implementation of exponential of potential evaluation
         * 
         * A matrix potential inheriting an implementation of this module
         * can evaluate its exponential potential and also inherits the evaluation module
         * given by EvalImpl
         * 
         * \tparam EvalImpl The implementation of the evaluation module used
         * \tparam Basis
         * Which basis (bases::Eigen or bases::Canonical) the potential is given in
         * \tparam N
         * Number of levels (dimension of square matrix when evaluated)
         * \tparam D
         * Dimension of argument space
         */
        template <class EvalImpl, template <int, int> class Basis, int N, int D>
        struct Standard : public Abstract<Standard<EvalImpl, Basis, N, D>, Basis, N, D>,
          public EvalImpl {
          IMPORT_TYPES_FROM( Basis, N, D );
          
          potential_evaluation_type evaluate_exponential_at_implementation(
            const RVector<D> &arg,
            real_t factor ) {
            // Compute matrix
            auto values = ( arg );
            potential_evaluation_type result;
            
            // Compute exponential
            Eigen::MatrixExponential<potential_evaluation_type> m_exp( factor * values );
            
            m_exp.compute( result );
            return result;
          }
        };
      }
      
      template <class EvalImpl, template <int, int> class Basis, int N, int D>
      using Exponential = exponential::Standard<EvalImpl, Basis, N, D>;
    }
  }
}
