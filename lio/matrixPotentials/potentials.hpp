#pragma once
#include "bases.hpp"
#include "matrixPotentials/modules/localRemainder.hpp"
#include "matrixPotentials/modules/evaluation.hpp"
#include "matrixPotentials/modules/exponential.hpp"

namespace lio
{
  namespace matrixPotentials
  {
    namespace potentials
    {
      /**
       * \brief Matrix potential in canonical basis for which homogenous local remainder can be computed
       * 
       * \tparam N
       * Number of levels (dimension of square matrix when evaluated)
       * \tparam D
       * Dimension of argument space
       */
      template <int N, int D>
      using Standard = modules::localRemainder::Homogenous <
                       modules::Evaluation<bases::Canonical, N, D>,
                       modules::Evaluation<bases::Eigen, 1, D>,
                       N,
                       D >;
        
     /**
       * \brief Matrix potential for which the exponential can be computed
       * 
       * \tparam Basis
       * Which basis (bases::Eigen or bases::Canonical) the potential is given in
       * \tparam N
       * Number of levels (dimension of square matrix when evaluated)
       * \tparam D
       * Dimension of argument space
       */                
      template <template <int, int> class Basis, int N, int D>
      using Exponential =
        modules::Exponential<modules::Evaluation<Basis, N, D>, Basis, N, D>;
    }
  }
  // More convenient typedefs
  template <int N, int D>
  using MatrixPotential = matrixPotentials::potentials::Standard<N, D>;
  
  template <template <int, int> class Basis, int N, int D>
  using MatrixPotentialExponential =
    matrixPotentials::potentials::Exponential<Basis, N, D>;
    
  template <template <int, int> class Basis, int N, int D>
  using MatrixPotentialEvaluate =
    matrixPotentials::modules::Evaluation<Basis, N, D>;
}
