#pragma once

#include "bases.hpp"
#include "modules/localRemainder.hpp"
#include "modules/evaluation.hpp"
#include "modules/exponential.hpp"

namespace waveblocks
{
  namespace matrixPotentials
  {
    namespace potentials
    {

      template<int N, int D>
      using InhomogenousMatrixPotential =
        modules::Inhomogenous<N,D>;

      template<int N, int D>
      using HomogenousMatrixPotential =
        modules::Homogenous<N,D>;

      template<int D>
      using ScalarMatrixPotential = HomogenousMatrixPotential<1,D>;
    }
  }

  using matrixPotentials::potentials::InhomogenousMatrixPotential;
  using matrixPotentials::potentials::HomogenousMatrixPotential;
  using matrixPotentials::potentials::ScalarMatrixPotential;


  template<int N, int D>
  using InhomogenousLeadingLevel = matrixPotentials::modules::LocalQuadratic<matrixPotentials::bases::Eigen<N,D>>;

  template<int N, int D>
  using HomogenousLeadingLevel = matrixPotentials::modules::LocalQuadratic<matrixPotentials::bases::Eigen<1,D>>;

  template<int D>
  using ScalarLeadingLevel = HomogenousLeadingLevel<1,D>;


}
