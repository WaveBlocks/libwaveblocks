#pragma once

#include "bases.hpp"
#include "modules/localRemainder.hpp"
#include "modules/evaluation.hpp"
#include "modules/exponential.hpp"


namespace waveblocks
{
  namespace potentials
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

  using potentials::potentials::InhomogenousMatrixPotential;
  using potentials::potentials::HomogenousMatrixPotential;
  using potentials::potentials::ScalarMatrixPotential;


  template<int N, int D>
  using InhomogenousLeadingLevel = potentials::modules::LocalQuadratic<potentials::bases::Eigen<N,D>>;

  template<int N, int D>
  using HomogenousLeadingLevel = potentials::modules::LocalQuadratic<potentials::bases::Eigen<1,D>>;

  template<int D>
  using ScalarLeadingLevel = HomogenousLeadingLevel<1,D>;
}
