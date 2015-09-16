#pragma once
#include "bases.hpp"
#include "modules/localRemainder.hpp"
#include "modules/evaluation.hpp"

namespace matrixPotentials {
  namespace potentials {
    template<int N, int D>
    using Standard = modules::localRemainder::Homogenous<
      modules::Evaluation<bases::Canonical,N,D>,
        modules::Evaluation<bases::Eigen,1,D>,N,D>;    
  }
  
  template<int N, int D>
  using MatrixPotential = potentials::Standard<N,D>;
}


