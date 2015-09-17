#pragma once
#include "bases.hpp"
#include "matrixPotentials/modules/localRemainder.hpp"
#include "matrixPotentials/modules/evaluation.hpp"
#include "matrixPotentials/modules/exponential.hpp"

namespace lio { 
  namespace matrixPotentials {
    namespace potentials {
      template<int N, int D>
      using Standard = modules::localRemainder::Homogenous<
        modules::Evaluation<bases::Canonical,N,D>,
          modules::Evaluation<bases::Eigen,1,D>,N,D>; 
          
        template<template<int,int> class Basis, int N, int D>
        using Exponential = modules::Exponential<
          modules::Evaluation<Basis,N,D>,Basis,N,D>;

    }
  }
  template<int N, int D>
  using MatrixPotential = matrixPotentials::potentials::Standard<N,D>;
  
  template<template<int,int> class Basis, int N, int D>
  using MatrixPotentialExponential = matrixPotentials::potentials::Exponential<Basis,N,D>;
  
  template<template<int,int> class Basis, int N, int D>
  using MatrixPotentialEvaluate = matrixPotentials::modules::Evaluation<Basis,N,D>;
  
  
}
