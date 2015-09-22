#pragma once
#include "bases.hpp"
#include "matrixPotentials/modules/localRemainder.hpp"
#include "matrixPotentials/modules/evaluation.hpp"
#include "matrixPotentials/modules/exponential.hpp"

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

      template<template<int,int> class Basis, int N, int D>
      using MatrixPotential = modules::Taylor<Basis,N,D>;
    }
  }

  template<int N, int D>
  using InhomogenousMatrixPotential =
    matrixPotentials::potentials::InhomogenousMatrixPotential<N,D>;

  template<int N, int D>
  using HomogenousMatrixPotential =
    matrixPotentials::potentials::HomogenousMatrixPotential<N,D>;

  template<template<int,int> class Basis, int N, int D>
    using MatrixPotential = matrixPotentials::potentials::MatrixPotential<Basis,N,D>;
}
