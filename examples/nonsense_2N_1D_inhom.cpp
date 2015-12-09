#include "propagators/Hagedorn.hpp"
#include "matrixPotentials/potentials.hpp"
#include "matrixPotentials/bases.hpp"
#include "types.hpp"
#include "hawp_commons.hpp"
#include "tiny_multi_index.hpp"
#include "shape_enumerator.hpp"
#include "shape_hypercubic.hpp"
#include "hawp_paramset.hpp"
#include <iostream>
#include <fstream>


using namespace waveblocks;
int main() {
  const int N = 2;
  const int D = 1;
  const int K = 3;
  const real_t sigma_x = 0.5;
  const real_t sigma_y = 0.5;

  const real_t T = 12;
  const real_t dt = 0.01;

  const real_t eps = 0.1;

  using MultiIndex = TinyMultiIndex<unsigned short, D>;

  // The parameter set of the initial wavepacket
  CMatrix<D,D> Q = CMatrix<D,D>::Identity();
  CMatrix<D,D> P = complex_t(0,1)*CMatrix<D,D>::Identity();
  RVector<D> q;
  RVector<D> p;
  CVector<N> S;

  // Setting up the wavepacket
  ShapeEnumerator<D, MultiIndex> enumerator;
  ShapeEnum<D, MultiIndex> shape_enum =
    enumerator.generate(HyperCubicShape<D>(K));
  HaWpParamSet<D> param_set(q,p,Q,P);
  Coefficients coeffs = Coefficients::Ones(std::pow(K, D), 1);
  InhomogeneousHaWp<D,MultiIndex> packet(N);


  packet.eps() = eps;
  packet.component(0).parameters() = param_set;

  packet.component(0).shape() = std::make_shared<ShapeEnum<D,MultiIndex>>(shape_enum);
  packet.component(0).coefficients() = coeffs;

  packet.component(1).parameters() = param_set;
  packet.component(1).shape() = std::make_shared<ShapeEnum<D,MultiIndex>>(shape_enum);
  packet.component(1).coefficients() = coeffs;

  // Defining the potential
  typename CanonicalBasis<N,D>::potential_type potential;

  potential(0,0) = [sigma_x,sigma_y,N](complex_t x) {
    return 0.5*(sigma_x*x).real();
  };

  potential(0,1) = [sigma_x,sigma_y,N](complex_t x) {
    return 0.5*(sigma_x*x).real();
  };
  potential(1,0) = [sigma_x,sigma_y,N](complex_t x) {
    return 0.5*(sigma_x*x).real();
  };
  potential(1,1) = [sigma_x,sigma_y,N](complex_t x) {
    return 0.5*(sigma_x*x).real();
  };

  typename InhomogenousLeadingLevel<N,D>::potential_type leading_level;
  leading_level(0) = [sigma_x,sigma_y,N](complex_t x) {
    return 0.5*(sigma_x*x).real();
  };


  typename InhomogenousLeadingLevel<N,D>::jacobian_type leading_jac;
  leading_jac(0) = [D,sigma_x,sigma_y,N](complex_t x) {
    return 0.5*(sigma_x*x).real();
  };

  typename InhomogenousLeadingLevel<N,D>::hessian_type leading_hess;
  leading_hess(0) =
    [D,sigma_x,sigma_y,N](complex_t x) {
    return 0.5*(sigma_x*x).real();
    };

    leading_level(1) = leading_level(0);
    leading_jac(1) = leading_jac(0);
    leading_hess(1) = leading_hess(0);

  InhomogenousMatrixPotential<N,D> V(potential,leading_level,leading_jac,leading_hess);

  // Quadrature rules
  using TQR = waveblocks::innerproducts::TensorProductQR<waveblocks::innerproducts::GaussHermiteQR<3>>;

  // Defining the propagator
  propagators::Hagedorn<N,D,MultiIndex, TQR> propagator;

  // Preparing the file

  // Propagation
  for (real_t t = 0; t < T; t += dt) {
    propagator.propagate(packet,dt,V,S);
  }
}
