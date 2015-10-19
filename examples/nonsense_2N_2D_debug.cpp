#include <iostream>
#include <fstream>

#include "waveblocks/propagators/Hagedorn.hpp"
#include "waveblocks/matrixPotentials/potentials.hpp"
#include "waveblocks/matrixPotentials/bases.hpp"
#include "waveblocks/types.hpp"
#include "waveblocks/hawp_commons.hpp"
#include "waveblocks/tiny_multi_index.hpp"
#include "waveblocks/shape_enumerator.hpp"
#include "waveblocks/shape_hypercubic.hpp"
#include "waveblocks/hawp_paramset.hpp"


using namespace waveblocks;
int main() {
  const int N = 2;
  const int D = 2;
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
  RVector<D> q = {-3.0,0.0};
  RVector<D> p = {0.0,0.5};

  // Setting up the wavepacket
  ShapeEnumerator<D, MultiIndex> enumerator;
  ShapeEnum<D, MultiIndex> shape_enum =
    enumerator.generate(HyperCubicShape<D>(K));
  HaWpParamSet<D> param_set(q,p,Q,P,0.0);
  Coefficients coeffs = Coefficients::Ones(std::pow(K, D), 1);
  HomogeneousHaWp<D,MultiIndex> packet(N);


  packet.eps() = eps;
  packet.parameters() = param_set;

  packet.component(0).shape() = std::make_shared<ShapeEnum<D,MultiIndex>>(shape_enum);
  packet.component(0).coefficients() = coeffs;

  packet.component(1).shape() = std::make_shared<ShapeEnum<D,MultiIndex>>(shape_enum);
  packet.component(1).coefficients() = coeffs;

  // Defining the potential
  typename CanonicalBasis<N,D>::potential_type potential;

  potential(0,0) = [sigma_x,sigma_y,N](CVector<D> x) {
    return 0.5*(sigma_x*x[0]*x[0] + sigma_y*x[1]*x[1]).real();
  };

  potential(0,1) = [sigma_x,sigma_y,N](CVector<D> x) {
    return 0.5*(sigma_x*x[0]*x[0] + sigma_y*x[1]*x[1]).real();
  };
  potential(1,0) = [sigma_x,sigma_y,N](CVector<D> x) {
    return 0.5*(sigma_x*x[0]*x[0] + sigma_y*x[1]*x[1]).real();
  };
  potential(1,1) = [sigma_x,sigma_y,N](CVector<D> x) {
    return 0.5*(sigma_x*x[0]*x[0] + sigma_y*x[1]*x[1]).real();
  };
  
  typename HomogenousLeadingLevel<N,D>::potential_type leading_level;
  leading_level = [sigma_x,sigma_y,N](CVector<D> x) {
    return 0.5*(sigma_x*x[0]*x[0] + sigma_y*x[1]*x[1]).real();
  };

 
  typename HomogenousLeadingLevel<N,D>::jacobian_type leading_jac;
  leading_jac = [sigma_x,sigma_y,N](CVector<D> x) {
    return CVector<D>{sigma_x*x[0], sigma_y*x[1]};
  };

  typename HomogenousLeadingLevel<N,D>::hessian_type leading_hess;
  leading_hess =
      [sigma_x,sigma_y,N](CVector<D> /*x*/) {
      CMatrix<D,D> res;
      res(0,0) = sigma_x;
      res(1,1) = sigma_y;
      return res;
    };
    
  HomogenousMatrixPotential<N,D> V(potential,leading_level,leading_jac,leading_hess);

  // Quadrature rules
  using TQR = waveblocks::TensorProductQR < waveblocks::GaussHermiteQR<3>,
              waveblocks::GaussHermiteQR<4>>;
  // Defining the propagator
  propagators::Hagedorn<N,D,MultiIndex, TQR> propagator;


  // Preparing the file
  
  
  // Output before
  std::cout << "Before:\n-------\n";
  for (size_t i = 0; i < packet.n_components(); ++i) {
    std::cout << "Component " << i << ":\n";
    std::cout << packet.component(i).parameters();
    std::cout << packet.component(i).coefficients() << "\n";
  }
  
  // Propagation
  for (real_t t = 0; t < T; t += dt) {
    propagator.propagate(packet,dt,V);
  }
  
  // Output after
  std::cout << "\nAfter:\n-------\n";
  for (size_t i = 0; i < packet.n_components(); ++i) {
    std::cout << "Component " << i << ":\n";
    std::cout << packet.component(i).parameters();
    std::cout << packet.component(i).coefficients() << "\n";
  }
}
