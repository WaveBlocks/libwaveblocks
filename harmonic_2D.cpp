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



using namespace waveblocks;
int main() {
  const int N = 1;
  const int D = 2;

  const real_t T = 12;
  const real_t dt = 0.01;

  const real_t eps = 0.1;

  using MultiIndex = TinyMultiIndex<unsigned short, D>;
  
  // The parameter set of the initial wavepacket
  CMatrix<D,D> Q = CMatrix<D,D>::Identity();
  CMatrix<D,D> P = complex_t(0,1)*CMatrix<D,D>::Identity();
  RVector<D> q = {-3.0,0.0};
  RVector<D> p = {0.0,0.5};
  CVector<N> S; S[0] = 0.;

  // Setting up the wavepacket
  ShapeEnumerator<D, MultiIndex> enumerator;
  ShapeEnum<D, MultiIndex> shape_enum =
    enumerator.generate(HyperCubicShape<D>(N));
  HaWpParamSet<D> param_set(q,p,Q,P);
  std::vector<complex_t> coeffs(std::pow(N, D), 1.0);
  InhomogeneousHaWp<D,MultiIndex> packet(N);
  auto& component = packet.component(0);


  packet.eps() = eps;
  component.parameters() = param_set;
  component.shape() = std::make_shared<ShapeEnum<D,MultiIndex>>(shape_enum);
  component.coefficients() = coeffs;

  // Defining the potential
  typename CanonicalBasis<N,D>::potential_type potential = [](CVector<D> x) {
    return 0.25*(x[0]*x[0] + x[1]*x[1]).real();
  };
  typename InhomogenousLeadingLevel<N,D>::potential_type leading_level = potential;
  typename InhomogenousLeadingLevel<N,D>::jacobian_type leading_jac =
    [D](CVector<D> x) {
      return 0.5*CVector<D>{x[0], x[1]};
    };
  typename InhomogenousLeadingLevel<N,D>::hessian_type leading_hess =
    [D](CVector<D> x) {
      return 0.5*CMatrix<D,D>::Identity();
    };

    
    
  InhomogenousMatrixPotential<N,D> V(potential,leading_level,leading_jac,leading_hess);

  // Quadrature rules
  using TQR = waveblocks::TensorProductQR < waveblocks::GaussHermiteQR<3>,
              waveblocks::GaussHermiteQR<4>>;
  // Defining the propagator
  propagators::Hagedorn<N,D,MultiIndex, TQR> propagator;

  for (auto& coeff: coeffs){
  std::cout <<coeff << " ";}
  std::cout << std::endl;
  
  // Propagation
  for (int t = 0; t < T; t += dt) {
    propagator.propagate(packet,dt,V,S);
  }
  for (auto& coeff: packet.component(0).coefficients()){
  std::cout <<coeff << " ";}
  
  std::cout << std::endl;

  
  

  




}
