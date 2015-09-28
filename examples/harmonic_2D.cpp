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
  const int N = 1;
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
  complex_t S = 0.;

  // Setting up the wavepacket
  ShapeEnumerator<D, MultiIndex> enumerator;
  ShapeEnum<D, MultiIndex> shape_enum =
    enumerator.generate(HyperCubicShape<D>(K));
  HaWpParamSet<D> param_set(q,p,Q,P);
  std::vector<complex_t> coeffs(std::pow(K, D), 1.0);
  ScalarHaWp<D,MultiIndex> packet;


  packet.eps() = eps;
  packet.parameters() = param_set;
  packet.shape() = std::make_shared<ShapeEnum<D,MultiIndex>>(shape_enum);
  packet.coefficients() = coeffs;

  // Defining the potential
  typename CanonicalBasis<N,D>::potential_type potential = [sigma_x,sigma_y](CVector<D> x) {
    return 0.5*(sigma_x*x[0]*x[0] + sigma_y*x[1]*x[1]).real();
  };
  typename ScalarLeadingLevel<D>::potential_type leading_level = potential;
  typename ScalarLeadingLevel<D>::jacobian_type leading_jac =
    [sigma_x,sigma_y](CVector<D> x) {
      return CVector<D>{sigma_x*x[0], sigma_y*x[1]};;
    };
  typename ScalarLeadingLevel<D>::hessian_type leading_hess =
    [sigma_x,sigma_y](CVector<D> x) {
      CMatrix<D,D> res;
      res(0,0) = sigma_x;
      res(1,1) = sigma_y;
      return res;
    };
    
    
  ScalarMatrixPotential<D> V(potential,leading_level,leading_jac,leading_hess);

  // Quadrature rules
  using TQR = waveblocks::TensorProductQR < waveblocks::GaussHermiteQR<3>,
              waveblocks::GaussHermiteQR<4>>;
  // Defining the propagator
  propagators::Hagedorn<N,D,MultiIndex, TQR> propagator;


  // Preparing the file
  std::ofstream csv;
  csv.open ("harmonic_2D.out");
  csv << "t, p, q, P, Q, S";

  
  csv << 0 << "," << param_set.p << "," << param_set.q << "," << param_set.P << "," << param_set.Q << "," << S << std::endl;
  
  // Propagation
  for (real_t t = 0; t < T; t += dt) {
    propagator.propagate(packet,dt,V,S);
    const auto& params = packet.parameters();
    csv << t << "," << params.p << "," <<
     params.q << "," << params.P << "," <<
     params.Q << "," << S  << std::endl;
  }
  
  

  csv.close();


  
  

  




}
