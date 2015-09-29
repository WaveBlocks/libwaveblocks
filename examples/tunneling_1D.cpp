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
#include "waveblocks/utilities/packetWriter.hpp"


using namespace waveblocks;
int main() {
  const int N = 1;
  const int D = 1;
  const int K = 12;
  
  const real_t sigma = 0.05;
  const real_t T = 70;
  const real_t dt = 0.005;

  const real_t eps = 0.1530417681822;

  using MultiIndex = TinyMultiIndex<unsigned short, D>;

  // The parameter set of the initial wavepacket
  CMatrix<D,D> Q; Q(0,0) = 3.5355339059327;
  CMatrix<D,D> P; P(0,0) = complex_t(0,0.2828427124746);
  RVector<D> q; q[0] = -7.5589045088306;
  RVector<D> p; p[0] = 0.2478854736792;
  complex_t S = 0.;

  // Setting up the wavepacket
  ShapeEnumerator<D, MultiIndex> enumerator;
  ShapeEnum<D, MultiIndex> shape_enum =
    enumerator.generate(HyperCubicShape<D>(K));
  HaWpParamSet<D> param_set(q,p,Q,P);
  std::vector<complex_t> coeffs(std::pow(K, D), 1.);
  ScalarHaWp<D,MultiIndex> packet;


  packet.eps() = eps;
  packet.parameters() = param_set;
  packet.shape() = std::make_shared<ShapeEnum<D,MultiIndex>>(shape_enum);
  packet.coefficients() = coeffs;

  // Defining the potential
  typename CanonicalBasis<N,D>::potential_type potential = [sigma](complex_t x) {
    return 0.25 * sigma * std::pow(x,4);
  };
  typename ScalarLeadingLevel<D>::potential_type leading_level = potential;
  typename ScalarLeadingLevel<D>::jacobian_type leading_jac =
    [sigma](complex_t x) {
      return sigma*x*x*x;
    };
  typename ScalarLeadingLevel<D>::hessian_type leading_hess =
    [sigma](complex_t x) {
      return 3*sigma*x*x;
    };
    
    
  ScalarMatrixPotential<D> V(potential,leading_level,leading_jac,leading_hess);

  // Quadrature rules
  using TQR = waveblocks::TensorProductQR <waveblocks::GaussHermiteQR<4>>;
  // Defining the propagator
  propagators::Hagedorn<N,D,MultiIndex, TQR> propagator;


  // Preparing the file
  utilities::PacketWriter<ScalarHaWp<D,MultiIndex>> writer("tunneling_1D.out");

  // Propagation
  for (real_t t = 0; t < T; t += dt) {
    writer.store_packet(t,packet,S);
    propagator.propagate(packet,dt,V,S);
  }

  writer.store_packet(T,packet,S);
  
  
  

  




}
