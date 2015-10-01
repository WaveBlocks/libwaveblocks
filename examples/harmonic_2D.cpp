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
#include "waveblocks/hawp_gradient_operator.hpp"
#include "waveblocks/hawp_gradient_evaluator.hpp"
#include "waveblocks/utilities/energy.hpp"

using namespace waveblocks;
int main() {
  const int N = 1;
  const int D = 2;
  const int K = 5;

  const real_t sigma_x = 0.5;
  const real_t sigma_y = 0.5;

  const real_t tol = 1e-10;

  const real_t T = 12;
  const real_t dt = 0.01;

  const real_t eps = 0.1;

  using MultiIndex = TinyMultiIndex<unsigned long, D>;

  // Parameter set PI
  CMatrix<D,D> Q = CMatrix<D,D>::Identity();
  CMatrix<D,D> P = complex_t(0,1) * CMatrix<D,D>::Identity();
  RVector<D> q = {-3.0, 0.0};
  RVector<D> p = { 0.0, 0.5};
  complex_t S = 0.;
  HaWpParamSet<D> param_set(q,p,Q,P,S);

  // Basis shape
  ShapeEnumerator<D, MultiIndex> enumerator;
  ShapeEnum<D, MultiIndex> shape_enum = enumerator.generate(HyperCubicShape<D>(K));

  // Gaussian Wavepacket phi_00 with c_00 = 1
  Coefficients coeffs = Coefficients::Ones(std::pow(K, D), 1);
  coeffs[0] = 1.0;
  Coefficients coefforig = Coefficients(coeffs);

  // Assemble packet
  ScalarHaWp<D,MultiIndex> packet;
  packet.eps() = eps;
  packet.parameters() = param_set;
  packet.shape() = std::make_shared<ShapeEnum<D,MultiIndex>>(shape_enum);
  packet.coefficients() = coeffs;

  // Defining the potential
  typename CanonicalBasis<N,D>::potential_type potential =
      [sigma_x,sigma_y](CVector<D> x) {
      return 0.5*(sigma_x*x[0]*x[0] + sigma_y*x[1]*x[1]).real();
  };
  typename ScalarLeadingLevel<D>::potential_type leading_level = potential;
  typename ScalarLeadingLevel<D>::jacobian_type leading_jac =
      [sigma_x,sigma_y](CVector<D> x) {
      return CVector<D>{sigma_x*x[0], sigma_y*x[1]};
  };
  typename ScalarLeadingLevel<D>::hessian_type leading_hess =
      [sigma_x,sigma_y](CVector<D>) {
      CMatrix<D,D> res;
      res(0,0) = sigma_x;
      res(1,1) = sigma_y;
      return res;
  };

  ScalarMatrixPotential<D> V(potential,leading_level,leading_jac,leading_hess);

  // Quadrature rules
  using TQR = waveblocks::TensorProductQR<
      waveblocks::GaussHermiteQR<5>,
      waveblocks::GaussHermiteQR<5>>;

  // Defining the propagator
  propagators::Hagedorn<N,D,MultiIndex, TQR> propagator;

  // Preparing the file
  utilities::PacketWriter<ScalarHaWp<D,MultiIndex>> writer("harmonic_2D.hdf5");

  // Propagation
  for (real_t t = 0; t < T; t += dt) {
    std::cout << "Time: " << t << std::endl;

    propagator.propagate(packet,dt,V);
    writer.store_packet(t,packet);

    real_t ekin = kinetic_energy<D,MultiIndex>(packet);
    real_t epot = potential_energy<ScalarMatrixPotential<D>,D,MultiIndex, TQR>(packet,V);
    writer.store_energies(t,epot,ekin);

    // Assure constant coefficients
    auto diff = (packet.coefficients() - coefforig).array().abs();
    auto norm = diff.matrix().template lpNorm<Eigen::Infinity>();
    bool flag = norm > tol ? false : true;

    std::cout << "coefficients constant? " << (flag ? "yes" : "no") << std::endl;
    std::cout << packet.parameters() << std::endl;
  }
}
