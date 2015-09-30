#include <iostream>
#include <fstream>
#include <complex>

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


using namespace waveblocks;
struct Level : public matrixPotentials::modules::taylor::Abstract<Level,CanonicalBasis<1,1>> {
  template <template <typename...> class Tuple = std::tuple>
  Tuple<potential_evaluation_type, jacobian_evaluation_type, hessian_evaluation_type> taylor_at_implementation( const argument_type &x ) const {
        const real_t sigma = 0.038088;
    const real_t a =     0.944858;
               return Tuple<potential_evaluation_type,jacobian_evaluation_type,hessian_evaluation_type>(
       sigma / std::pow(std::cosh(x/a),2),
        - (2.0*sigma*std::tanh(x/a)) / (a*std::pow(std::cosh(x/a),2)),
         (2*sigma*std::cosh(2.0*x/a) - 4*sigma) / (a*a*std::pow(std::cosh(x/a),4)));
  }
};

struct Potential : public matrixPotentials::modules::evaluation::Abstract<Potential,CanonicalBasis<1,1>> {
  complex_t evaluate_at_implementation(const complex_t& x) const {
            const real_t sigma = 0.038088;
    const real_t a =     0.944858;
    return sigma / std::pow(std::cosh(x/a),2);
  }
};

struct Remain : public matrixPotentials::modules::localRemainder::Abstract<Remain, 1,1>, public Potential, public LeadingLevelOwner<Level> {
  complex_t evaluate_local_remainder_at( const complex_t &x,
                    const complex_t &q ) const {
    const real_t sigma = 0.038088;
    const real_t a =     0.944858;
    const auto xmq = x - q;
    
    const  auto J = - (2.0*sigma*std::tanh(q/a)) / (a*std::pow(std::cosh(q/a),2));
    const auto H = (2*sigma*std::cosh(2.0*q/a) - 4*sigma) / (a*a*std::pow(std::cosh(q/a),4));

    return -J*xmq - 0.5*xmq*H*xmq;
  }
};

using MultiIndex = TinyMultiIndex<unsigned short, 1>;

template<class Potential>
real_t potential_energy(const ScalarHaWp<1, MultiIndex>& packet, const Potential& V) {
  using TQR = waveblocks::TensorProductQR <waveblocks::GaussHermiteQR<4>>;
  HomogeneousInnerProduct<1, MultiIndex, TQR> ip;
  return ip.quadrature(packet, [&V] (CMatrix<1,Eigen::Dynamic> nodes, CMatrix<1,1> pos) -> CMatrix<1,Eigen::Dynamic> {
    const dim_t n_nodes = nodes.cols();
    CMatrix<1,Eigen::Dynamic> result(1, n_nodes);
    for(int i = 0; i < n_nodes; ++i)  {
      result(0,i) = V.evaluate_at(pos[0]);
    }
    return result;
  }).real();
}

real_t kinetic_energy(const ScalarHaWp<1, MultiIndex>& packet) {
  HaWpGradientOperator<1,MultiIndex> nabla;
  HaWpGradient<1,MultiIndex> gradwp = nabla(packet);
  complex_t result(0,0);
  for (size_t i = 0 ; i < gradwp.n_components(); ++i) {
    for (const auto& c : gradwp.component(i).coefficients()) {
      result += std::conj(c) * c;
    }
  }
  return result.real();
}



int main() {
    const int N = 1;
    const int D = 1;
    const int K = 512;

    const real_t T = 70;
    const real_t dt = 0.005;

    const real_t eps = 0.1530417681822;


    // The parameter set of the initial wavepacket
    CMatrix<D,D> Q; Q(0,0) = 3.5355339059327;
    CMatrix<D,D> P; P(0,0) = complex_t(0,0.2828427124746);
    RVector<D> q; q[0] = -7.5589045088306;
    RVector<D> p; p[0] = 0.2478854736792;
    complex_t S = 0.0;

    // Setting up the wavepacket
    ShapeEnumerator<D, MultiIndex> enumerator;
    ShapeEnum<D, MultiIndex> shape_enum = enumerator.generate(HyperCubicShape<D>(K));
    HaWpParamSet<D> param_set(q,p,Q,P);
    Coefficients coeffs = Coefficients::Ones(std::pow(K, D), 1);
    ScalarHaWp<D,MultiIndex> packet;

    packet.eps() = eps;
    packet.parameters() = param_set;
    packet.shape() = std::make_shared<ShapeEnum<D,MultiIndex>>(shape_enum);
    packet.coefficients() = coeffs;

    // Defining the potential
    const real_t sigma = 0.038088;
    const real_t a =     0.944858;

    typename CanonicalBasis<N,D>::potential_type potential = [sigma,a](complex_t x) {
        return sigma / std::pow(std::cosh(x/a),2);
    };
    typename ScalarLeadingLevel<D>::jacobian_type leading_jac = [sigma,a](complex_t x) {
        return - (2.0*sigma*std::tanh(x/a)) / (a*std::pow(std::cosh(x/a),2));
    };
    typename ScalarLeadingLevel<D>::hessian_type leading_hess = [sigma,a](complex_t x) {
        return (2*sigma*std::cosh(2.0*x/a) - 4*sigma) / (a*a*std::pow(std::cosh(x/a),4));
    };
    typename ScalarLeadingLevel<D>::potential_type leading_level = potential;

    //~ ScalarMatrixPotential<D> V(potential, leading_level, leading_jac, leading_hess);
    Remain V;
    
    // Quadrature rules
    using TQR = waveblocks::TensorProductQR <waveblocks::GaussHermiteQR<4>>;

    // Defining the propagator
    propagators::Hagedorn<N,D,MultiIndex, TQR> propagator;

    // Preparing the file
    utilities::PacketWriter<ScalarHaWp<D,MultiIndex>> writer("tunneling_1D.out");

    // Propagation
    for (real_t t = 0; t < T; t += dt) {
        std::cout << t << std::endl;

        writer.store_packet(t,packet,S);
        propagator.propagate(packet,dt,V,S);
        std::cout << kinetic_energy(packet) << "," << potential_energy(packet,V); 
        //~ std::cout << packet.parameters() << std::endl;
    }

    writer.store_packet(T,packet,S);
}
