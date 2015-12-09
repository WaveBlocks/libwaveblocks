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
#include "waveblocks/innerproducts/gauss_hermite_qr.hpp"
#include "waveblocks/innerproducts/tensor_product_qr.hpp"
#include "waveblocks/observables/energy.hpp"
#include "waveblocks/utilities/packetWriter.hpp"


using namespace waveblocks;
struct Level : public matrixPotentials::modules::taylor::Abstract<Level,CanonicalBasis<1,1>> {
    template <template <typename...> class Tuple = std::tuple>
        Tuple<potential_evaluation_type, jacobian_evaluation_type, hessian_evaluation_type> taylor_at_implementation(const argument_type &x ) const {
        const real_t sigma = 0.038088;
        const real_t a =     0.944858;
        const complex_t cx = std::pow(std::cosh(x/a),2);
        return Tuple<potential_evaluation_type,jacobian_evaluation_type,hessian_evaluation_type>(sigma / cx,
                                                                                                 - (2.0*sigma*std::tanh(x/a)) / (a*cx),
                                                                                                 (2*sigma*std::cosh(2.0*x/a) - 4*sigma) / (a*a*cx*cx));
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
    complex_t evaluate_local_remainder_at(const complex_t &x,
                                          const complex_t &q ) const {
        const real_t sigma = 0.038088;
        const real_t a =     0.944858;
        const complex_t cq = std::pow(std::cosh(q/a),2);
        const auto xmq = x - q;
        const auto V = sigma / std::pow(std::cosh(x/a),2);
        const auto U = sigma / cq;
        const auto J = - (2.0*sigma*std::tanh(q/a)) / (a*cq);
        const auto H = (2*sigma*std::cosh(2.0*q/a) - 4*sigma) / (a*a*cq*cq);
        return V - U - J*xmq - 0.5*xmq*H*xmq;
    }
};


int main() {
    const int N = 1;
    const int D = 1;
    const int K = 512;

    const real_t T = 70;
    const real_t dt = 0.005;

    const real_t eps = 0.1530417681822;

    using MultiIndex = TinyMultiIndex<unsigned short, D>;

    // The parameter set of the initial wavepacket
    CMatrix<D,D> Q; Q(0,0) = 3.5355339059327;
    CMatrix<D,D> P; P(0,0) = complex_t(0,0.2828427124746);
    RVector<D> q; q[0] = -7.5589045088306;
    RVector<D> p; p[0] = 0.2478854736792;

    // Setting up the wavepacket
    ShapeEnumerator<D, MultiIndex> enumerator;
    ShapeEnum<D, MultiIndex> shape_enum = enumerator.generate(HyperCubicShape<D>(K));
    HaWpParamSet<D> param_set(q,p,Q,P,0.0);
    Coefficients coeffs = Coefficients::Zero(std::pow(K, D), 1);
    coeffs[0] = 1.0;
    ScalarHaWp<D,MultiIndex> packet;

    packet.eps() = eps;
    packet.parameters() = param_set;
    packet.shape() = std::make_shared<ShapeEnum<D,MultiIndex>>(shape_enum);
    packet.coefficients() = coeffs;

    Remain V;

    // Quadrature rules
    using TQR = waveblocks::innerproducts::TensorProductQR<waveblocks::innerproducts::GaussHermiteQR<K+4>>;

    // Defining the propagator
    propagators::Hagedorn<N,D,MultiIndex, TQR> propagator;

    // Preparing the file
    utilities::PacketWriter<ScalarHaWp<D,MultiIndex>> writer("tunneling_1D.hdf5");

    // Propagation
    for (real_t t = 0; t < T; t += dt) {
        std::cout << "Time: " << t << std::endl;

        // Propagate
        propagator.propagate(packet,dt,V);
        std::cout << packet.parameters() << std::endl;
        writer.store_packet(t,packet);

        // Compute energies
        real_t kinetic = kinetic_energy<D,MultiIndex>(packet);
        real_t potential = potential_energy<Remain,D,MultiIndex, TQR>(packet,V);
        real_t total = kinetic+potential;
        std::cout << "E: (p,k,t) " << potential << ", " << kinetic << ", " << total << std::endl;
        std::cout << potential << "," << kinetic << ", "<< total << std::endl;
        writer.store_energies(t,potential,kinetic);
    }
}
