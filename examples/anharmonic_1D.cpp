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
#include "waveblocks/utilities/energy.hpp"


using namespace waveblocks;
struct Level : public matrixPotentials::modules::taylor::Abstract<Level,CanonicalBasis<1,1>> {
    template <template <typename...> class Tuple = std::tuple>
        Tuple<potential_evaluation_type, jacobian_evaluation_type, hessian_evaluation_type> taylor_at_implementation( const argument_type &x ) const {
        return Tuple<potential_evaluation_type,jacobian_evaluation_type,hessian_evaluation_type>(
                                                                                                 1.0 + std::pow(x,4),
                                                                                                 4.0*std::pow(x,3),
                                                                                                 12.0*x*x);
    }
};

struct Potential : public matrixPotentials::modules::evaluation::Abstract<Potential,CanonicalBasis<1,1>> {
    complex_t evaluate_at_implementation(const complex_t& x) const {
        return 1.0 + std::pow(x,4);
    }
};

struct Remain : public matrixPotentials::modules::localRemainder::Abstract<Remain, 1,1>, public Potential, public LeadingLevelOwner<Level> {
    complex_t evaluate_local_remainder_at( const complex_t &x,
                                           const complex_t &q ) const {
        const auto xmq = x - q;
        const auto J = 4.0*std::pow(q,3);
        const auto H = 12.0*x*x;
        return -J*xmq - 0.5*xmq*H*xmq;
    }
};



int main() {
    const int N = 1;
    const int D = 1;
    const int K = 128;

    const real_t T = 6;
    const real_t dt = 0.01;

    const real_t eps = 0.1;

    using MultiIndex = TinyMultiIndex<unsigned short, D>;

    // The parameter set of the initial wavepacket
    CMatrix<D,D> Q; Q(0,0) = 1.0;
    CMatrix<D,D> P; P(0,0) = complex_t(0,1);
    RVector<D> q; q[0] = 0.0;
    RVector<D> p; p[0] = 1.0;
    complex_t S = 0.0;

    // Setting up the wavepacket
    ShapeEnumerator<D, MultiIndex> enumerator;
    ShapeEnum<D, MultiIndex> shape_enum = enumerator.generate(HyperCubicShape<D>(K));
    HaWpParamSet<D> param_set(q,p,Q,P,S);
    Coefficients coeffs = Coefficients::Zero(std::pow(K, D), 1);
    coeffs[0] = 1.0;
    ScalarHaWp<D,MultiIndex> packet;

    packet.eps() = eps;
    packet.parameters() = param_set;
    packet.shape() = std::make_shared<ShapeEnum<D,MultiIndex>>(shape_enum);
    packet.coefficients() = coeffs;

    Remain V;

    // Quadrature rules
    using TQR = waveblocks::TensorProductQR<waveblocks::GaussHermiteQR<K+4>>;

    // Defining the propagator
    propagators::Hagedorn<N,D,MultiIndex, TQR> propagator;

    // Preparing the file
    utilities::PacketWriter<ScalarHaWp<D,MultiIndex>> writer("anharmonic_1D.hdf5");

    // Propagation
    for (real_t t = 0; t < T; t += dt) {
        propagator.propagate(packet,dt,V);

        std::cout << "----------------------------" << std::endl;
        std::cout << "Time: " << t << std::endl;
        std::cout << packet.parameters() << std::endl;

        real_t kinetic = kinetic_energy<D,MultiIndex>(packet);
        real_t potential = potential_energy<Remain,D,MultiIndex, TQR>(packet,V);
        writer.store_energies(t,potential,kinetic);
        std::cout << "............................" << std::endl;
        writer.store_packet(t,packet);
    }
}
