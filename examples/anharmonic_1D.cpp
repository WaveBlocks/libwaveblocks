#include <iostream>
#include <fstream>
#include <complex>

#include "waveblocks/types.hpp"
#include "waveblocks/wavepackets/shapes/tiny_multi_index.hpp"
#include "waveblocks/potentials/potentials.hpp"
#include "waveblocks/potentials/bases.hpp"
#include "waveblocks/wavepackets/hawp_paramset.hpp"
#include "waveblocks/wavepackets/hawp_commons.hpp"
#include "waveblocks/wavepackets/shapes/shape_enumerator.hpp"
#include "waveblocks/wavepackets/shapes/shape_hypercubic.hpp"
#include "waveblocks/innerproducts/gauss_hermite_qr.hpp"
#include "waveblocks/innerproducts/tensor_product_qr.hpp"
#include "waveblocks/propagators/Hagedorn.hpp"
#include "waveblocks/observables/energy.hpp"
#include "waveblocks/io/hdf5writer.hpp"


using namespace waveblocks;

struct Level : public potentials::modules::taylor::Abstract<Level,CanonicalBasis<1,1>> {
    template <template <typename...> class Tuple = std::tuple>
        Tuple<potential_evaluation_type, jacobian_evaluation_type, hessian_evaluation_type> taylor_at_implementation( const argument_type &x ) const {
        return Tuple<potential_evaluation_type,jacobian_evaluation_type,hessian_evaluation_type>(1.0 + std::pow(x,4),
                                                                                                 4.0*std::pow(x,3),
                                                                                                 12.0*x*x);
    }
};

struct Potential : public potentials::modules::evaluation::Abstract<Potential,CanonicalBasis<1,1>> {
    complex_t evaluate_at_implementation(const complex_t& x) const {
        return 1.0 + std::pow(x,4);
    }
};

struct Remain : public potentials::modules::localRemainder::Abstract<Remain, 1,1>, public Potential, public LeadingLevelOwner<Level> {
    complex_t evaluate_local_remainder_at( const complex_t &x,
                                           const complex_t &q ) const {
        const auto xmq = x - q;
        const auto V = 1.0 + std::pow(x,4);
        const auto U = 1.0 + std::pow(q,4);
        const auto J = 4.0*std::pow(q,3);
        const auto H = 12.0*q*q;
        return V - U - J*xmq - 0.5*xmq*H*xmq;
    }
};


int main() {
    const int N = 1;
    const int D = 1;
    const int K = 128;

    const real_t T = 6;
    const real_t dt = 0.01;

    const real_t eps = 0.1;

    using MultiIndex = wavepackets::shapes::TinyMultiIndex<unsigned short, D>;

    // The parameter set of the initial wavepacket
    CMatrix<D,D> Q; Q(0,0) = 1.0;
    CMatrix<D,D> P; P(0,0) = complex_t(0,1);
    RVector<D> q; q[0] = 0.0;
    RVector<D> p; p[0] = 1.0;
    complex_t S = 0.0;

    // Setting up the wavepacket
    wavepackets::shapes::ShapeEnumerator<D, MultiIndex> enumerator;
    wavepackets::shapes::ShapeEnum<D, MultiIndex> shape_enum = enumerator.generate(wavepackets::shapes::HyperCubicShape<D>(K));
    wavepackets::HaWpParamSet<D> param_set(q,p,Q,P,S);
    Coefficients coeffs = Coefficients::Zero(std::pow(K, D), 1);
    coeffs[0] = 1.0;
    wavepackets::ScalarHaWp<D,MultiIndex> packet;

    packet.eps() = eps;
    packet.parameters() = param_set;
    packet.shape() = std::make_shared<wavepackets::shapes::ShapeEnum<D,MultiIndex>>(shape_enum);
    packet.coefficients() = coeffs;

    Remain V;

    // Quadrature rules
    using QR = innerproducts::GaussHermiteQR<K+4>;

    // Defining the propagator
    propagators::Hagedorn<N,D,MultiIndex,QR> propagator;

    // Preparing the file
    io::hdf5writer<D> mywriter("anharmonic_1D_cpp.hdf5");
    mywriter.set_write_energies(true);
    mywriter.set_write_norm(true);
    mywriter.prestructuring<MultiIndex>(packet,dt);

    //time=0
    std::cout << "Time: " << 0 << std::endl;
    mywriter.store_packet(packet);
    std::cout << packet.parameters() << std::endl;
    real_t ekin = observables::kinetic_energy<D,MultiIndex>(packet);
    real_t epot = observables::potential_energy<Remain,D,MultiIndex,QR>(packet,V);
    mywriter.store_energies(epot,ekin);
    mywriter.store_norm(packet);
    std::cout << "E: (p,k,t) " << epot << ", " << ekin << ", " << epot+ekin << std::endl;


    // Propagation
    for (real_t t = 0; t < T; t += dt) {
        std::cout << "Time: " << t+dt << std::endl;

        // Propagate
        propagator.propagate(packet,dt,V);
        std::cout << packet.parameters() << std::endl;
        mywriter.store_packet(packet);

        // Compute energies
        real_t ekin = observables::kinetic_energy<D,MultiIndex>(packet);
        real_t epot = observables::potential_energy<Remain,D,MultiIndex,QR>(packet,V);
        real_t etot = ekin + epot;
        std::cout << "E: (p,k,t) " << epot << ", " << ekin << ", " << etot << std::endl;
        mywriter.store_energies(epot,ekin);
        mywriter.store_norm(packet);
    }
    mywriter.poststructuring();
}
