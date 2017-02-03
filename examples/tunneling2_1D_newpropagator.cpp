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
        Tuple<potential_evaluation_type, jacobian_evaluation_type, hessian_evaluation_type> taylor_at_implementation(const argument_type &x ) const {
        const real_t sigma = 0.038088;
        const real_t a =     0.944858;
        const double d =     4.0;
        const complex_t cxmd = std::pow(std::cosh((x-d)/a),2);
        const complex_t cxpd = std::pow(std::cosh((x+d)/a),2);

        return Tuple<potential_evaluation_type,jacobian_evaluation_type,hessian_evaluation_type>(
          sigma / cxmd +
          sigma / cxpd,
          -(2*sigma * std::tanh((x-d)/a) / cxmd)
          -(2*sigma * std::tanh((x+d)/a) / cxpd),
          (2*sigma*std::cosh(2.0*(x-d)/a) - 4*sigma) / (a*a*cxmd*cxmd) +
          (2*sigma*std::cosh(2.0*(x+d)/a) - 4*sigma) / (a*a*cxpd*cxpd));
    }
};

struct Potential : public potentials::modules::evaluation::Abstract<Potential,CanonicalBasis<1,1>> {
    complex_t evaluate_at_implementation(const complex_t& x) const {
        const real_t sigma = 0.038088;
        const real_t a =     0.944858;
        const double d =     4.0;
        const complex_t cxmd = std::pow(std::cosh((x-d)/a),2);
        const complex_t cxpd = std::pow(std::cosh((x+d)/a),2);
        return
            sigma / cxmd +
            sigma / cxpd;
    }
};

struct Remain : public potentials::modules::localRemainder::Abstract<Remain, 1,1>, public Potential, public LeadingLevelOwner<Level> {
    complex_t evaluate_local_remainder_at(const complex_t &x,
                                          const complex_t &q ) const {
        const real_t sigma = 0.038088;
        const real_t a =     0.944858;
        const double d = 4.0;
        const complex_t cxmd = std::pow(std::cosh((x-d)/a),2);
        const complex_t cxpd = std::pow(std::cosh((x+d)/a),2);
        const complex_t cqmd = std::pow(std::cosh((q-d)/a),2);
        const complex_t cqpd = std::pow(std::cosh((q+d)/a),2);
        const auto xmq = x - q;
        const auto V =
            sigma / cxmd +
            sigma / cxpd;
        const auto U =
            sigma / cqmd +
            sigma / cqpd;
        const auto J =
            -(2*sigma * std::tanh((q-d)/a) / cqmd)
            -(2*sigma * std::tanh((q+d)/a) / cqpd);
        const auto H =
            (2*sigma*std::cosh(2.0*(q-d)/a) - 4*sigma) / (a*a*cqmd*cqmd) +
            (2*sigma*std::cosh(2.0*(q+d)/a) - 4*sigma) / (a*a*cqpd*cqpd);

        return V - U - J*xmq - 0.5*xmq*H*xmq;
    }
};


int main() {

	// TODO: implement with new propagators

    const int N = 1;
    const int D = 1;
    const int K = 700;

    const real_t T = 100;
    const real_t dt = 0.005;

    const real_t eps = 0.1530417681822;

    using MultiIndex = wavepackets::shapes::TinyMultiIndex<unsigned short, D>;

    // The parameter set of the initial wavepacket
    CMatrix<D,D> Q; Q(0,0) = 3.5355339059327;
    CMatrix<D,D> P; P(0,0) = complex_t(0,0.2828427124746);
    RVector<D> q; q[0] = -7.5589045088306 - 4.0;
    RVector<D> p; p[0] = 0.2478854736792;

    // Setting up the wavepacket
    wavepackets::shapes::ShapeEnumerator<D, MultiIndex> enumerator;
    wavepackets::shapes::ShapeEnum<D, MultiIndex> shape_enum = enumerator.generate(wavepackets::shapes::HyperCubicShape<D>(K));
    wavepackets::HaWpParamSet<D> param_set(q,p,Q,P,0.0);
    Coefficients coeffs = Coefficients::Zero(std::pow(K, D), 1);
    coeffs[0] = 1.0;
    wavepackets::ScalarHaWp<D,MultiIndex> packet;

    packet.eps() = eps;
    packet.parameters() = param_set;
    packet.shape() = std::make_shared<wavepackets::shapes::ShapeEnum<D,MultiIndex>>(shape_enum);
    packet.coefficients() = coeffs;

    Remain V;

    // Quadrature rules
    using TQR = innerproducts::TensorProductQR<innerproducts::GaussHermiteQR<K+4>>;

    // Defining the propagator
    propagators::Hagedorn<N,D,MultiIndex, TQR> propagator;

    // Preparing the file
    io::hdf5writer<D> mywriter("tunneling2_1D_cpp.hdf5");
    mywriter.set_write_energies(true);
    mywriter.prestructuring<MultiIndex>(packet,dt);

    //time 0
    std::cout << "Time: " << 0 << std::endl;
    std::cout << packet.parameters() << std::endl;

    mywriter.store_packet(packet);
    real_t ekin = observables::kinetic_energy<D,MultiIndex>(packet);
    real_t epot = observables::potential_energy<Remain,D,MultiIndex, TQR>(packet,V);
    mywriter.store_energies(epot,ekin);
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
        real_t epot = observables::potential_energy<Remain,D,MultiIndex, TQR>(packet,V);
        std::cout << "E: (p,k,t) " << epot << ", " << ekin << ", " << epot+ekin << std::endl;
        mywriter.store_energies(epot,ekin);
    }
    mywriter.poststructuring();
}
