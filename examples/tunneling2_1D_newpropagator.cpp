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
#include "waveblocks/propagators/HagedornPropagator.hpp"
#include "waveblocks/propagators/SemiclassicalPropagator.hpp"
#include "waveblocks/propagators/MG4Propagator.hpp"
#include "waveblocks/propagators/Pre764Propagator.hpp"
#include "waveblocks/propagators/McL42Propagator.hpp"
#include "waveblocks/propagators/McL84Propagator.hpp"
#include "waveblocks/observables/energy.hpp"
#include "waveblocks/io/hdf5writer.hpp"


using namespace waveblocks;

namespace split = propagators::splitting_parameters;

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

	////////////////////////////////////////////////////
    // General parameters
	////////////////////////////////////////////////////

    const int N = 1;
    const int D = 1;
    const int K = 700;

    const real_t T = 100;
    const real_t Dt = 0.005;

    const real_t eps = 0.1530417681822;


	////////////////////////////////////////////////////
    // Composed Types
	////////////////////////////////////////////////////

    using MultiIndex = wavepackets::shapes::TinyMultiIndex<unsigned short, D>;
    using TQR = innerproducts::TensorProductQR<innerproducts::GaussHermiteQR<K+4>>;
    using Packet_t = wavepackets::ScalarHaWp<D,MultiIndex>;
    using Potential_t = Remain;


	/////////////////////////////////////////////////////
	// Building the WavePacket
	/////////////////////////////////////////////////////

    // basis shapes
    wavepackets::shapes::ShapeEnumerator<D, MultiIndex> enumerator;
    wavepackets::shapes::ShapeEnum<D, MultiIndex> shape_enum = enumerator.generate(wavepackets::shapes::HyperCubicShape<D>(K));

    // initial parameters
    CMatrix<D,D> Q; Q(0,0) = 3.5355339059327;
    CMatrix<D,D> P; P(0,0) = complex_t(0,0.2828427124746);
    RVector<D> q; q[0] = -7.5589045088306 - 4.0;
    RVector<D> p; p[0] = 0.2478854736792;
    wavepackets::HaWpParamSet<D> param_set(q,p,Q,P,0.0);

	// initial coefficients
    Coefficients coeffs = Coefficients::Zero(std::pow(K, D), 1);
    coeffs[0] = 1.0;

    // assemble packet
    Packet_t packet;
    packet.eps() = eps;
    packet.parameters() = param_set;
    packet.shape() = std::make_shared<wavepackets::shapes::ShapeEnum<D,MultiIndex>>(shape_enum);
    packet.coefficients() = coeffs;


	/////////////////////////////////////////////////////
    // Defining the potential
	/////////////////////////////////////////////////////

    Potential_t V;


	////////////////////////////////////////////////////
	// Defining the Propagator
	////////////////////////////////////////////////////

	propagators::HagedornPropagator<N,D,MultiIndex,TQR,Potential_t,Packet_t> pHagedorn(packet,V);
	propagators::SemiclassicalPropagator<N,D,MultiIndex,TQR,Potential_t,Packet_t,propagators::SplitCoefs<1,1>> pSemiclassical(packet,V,split::coefLT);
	propagators::MG4Propagator<N,D,MultiIndex,TQR,Potential_t,Packet_t,propagators::SplitCoefs<1,1>> pMG4(packet,V,split::coefLT);
	propagators::Pre764Propagator<N,D,MultiIndex,TQR,Potential_t,Packet_t,propagators::SplitCoefs<1,1>> pPre764(packet,V,split::coefLT);
	propagators::McL42Propagator<N,D,MultiIndex,TQR,Potential_t,Packet_t,propagators::SplitCoefs<1,1>> pMcL42(packet,V,split::coefLT);
	propagators::McL84Propagator<N,D,MultiIndex,TQR,Potential_t,Packet_t,propagators::SplitCoefs<1,1>> pMcL84(packet,V,split::coefLT);


	////////////////////////////////////////////////////
	// Defining Callback Function
	////////////////////////////////////////////////////

	// set up writer: preparing the file and I/O writer
    io::hdf5writer<D> mywriter("tunneling2_1D_newpropagator.hdf5");
    mywriter.set_write_energies(true);

	// Write Data Callback
	std::function<void(unsigned,real_t)> writeenergies = [&](unsigned,real_t) {
		real_t ekin = observables::kinetic_energy<D,MultiIndex>(packet);
		real_t epot = observables::potential_energy<Remain,D,MultiIndex, TQR>(packet,V);
		mywriter.store_energies(epot,ekin);
		mywriter.store_packet(packet);
	};


	////////////////////////////////////////////////////
	// Propagate
	////////////////////////////////////////////////////

    mywriter.prestructuring<MultiIndex>(packet,Dt);

	pHagedorn.evolve(T,Dt,writeenergies);
	pSemiclassical.evolve(T,Dt);
	pMG4.evolve(T,Dt);
	pPre764.evolve(T,Dt);
	pMcL42.evolve(T,Dt);
	pMcL84.evolve(T,Dt);

    mywriter.poststructuring();

    return 0;

}
