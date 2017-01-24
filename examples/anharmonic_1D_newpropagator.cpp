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
#include "waveblocks/propagators/HagedornPropagator.hpp"
#include "waveblocks/propagators/SemiclassicalPropagator.hpp"
#include "waveblocks/propagators/MagnusPropagator.hpp"
#include "waveblocks/propagators/Pre764Propagator.hpp"
#include "waveblocks/propagators/McL42Propagator.hpp"
#include "waveblocks/propagators/McL84Propagator.hpp"
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

namespace split = propagators::splitting_parameters;

int main() {

    // General parameters
    const int N = 1;
    const int D = 1;
    const int K = 128;

    const real_t T = 1;
    const real_t Dt = 0.01;

    const real_t eps = 0.1;

    using MultiIndex = wavepackets::shapes::TinyMultiIndex<unsigned short, D>;


	/////////////////////////////////////////////////////
	// Building the WavePacket
	/////////////////////////////////////////////////////

    wavepackets::ScalarHaWp<D,MultiIndex> packet;
    packet.eps() = eps;

    // Initial Parameter Set
	wavepackets::HaWpParamSet<D> param_set;
	param_set.p(RVector<D>::Ones()); // set p - give momentum
    packet.parameters() = param_set;

    // Basis Shapes
    wavepackets::shapes::ShapeEnumerator<D, MultiIndex> enumerator;
    wavepackets::shapes::ShapeEnum<D, MultiIndex> shape_enum = enumerator.generate(wavepackets::shapes::HyperCubicShape<D>(K));
    packet.shape() = std::make_shared<wavepackets::shapes::ShapeEnum<D,MultiIndex>>(shape_enum);

    // Set Coefficients: Gaussian Wavepacket phi_0 with c_0 = 1
    Coefficients coefs = Coefficients::Zero(std::pow(K, D), 1);
    coefs[0] = 1.0;
    packet.coefficients() = coefs;


	/////////////////////////////////////////////////////
	// Defining the Propagator
	/////////////////////////////////////////////////////

    Remain V;                                               // Potential
    using QR = innerproducts::GaussHermiteQR<K+4>;          // Quadrature rule

	// TODO: consider having Potential as an argument to propagate
    // propagators::Hagedorn<N,D,MultiIndex,QR> propagator;


	// Define propagator

	using Packet_t = wavepackets::ScalarHaWp<D,MultiIndex>;
	propagators::HagedornPropagator<N,D,MultiIndex,QR,Remain,Packet_t> pHagedorn(packet,V);
	propagators::SemiclassicalPropagator<N,D,MultiIndex,QR,Remain,Packet_t,propagators::SplitCoefs<1,1>> pSemiclassical(packet,V,split::coefLT);
	propagators::MagnusPropagator<N,D,MultiIndex,QR,Remain,Packet_t,propagators::SplitCoefs<1,1>> pMagnus(packet,V,split::coefLT);
	propagators::Pre764Propagator<N,D,MultiIndex,QR,Remain,Packet_t,propagators::SplitCoefs<1,1>> pPre764(packet,V,split::coefLT);
	propagators::McL42Propagator<N,D,MultiIndex,QR,Remain,Packet_t,propagators::SplitCoefs<1,1>> pMcL42(packet,V,split::coefLT);
	propagators::McL84Propagator<N,D,MultiIndex,QR,Remain,Packet_t,propagators::SplitCoefs<1,1>> pMcL84(packet,V,split::coefLT);


	////////////////////////////////////////////////////
	// Callback Function - writing to file
	////////////////////////////////////////////////////
	// set up writer
	// Preparing the file and I/O writer
	io::hdf5writer<D> mywriter("data.hdf5");
	mywriter.set_write_norm(true);
	mywriter.set_write_energies(true);
	//////////////////////////////////////////////////////////////////////////////
	
	std::function<void(unsigned,real_t)> writeenergies = [&](unsigned i, real_t t) {
		// (void) i; // avoid unused variable warning
		// (void) t; // avoid unused variable warning 

		std::cout << "i: " << i << "  t: " << t << "\n";
		// Write Data
		// TODO: provide getekin, getepot functions here in propagator
		real_t ekin = observables::kinetic_energy<D,MultiIndex>(packet);
		real_t epot = observables::potential_energy<Remain,D,MultiIndex,QR>(packet,V);
		// TODO: pre/post processing requires applying transformations before measuring

		mywriter.store_packet(packet);
		mywriter.store_norm(packet);
		mywriter.store_energies(epot,ekin);
	};
	
	std::function<void(unsigned,real_t)> emptycallback = [&](unsigned,real_t) {};

	//////////////////////////////////////////////////////////////////////////////

	mywriter.prestructuring<MultiIndex>(packet,Dt);
	pHagedorn.evolve(T,Dt,writeenergies);
	pSemiclassical.evolve(T,Dt,emptycallback);
	pMagnus.evolve(T,Dt,emptycallback);
	pPre764.evolve(T,Dt,emptycallback);
	pMcL42.evolve(T,Dt,emptycallback);
	pMcL84.evolve(T,Dt,emptycallback);
	mywriter.poststructuring();

	return 0;

}
