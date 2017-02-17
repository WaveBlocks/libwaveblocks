#include <iostream>
#include <fstream>
#include <functional>

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
using namespace std::placeholders;

int main() {

	////////////////////////////////////////////////////
    // General parameters
	////////////////////////////////////////////////////

    const int N = 1;
    const int D = 1;
    const int K = 16;

    const real_t C = 0.004164;
    const real_t a = 0.896696;
    const real_t x0 = 5.542567;

    const real_t T = 50;
    const real_t Dt = 0.1;

    const real_t eps = 0.0484;


	////////////////////////////////////////////////////
    // Composed Types
	////////////////////////////////////////////////////

    using MultiIndex = wavepackets::shapes::TinyMultiIndex<unsigned short, D>;
    using QR = innerproducts::GaussHermiteQR<K+4>;
	using Packet_t = wavepackets::ScalarHaWp<D,MultiIndex>;
	using Potential_t = ScalarMatrixPotential<D>;


	/////////////////////////////////////////////////////
	// Building the WavePacket
	/////////////////////////////////////////////////////

    // basis shapes
    wavepackets::shapes::ShapeEnumerator<D, MultiIndex> enumerator;
    wavepackets::shapes::ShapeEnum<D, MultiIndex> shape_enum = enumerator.generate(wavepackets::shapes::HyperCubicShape<D>(K));

    // initial parameters
    CMatrix<D,D> Q = complex_t(3.4957,0.) * CMatrix<D,D>::Identity();
    CMatrix<D,D> P = complex_t(0.,0.2861) * CMatrix<D,D>::Identity();
    RVector<D> q = 6.0426 * RVector<D>::Ones();
    RVector<D> p = -0.1100 * RVector<D>::Ones();
    complex_t S = 0.0;
    wavepackets::HaWpParamSet<D> param_set(q,p,Q,P,S);

	// initial coefficients
    // Gaussian wave packet phi_00 with c_00 = 1
    Coefficients coeffs = Coefficients::Zero(std::pow(K,D),1);
    coeffs[0] = 1.0;

	// assemble packet
    Packet_t packet;
    packet.eps() = eps;
    packet.shape() = std::make_shared<wavepackets::shapes::ShapeEnum<D,MultiIndex>>(shape_enum);
    packet.parameters() = param_set;
    packet.coefficients() = coeffs;


	/////////////////////////////////////////////////////
    // Defining the potential
	/////////////////////////////////////////////////////

    typename CanonicalBasis<N,D>::potential_type potential = [C,a,x0](complex_t x) {
		return C * (std::exp(-2.*a*(x-x0)) - 2.*std::exp(-a*(x-x0)));
    };
    typename ScalarLeadingLevel<D>::potential_type leading_level = potential;
    typename ScalarLeadingLevel<D>::jacobian_type leading_jac = [C,a,x0](complex_t x) {
		return C * (-2.*a*std::exp(-2.*a*(x-x0)) + 2.*a*std::exp(-a*(x-x0)));
    };
    typename ScalarLeadingLevel<D>::hessian_type leading_hess = [C,a,x0](complex_t x) {
		return C * (4.*a*a*std::exp(-2.*a*(x-x0)) - 2.*a*a*std::exp(-a*(x-x0)));
    };

    Potential_t V(potential,leading_level,leading_jac,leading_hess);


	////////////////////////////////////////////////////
	// Defining the Propagator
	////////////////////////////////////////////////////

	using SplitCoefs_t = propagators::SplitCoefs<4,4>;
	const SplitCoefs_t splitCoefs = propagators::splitting_parameters::coefY4;

	propagators::HagedornPropagator<N,D,MultiIndex,QR,Potential_t,Packet_t> pHagedorn(packet,V);
	propagators::SemiclassicalPropagator<N,D,MultiIndex,QR,Potential_t,Packet_t,SplitCoefs_t> pSemiclassical(packet,V,splitCoefs);
	propagators::MG4Propagator<N,D,MultiIndex,QR,Potential_t,Packet_t,SplitCoefs_t> pMG4(packet,V,splitCoefs);
	propagators::Pre764Propagator<N,D,MultiIndex,QR,Potential_t,Packet_t,SplitCoefs_t> pPre764(packet,V,splitCoefs);
	propagators::McL42Propagator<N,D,MultiIndex,QR,Potential_t,Packet_t,SplitCoefs_t> pMcL42(packet,V,splitCoefs);
	propagators::McL84Propagator<N,D,MultiIndex,QR,Potential_t,Packet_t,SplitCoefs_t> pMcL84(packet,V,splitCoefs);


	////////////////////////////////////////////////////
	// Defining Callback Function
	////////////////////////////////////////////////////

	auto write = [&](io::hdf5writer<D>& writer,unsigned,real_t) {
		real_t ekin = observables::kinetic_energy<D,MultiIndex>(packet);
		real_t epot = observables::potential_energy<Potential_t,D,MultiIndex,QR>(packet,V);
		writer.store_packet(packet);
		writer.store_norm(packet);
		writer.store_energies(epot,ekin);
	};
	

	////////////////////////////////////////////////////
	// Propagate
	////////////////////////////////////////////////////

	// Hagedorn
	{
    io::hdf5writer<D> writerHagedorn("morse_1D_Hagedorn.hdf5"); writerHagedorn.set_write_norm(true); writerHagedorn.set_write_energies(true);
    writerHagedorn.prestructuring<MultiIndex>(packet,Dt);
	pHagedorn.evolve(T,Dt,std::bind(write,writerHagedorn,_1,_2));
    writerHagedorn.poststructuring();
	}

	// Semiclassical
	{
    io::hdf5writer<D> writerSemiclassical("morse_1D_Semiclassical.hdf5"); writerSemiclassical.set_write_norm(true); writerSemiclassical.set_write_energies(true);
    writerSemiclassical.prestructuring<MultiIndex>(packet,Dt);
	pSemiclassical.evolve(T,Dt,std::bind(write,writerSemiclassical,_1,_2));
    writerSemiclassical.poststructuring();
	}

	// MG4
	{
    io::hdf5writer<D> writerMG4("morse_1D_MG4.hdf5"); writerMG4.set_write_norm(true); writerMG4.set_write_energies(true);
    writerMG4.prestructuring<MultiIndex>(packet,Dt);
	pMG4.evolve(T,Dt,std::bind(write,writerMG4,_1,_2));
    writerMG4.poststructuring();
	}

	// McL42
	{
    io::hdf5writer<D> writerMcL42("morse_1D_McL42.hdf5"); writerMcL42.set_write_norm(true); writerMcL42.set_write_energies(true);
    writerMcL42.prestructuring<MultiIndex>(packet,Dt);
	pMcL42.evolve(T,Dt,std::bind(write,writerMcL42,_1,_2));
    writerMcL42.poststructuring();
	}

	// McL84
	{
    io::hdf5writer<D> writerMcL84("morse_1D_McL84.hdf5"); writerMcL84.set_write_norm(true); writerMcL84.set_write_energies(true);
    writerMcL84.prestructuring<MultiIndex>(packet,Dt);
	pMcL84.evolve(T,Dt,std::bind(write,writerMcL84,_1,_2));
    writerMcL84.poststructuring();
	}

	// Pre764
	{
    io::hdf5writer<D> writerPre764("morse_1D_Pre764.hdf5"); writerPre764.set_write_norm(true); writerPre764.set_write_energies(true);
    writerPre764.prestructuring<MultiIndex>(packet,Dt);
	pPre764.evolve(T,Dt,std::bind(write,writerPre764,_1,_2));
    writerPre764.poststructuring();
	}

    return 0;
}