#include <iostream>
#include <fstream>

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

int main() {

	////////////////////////////////////////////////////
    // General parameters
	////////////////////////////////////////////////////

    const int N = 1;
    const int D = 2;
    const int K = 4;

    const real_t T = 12;
    const real_t Dt = 0.01;

    const real_t eps = 0.1;


	////////////////////////////////////////////////////
    // Composed Types
	////////////////////////////////////////////////////

    using MultiIndex = wavepackets::shapes::TinyMultiIndex<unsigned long, D>;
    using TQR = innerproducts::TensorProductQR<innerproducts::GaussHermiteQR<K+4>,
                                               innerproducts::GaussHermiteQR<K+4>>;
	using Packet_t = wavepackets::ScalarHaWp<D,MultiIndex>;
	using Potential_t = ScalarMatrixPotential<D>;


	/////////////////////////////////////////////////////
	// Building the WavePacket
	/////////////////////////////////////////////////////

    // basis shapes
    wavepackets::shapes::ShapeEnumerator<D, MultiIndex> enumerator;
    wavepackets::shapes::ShapeEnum<D, MultiIndex> shape_enum = enumerator.generate(wavepackets::shapes::HyperCubicShape<D>(K));

    // initial parameters
    CMatrix<D,D> Q = CMatrix<D,D>::Identity();
    CMatrix<D,D> P = complex_t(0,1) * CMatrix<D,D>::Identity();
    // TODO: give right values
    const float_t norm_q = .5;
    const float_t norm_p = .8;
    RVector<D> q = RVector<D>::Zero();
    RVector<D> p = RVector<D>::Zero();
	for(unsigned i=0; i<D; ++i) {
		q[i] = 1. * std::sqrt(norm_q/D); // scale the 2-norm of q to norm_q
		p[i] = (i<D/2 ? +1. : -1.) * std::sqrt(norm_p/D); // scale the 2-norm of p to norm_p
	}
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

    typename CanonicalBasis<N,D>::potential_type potential =
        [](CVector<D> x) {
			complex_t sum = 0;
			for(unsigned i=0; i<D; ++i) {
				sum += x[i]*x[i];
			}
        return sum;
    };
    typename ScalarLeadingLevel<D>::potential_type leading_level = potential;
    typename ScalarLeadingLevel<D>::jacobian_type leading_jac =
        [](CVector<D> x) {
        return 2.*x;
    };
    typename ScalarLeadingLevel<D>::hessian_type leading_hess =
        [](CVector<D> /*x*/) {
        return 2.*CMatrix<D,D>::Identity();
    };

    Potential_t V(potential,leading_level,leading_jac,leading_hess);


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
    io::hdf5writer<D> mywriter("harmonic_XD_newpropagator.hdf5");
    mywriter.set_write_norm(true);
    mywriter.set_write_energies(true);
	
	// Write Data Callback
	std::function<void(unsigned,real_t)> writeenergies = [&](unsigned,real_t) {
		real_t ekin = observables::kinetic_energy<D,MultiIndex>(packet);
		real_t epot = observables::potential_energy<Potential_t,D,MultiIndex,TQR>(packet,V);
		mywriter.store_packet(packet);
		mywriter.store_norm(packet);
		mywriter.store_energies(epot,ekin);
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
