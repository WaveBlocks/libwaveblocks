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

    const real_t sigma_x = 0.5;
    const real_t sigma_y = 0.5;

    const real_t T = 100;
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
    RVector<D> q = {-3.0, 0.0};
    RVector<D> p = { 0.0, 0.5};
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
        [sigma_x,sigma_y](CVector<D> x) {
        return 0.5*(sigma_x*x[0]*x[0] + sigma_y*x[1]*x[1]).real();
    };
    typename ScalarLeadingLevel<D>::potential_type leading_level = potential;
    typename ScalarLeadingLevel<D>::jacobian_type leading_jac =
        [sigma_x,sigma_y](CVector<D> x) {
        return CVector<D>{sigma_x*x[0], sigma_y*x[1]};
    };
    typename ScalarLeadingLevel<D>::hessian_type leading_hess =
        [sigma_x,sigma_y](CVector<D> /*x*/) {
        CMatrix<D,D> res;
        res(0,0) = sigma_x;
        res(1,1) = sigma_y;
        return res;
    };

    Potential_t V(potential,leading_level,leading_jac,leading_hess);


	////////////////////////////////////////////////////
	// Defining the Propagator
	////////////////////////////////////////////////////

	// propagators::HagedornPropagator<N,D,MultiIndex,TQR,Potential_t,Packet_t> pHagedorn(packet,V);
	propagators::SemiclassicalPropagator<N,D,MultiIndex,TQR,Potential_t,Packet_t,propagators::SplitCoefs<1,1>> pSemiclassical(packet,V,split::coefLT);
	// propagators::MG4Propagator<N,D,MultiIndex,TQR,Potential_t,Packet_t,propagators::SplitCoefs<1,1>> pMG4(packet,V,split::coefLT);
	// propagators::Pre764Propagator<N,D,MultiIndex,TQR,Potential_t,Packet_t,propagators::SplitCoefs<1,1>> pPre764(packet,V,split::coefLT);
	// propagators::McL42Propagator<N,D,MultiIndex,TQR,Potential_t,Packet_t,propagators::SplitCoefs<1,1>> pMcL42(packet,V,split::coefLT);
	// propagators::McL84Propagator<N,D,MultiIndex,TQR,Potential_t,Packet_t,propagators::SplitCoefs<1,1>> pMcL84(packet,V,split::coefLT);


	////////////////////////////////////////////////////
	// Defining Callback Function
	////////////////////////////////////////////////////

	// set up writer: preparing the file and I/O writer
    io::hdf5writer<D> mywriter("harmonic_2D_newpropagator.hdf5");
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

	// pHagedorn.evolve(T,Dt,writeenergies);
	pSemiclassical.evolve(T,Dt);
	// pMG4.evolve(T,Dt,writeenergies);
	// pPre764.evolve(T,Dt,writeenergies);
	// pMcL42.evolve(T,Dt,writeenergies);
	// pMcL84.evolve(T,Dt,writeenergies);

    mywriter.poststructuring();

    return 0;
}
