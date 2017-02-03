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

    const real_t sigma_harm = 1.;
    const real_t sigma_mix = 0.2;

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
    CMatrix<D,D> Q = CMatrix<D,D>::Zero();
    float_t a = std::sqrt(2.*0.56);
    float_t b = std::sqrt(2.*0.24); 
    Q(0,0) = complex_t(a);
    Q(1,1) = complex_t(b);
    CMatrix<D,D> P = CMatrix<D,D>::Zero();
    P(0,0) = complex_t(0,1./a);
    P(1,1) = complex_t(0,1./b);
    RVector<D> q = {1.8, 0.0};
    RVector<D> p = { 0.0, 1.2};
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
        [sigma_harm,sigma_mix](CVector<D> x) {
			complex_t sum_harm = 0;
        	complex_t sum_mix = 0;
			for(unsigned j=0; j<D; ++j) {
				sum_harm += sigma_harm * x[j] * x[j];
			}
			for(unsigned j=0; j<D-1; ++j) {
				sum_mix += (
					x[j] * (x[j+1]*x[j+1] - 1./3.*x[j]*x[j]) +
					0.0625 * sigma_mix * (x[j]*x[j] + x[j+1]*x[j+1]) * (x[j]*x[j] + x[j+1]*x[j+1])
				);
			}
        return .5*sum_harm + sigma_mix*sum_mix;
    };
    typename ScalarLeadingLevel<D>::potential_type leading_level = potential;
    typename ScalarLeadingLevel<D>::jacobian_type leading_jac =
        [sigma_harm,sigma_mix](CVector<D> x) {
			auto jac = [&](unsigned i) {
				return sigma_harm * x[i] +
					sigma_mix * (
						( i==0 ? 0. : 2.*x[i-1]*x[i] ) -
						x[i]*x[i] +
						( i==D-1 ? 0. : x[i+1]*x[i+1]) +
						.25*sigma_mix*x[i] * (
							( i==0 ? 0. : x[i-1]*x[i-1] ) +
							2.*x[i]*x[i] +
							( i==D-1 ? 0. : x[i+1]*x[i+1] )
						)
					);
			};
			CVector<D> res;
			for(unsigned i=0; i<D; ++i) {
				res[i] = jac(i);
			}
        return res;
    };
    typename ScalarLeadingLevel<D>::hessian_type leading_hess =
        [sigma_harm,sigma_mix](CVector<D> x) {

		auto hess = [&](unsigned k, unsigned l) {
			if(k==l-1) {
				return sigma_mix * (
					k==D-1 ? 0. :
					2.*x[k+1] + .5*sigma_mix*x[k+1]*x[k]
				);
			}
			if(k==l) {
				return sigma_harm + sigma_mix * (
					( k==0 ? 0. : 2.*x[k-1] ) -
					2.*x[k] + 
					.25*sigma_mix * (
						( k==0 ? 0. : x[k-1]*x[k-1] ) +
						6.*x[k]*x[k] +
						( k==D-1 ? 0. : x[k+1]*x[k+1] )
					)
				);
			}
			if(k==(l+1)) {
				return sigma_mix * (
					2.*x[k] +
					(k==0 ? 0. : .5*sigma_mix*x[k-1]*x[k])
				);
			}
			return complex_t(0.,0.);
		};

        CMatrix<D,D> res;
		for(unsigned k=0; k<D; ++k) {
			for(unsigned l=0; l<D; ++l) {
				res(k,l) = hess(k,l);
			}
		}
        return res;
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

	pHagedorn.evolve(T,Dt,writeenergies);
	// pSemiclassical.evolve(T,Dt);
	// pMG4.evolve(T,Dt);
	// pPre764.evolve(T,Dt);
	// pMcL42.evolve(T,Dt);
	// pMcL84.evolve(T,Dt);

    mywriter.poststructuring();

    return 0;
}
