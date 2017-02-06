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

    const real_t T = 10;
    real_t Dt;

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
	// Convergence Study
	////////////////////////////////////////////////////
	
	// define a grid for evaluation
    const int G = 1000;
    RMatrix<D, G> grid(D,G);
    for (int d = 0; d < D; d++) {
        for (int i = 0; i < G; i++) {
            grid(d,i) = (-1.0 + 2.0*(d+1)/float(D)) * (i+1)/G;
        }
    }

	
	// compute reference solution

	Dt = 0.001;

	// build gold packet (reference solution)
	Packet_t gold_packet;
	gold_packet.eps() = eps;
	gold_packet.shape() = std::make_shared<wavepackets::shapes::ShapeEnum<D,MultiIndex>>(shape_enum);
	gold_packet.parameters() = param_set;
	gold_packet.coefficients() = coeffs;

	// define the high-precision propagator and evolve
	propagators::SemiclassicalPropagator<N,D,MultiIndex,TQR,Potential_t,Packet_t,propagators::SplitCoefs<34,34>> pGold(gold_packet,V,split::coefKL10);
	pGold.evolve(T,Dt);

	// evaluate gold solution on grid
	const CArray<1,G> eval_gold = gold_packet.evaluate(grid);


	// Convergence study
	// Start with large stepsize, then gradually refine
    
	Dt = 1;
	for(unsigned i=0; i<10; ++i) {

		Dt /= 2;

		// assemble new packet
		Packet_t packet;
		packet.eps() = eps;
		packet.shape() = std::make_shared<wavepackets::shapes::ShapeEnum<D,MultiIndex>>(shape_enum);
		packet.parameters() = param_set;
		packet.coefficients() = coeffs;
	
		// Defining the Propagator
		propagators::HagedornPropagator<N,D,MultiIndex,TQR,Potential_t,Packet_t> pHagedorn(packet,V);
		propagators::SemiclassicalPropagator<N,D,MultiIndex,TQR,Potential_t,Packet_t,propagators::SplitCoefs<1,1>> pSemiclassical(packet,V,split::coefLT);
		propagators::MG4Propagator<N,D,MultiIndex,TQR,Potential_t,Packet_t,propagators::SplitCoefs<1,1>> pMG4(packet,V,split::coefLT);
		propagators::Pre764Propagator<N,D,MultiIndex,TQR,Potential_t,Packet_t,propagators::SplitCoefs<1,1>> pPre764(packet,V,split::coefLT);
		propagators::McL42Propagator<N,D,MultiIndex,TQR,Potential_t,Packet_t,propagators::SplitCoefs<1,1>> pMcL42(packet,V,split::coefLT);
		propagators::McL84Propagator<N,D,MultiIndex,TQR,Potential_t,Packet_t,propagators::SplitCoefs<1,1>> pMcL84(packet,V,split::coefLT);

		// evolve
		// pHagedorn.evolve(T,Dt);
		// pSemiclassical.evolve(T,Dt);
		// pMG4.evolve(T,Dt);
		// pPre764.evolve(T,Dt);
		// pMcL42.evolve(T,Dt);
		pMcL84.evolve(T,Dt);

		// Compute L2 error
		const CArray<1,G> eval_this = packet.evaluate(grid);
		real_t error = 0;
		for(unsigned i=0; i<G; ++i) {
			real_t diff = std::abs(eval_this[i] - eval_gold[i]);
			error += diff*diff;
		}
		error = std::sqrt(error);
		std::cout << "\n\tDt: " << Dt;
		std::cout << "\n\tError: " << error;

	}

    return 0;
}
