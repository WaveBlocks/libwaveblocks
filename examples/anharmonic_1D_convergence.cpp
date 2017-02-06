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

	////////////////////////////////////////////////////
    // General parameters
	////////////////////////////////////////////////////

    const int N = 1;
    const int D = 1;
    const int K = 128;

    const real_t T = 1;
    real_t Dt;

    const real_t eps = 0.1;
	

	////////////////////////////////////////////////////
    // Composed Types
	////////////////////////////////////////////////////

    using MultiIndex = wavepackets::shapes::TinyMultiIndex<unsigned short, D>;
    using QR = innerproducts::GaussHermiteQR<K+4>;
	using Packet_t = wavepackets::ScalarHaWp<D,MultiIndex>;
	using Potential_t = Remain;


	/////////////////////////////////////////////////////
	// Building the Wave Packet
	/////////////////////////////////////////////////////
	
    // basis shapes
    wavepackets::shapes::ShapeEnumerator<D, MultiIndex> enumerator;
    wavepackets::shapes::ShapeEnum<D, MultiIndex> shape_enum = enumerator.generate(wavepackets::shapes::HyperCubicShape<D>(K));

    // initial parameters
	wavepackets::HaWpParamSet<D> param_set;
	param_set.p(RVector<D>::Ones()); // set p - give momentum

	// initial coefficients
    // Gaussian wave packet phi_0 with c_0 = 1
    Coefficients coefs = Coefficients::Zero(std::pow(K,D),1);
    coefs[0] = 1.0;

	
	/////////////////////////////////////////////////////
    // Defining the potential
	/////////////////////////////////////////////////////

    Potential_t V;


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
	gold_packet.coefficients() = coefs;

	// define the high-precision propagator and evolve
	propagators::SemiclassicalPropagator<N,D,MultiIndex,QR,Potential_t,Packet_t,propagators::SplitCoefs<34,34>> pGold(gold_packet,V,split::coefKL10);
	pGold.evolve(T,Dt);

	// evaluate gold solution on grid
	const CArray<1,G> eval_gold = gold_packet.evaluate(grid);


	// Convergence study
	// Start with large stepsize, then gradually refine
    
	Dt = 1;
	for(unsigned i=0; i<4; ++i) {

		Dt /= 2;

		// assemble new packet
		Packet_t packet;
		packet.eps() = eps;
		packet.shape() = std::make_shared<wavepackets::shapes::ShapeEnum<D,MultiIndex>>(shape_enum);
		packet.parameters() = param_set;
		packet.coefficients() = coefs;

	
		// Defining the Propagator
		propagators::HagedornPropagator<N,D,MultiIndex,QR,Potential_t,Packet_t> pHagedorn(packet,V);
		propagators::SemiclassicalPropagator<N,D,MultiIndex,QR,Potential_t,Packet_t,propagators::SplitCoefs<1,1>> pSemiclassical(packet,V,split::coefLT);
		propagators::MG4Propagator<N,D,MultiIndex,QR,Potential_t,Packet_t,propagators::SplitCoefs<1,1>> pMG4(packet,V,split::coefLT);
		propagators::Pre764Propagator<N,D,MultiIndex,QR,Potential_t,Packet_t,propagators::SplitCoefs<1,1>> pPre764(packet,V,split::coefLT);
		propagators::McL42Propagator<N,D,MultiIndex,QR,Potential_t,Packet_t,propagators::SplitCoefs<1,1>> pMcL42(packet,V,split::coefLT);
		propagators::McL84Propagator<N,D,MultiIndex,QR,Potential_t,Packet_t,propagators::SplitCoefs<1,1>> pMcL84(packet,V,split::coefLT);

		// evolve
		pHagedorn.evolve(T,Dt);
		// pSemiclassical.evolve(T,Dt);
		// pMG4.evolve(T,Dt);
		// pPre764.evolve(T,Dt);
		// pMcL42.evolve(T,Dt);
		// pMcL84.evolve(T,Dt);

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

}
