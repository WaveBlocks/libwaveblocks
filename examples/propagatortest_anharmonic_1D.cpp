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
#include "waveblocks/propagators/Propagator.hpp"
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

    // General parameters
    const int N = 1;
    const int D = 1;
    const int K = 128;

    const real_t T = 6;
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


	/////////////////////////////////////////////////////
	// Propagate
	/////////////////////////////////////////////////////

	propagators::steps::HagedornPropagator<Remain> propagator(packet,V); // <typename Remain> propagator(packet,V);
	propagator.simulate(T,Dt);


	return 0;

}
