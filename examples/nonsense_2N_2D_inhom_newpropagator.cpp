#include <iostream>
#include <fstream>

#include "waveblocks/types.hpp"
#include "waveblocks/propagators/HagedornPropagator.hpp"
#include "waveblocks/propagators/SemiclassicalPropagator.hpp"
#include "waveblocks/propagators/MG4Propagator.hpp"
#include "waveblocks/propagators/Pre764Propagator.hpp"
#include "waveblocks/propagators/McL42Propagator.hpp"
#include "waveblocks/propagators/McL84Propagator.hpp"
#include "waveblocks/potentials/potentials.hpp"
#include "waveblocks/potentials/bases.hpp"
#include "waveblocks/wavepackets/shapes/tiny_multi_index.hpp"
#include "waveblocks/wavepackets/hawp_commons.hpp"
#include "waveblocks/wavepackets/hawp_paramset.hpp"
#include "waveblocks/wavepackets/shapes/shape_enumerator.hpp"
#include "waveblocks/wavepackets/shapes/shape_hypercubic.hpp"
#include "waveblocks/innerproducts/gauss_hermite_qr.hpp"
#include "waveblocks/innerproducts/tensor_product_qr.hpp"


using namespace waveblocks;

namespace split = propagators::splitting_parameters;

int main() {

	////////////////////////////////////////////////////
    // General parameters
	////////////////////////////////////////////////////

    const int N = 2;
    const int D = 2;
    const int K = 3;
    const real_t sigma_x = 0.5;
    const real_t sigma_y = 0.5;

    const real_t T = 12;
    const real_t Dt = 0.01;

    const real_t eps = 0.1;


	////////////////////////////////////////////////////
    // Composed Types
	////////////////////////////////////////////////////

    using MultiIndex = wavepackets::shapes::TinyMultiIndex<unsigned short, D>;
    using TQR = innerproducts::TensorProductQR<innerproducts::GaussHermiteQR<3>,
                                               innerproducts::GaussHermiteQR<4>>;
    using Packet_t = wavepackets::InhomogeneousHaWp<D,MultiIndex>;
    using Potential_t = InhomogenousMatrixPotential<N,D>;


	/////////////////////////////////////////////////////
	// Building the WavePacket
	/////////////////////////////////////////////////////

    // basis shapes
    wavepackets::shapes::ShapeEnumerator<D, MultiIndex> enumerator;
    wavepackets::shapes::ShapeEnum<D, MultiIndex> shape_enum = enumerator.generate(wavepackets::shapes::HyperCubicShape<D>(K));

    // initial parameters
    CMatrix<D,D> Q = CMatrix<D,D>::Identity();
    CMatrix<D,D> P = complex_t(0,1)*CMatrix<D,D>::Identity();
    RVector<D> q = {-3.0, 0.0};
    RVector<D> p = {0.0, 0.5};
    wavepackets::HaWpParamSet<D> param_set(q, p, Q, P, 0.0);
    wavepackets::HaWpParamSet<D> param_set2(2*q, 0.5*p, Q, P, 0.0);

    // initial coefficients
    Coefficients coeffs = Coefficients::Ones(std::pow(K, D), 1);
    Packet_t packet(N);

    // assemble packet
    packet.eps() = eps;
    packet.component(0).shape() = std::make_shared<wavepackets::shapes::ShapeEnum<D,MultiIndex>>(shape_enum);
    packet.component(1).shape() = std::make_shared<wavepackets::shapes::ShapeEnum<D,MultiIndex>>(shape_enum);
    packet.component(0).parameters() = param_set;
    packet.component(1).parameters() = param_set2;
    packet.component(0).coefficients() = coeffs;
    packet.component(1).coefficients() = coeffs;


	/////////////////////////////////////////////////////
    // Defining the potential
	/////////////////////////////////////////////////////

    typename CanonicalBasis<N,D>::potential_type potential;

    potential(0,0) = [sigma_x,sigma_y,N](CVector<D> x) {
        return 0.5*(sigma_x*x[0]*x[0] + sigma_y*x[1]*x[1]).real();
    };

    potential(0,1) = [sigma_x,sigma_y,N](CVector<D> x) {
        return 0.5*(sigma_x*x[0]*x[0] + sigma_y*x[1]*x[1]).real();
    };
    potential(1,0) = [sigma_x,sigma_y,N](CVector<D> x) {
        return 0.5*(sigma_x*x[0]*x[0] + sigma_y*x[1]*x[1]).real();
    };
    potential(1,1) = [sigma_x,sigma_y,N](CVector<D> x) {
        return 0.5*(sigma_x*x[0]*x[0] + sigma_y*x[1]*x[1]).real();
    };

    typename InhomogenousLeadingLevel<N,D>::potential_type leading_level;
    leading_level(0) = [sigma_x,sigma_y,N](CVector<D> x) {
        return 0.5*(sigma_x*x[0]*x[0] + sigma_y*x[1]*x[1]).real();
    };

    typename InhomogenousLeadingLevel<N,D>::jacobian_type leading_jac;
    leading_jac(0) = [sigma_x,sigma_y,N](CVector<D> x) {
        return CVector<D>{sigma_x*x[0], sigma_y*x[1]};
    };

    typename InhomogenousLeadingLevel<N,D>::hessian_type leading_hess;
    leading_hess(0) =
        [sigma_x,sigma_y,N](CVector<D> x) {
        (void) x; // Unused
        CMatrix<D,D> res;
        res(0,0) = sigma_x;
        res(1,1) = sigma_y;
        return res;
    };

    leading_level(1) = leading_level(0);
    leading_jac(1) = leading_jac(0);
    leading_hess(1) = leading_hess(0);

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
	// Propagate
	////////////////////////////////////////////////////

	pHagedorn.evolve(T,Dt);
	pSemiclassical.evolve(T,Dt);
	pMG4.evolve(T,Dt);
	pPre764.evolve(T,Dt);
	pMcL42.evolve(T,Dt);
	pMcL84.evolve(T,Dt);

    return 0;
}
