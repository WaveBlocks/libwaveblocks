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


class Parameters_Harmonic_1D {

	public:

		static const int N = 1;
		static const int D = 1;
		static const int K = 16;
		
		using MultiIndex = wavepackets::shapes::TinyMultiIndex<unsigned short, D>;
		using QR = innerproducts::GaussHermiteQR<K+4>;
		using Packet_t = wavepackets::ScalarHaWp<D,MultiIndex>;
		using Potential_t = ScalarMatrixPotential<D>;
		using SplitCoefs_t = propagators::SplitCoefs<4,4>;

		// general parameters
		const real_t sigma_x;
		const real_t T;
		const real_t Dt;
		const real_t eps;

		// wave packet
		Packet_t packet;

		// initial values
		CMatrix<D,D> Q;
		CMatrix<D,D> P;
		RVector<D> q;
		RVector<D> p;
		complex_t S;
		wavepackets::HaWpParamSet<D> param_set;

		// initial coefficients 
		Coefficients coeffs;

		const SplitCoefs_t splitCoefs;

		// basis shapes
		wavepackets::shapes::ShapeEnumerator<D, MultiIndex> enumerator;
		wavepackets::shapes::ShapeEnum<D, MultiIndex> shape_enum;


		Parameters_Harmonic_1D()
		 : sigma_x(1.)
		 , T(10)
		 , Dt(0.05)
		 , eps(0.01)
		 , Q(CMatrix<D,D>::Identity())
		 , P(complex_t(0,1) * CMatrix<D,D>::Identity())
		 , q(RVector<D>::Ones())
		 , p(RVector<D>::Zero())
		 , S(0.)
		 , param_set(q,p,Q,P,S)
		 , coeffs(Coefficients::Zero(std::pow(K,D),1))
		 , splitCoefs(propagators::splitting_parameters::coefY4)
		 , shape_enum(enumerator.generate(wavepackets::shapes::HyperCubicShape<D>(K)))
		{
			coeffs[0] = 1.0;
			packet.eps() = eps;
			packet.shape() = std::make_shared<wavepackets::shapes::ShapeEnum<D,MultiIndex>>(shape_enum);
			packet.parameters() = param_set;
			packet.coefficients() = coeffs;
		}


};

// TODO: move more stuff into parameter class for convenience!
// TODO: keep testing on the go

int main() {

	/////////////////////////////////////////////////////
    // Defining the potential
	/////////////////////////////////////////////////////

    typename CanonicalBasis<Parameters_Harmonic_1D::N,Parameters_Harmonic_1D::D>::potential_type potential = [/*sigma_x*/](complex_t x) {
		return 0.5*(/*sigma_x*/x*x).real();
    };
    typename ScalarLeadingLevel<Parameters_Harmonic_1D::D>::potential_type leading_level = potential;
    typename ScalarLeadingLevel<Parameters_Harmonic_1D::D>::jacobian_type leading_jac = [/*sigma_x*/](complex_t x) {
        return /*sigma_x*/x;
    };
    typename ScalarLeadingLevel<Parameters_Harmonic_1D::D>::hessian_type leading_hess = [/*sigma_x*/](complex_t /*x*/) {
        return 0; // sigma_x;
    };

	Parameters_Harmonic_1D::Potential_t V(potential,leading_level,leading_jac,leading_hess);


	////////////////////////////////////////////////////
	// Propagate
	////////////////////////////////////////////////////

	// Semiclassical
	{
	Parameters_Harmonic_1D param;
	auto write = [&](io::hdf5writer<Parameters_Harmonic_1D::D>& writer,unsigned,real_t) {
		real_t ekin = observables::kinetic_energy<Parameters_Harmonic_1D::D,Parameters_Harmonic_1D::MultiIndex>(param.packet);
		real_t epot = observables::potential_energy<Parameters_Harmonic_1D::Potential_t,Parameters_Harmonic_1D::D,Parameters_Harmonic_1D::MultiIndex,Parameters_Harmonic_1D::QR>(param.packet,V);
		writer.store_packet(param.packet);
		writer.store_norm(param.packet);
		writer.store_energies(epot,ekin);
	};
	propagators::SemiclassicalPropagator<Parameters_Harmonic_1D::N,Parameters_Harmonic_1D::D,Parameters_Harmonic_1D::MultiIndex,Parameters_Harmonic_1D::QR,Parameters_Harmonic_1D::Potential_t,Parameters_Harmonic_1D::Packet_t,Parameters_Harmonic_1D::SplitCoefs_t> pSemiclassical(param.packet,V,param.splitCoefs);
    io::hdf5writer<Parameters_Harmonic_1D::D> writerSemiclassical("harmonic_1D_Semiclassical.hdf5"); writerSemiclassical.set_write_norm(true); writerSemiclassical.set_write_energies(true);
    writerSemiclassical.prestructuring<Parameters_Harmonic_1D::MultiIndex>(param.packet,param.Dt);
	pSemiclassical.evolve(param.T,param.Dt,std::bind(write,std::ref(writerSemiclassical),_1,_2));
    writerSemiclassical.poststructuring();
	}

    return 0;
}
