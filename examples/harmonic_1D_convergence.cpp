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


using namespace waveblocks;
using namespace std::placeholders;


class Parameters_Harmonic_1D {

	public:

		static const int N = 1;
		static const int D = 1;
		static const int K = 16;
		
		class Potential : public potentials::modules::evaluation::Abstract<Potential,CanonicalBasis<N,D>>,		                   
		                   public potentials::modules::taylor::Abstract<Potential,CanonicalBasis<N,D>>,
		                   public potentials::modules::localRemainder::Abstract<Potential,N,D>,
		                   public LeadingLevelOwner<potentials::modules::taylor::Abstract<Potential,CanonicalBasis<N,D>>> {

			public:

				using Taylor = potentials::modules::taylor::Abstract<Potential,CanonicalBasis<N,D>>;

				inline Taylor::potential_evaluation_type evalV(const Taylor::argument_type& x) const {
					return 0.5*(x*x).real();
				}
				
				inline Taylor::jacobian_evaluation_type evalJ(const Taylor::argument_type& x) const {
					return x;
				}
				
				inline Taylor::hessian_evaluation_type evalH(const Taylor::argument_type& /*x*/) const {
					return 1.;
				}

				complex_t evaluate_at_implementation(const Taylor::argument_type& x) const {
					return evalV(x);
				}

				template <template <typename...> class Tuple = std::tuple>
					Tuple<Taylor::potential_evaluation_type, Taylor::jacobian_evaluation_type, Taylor::hessian_evaluation_type>
					taylor_at_implementation( const Taylor::argument_type &x ) const {
						return Tuple<Taylor::potential_evaluation_type,Taylor::jacobian_evaluation_type,Taylor::hessian_evaluation_type>
							(evalV(x), evalJ(x), evalH(x));
				}

				complex_t evaluate_local_remainder_at( const complex_t &x, const complex_t &q ) const {
					const auto xmq = x - q;
					return evalV(x) - evalV(q) - evalJ(q)*xmq - 0.5*xmq*evalH(q)*xmq;
				}

		};


		using MultiIndex = wavepackets::shapes::TinyMultiIndex<unsigned short, D>;
		using QR = innerproducts::GaussHermiteQR<K+4>;
		using Packet_t = wavepackets::ScalarHaWp<D,MultiIndex>;
		using SplitCoefs_t = propagators::SplitCoefs<1,1>;

		// general parameters
		const real_t sigma_x;
		const real_t T;
		const real_t eps;

		// wave packet
		Potential V;
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
		 , eps(0.01)
		 , Q(CMatrix<D,D>::Identity())
		 , P(complex_t(0,1) * CMatrix<D,D>::Identity())
		 , q(RVector<D>::Ones())
		 , p(RVector<D>::Zero())
		 , S(0.)
		 , param_set(q,p,Q,P,S)
		 , coeffs(Coefficients::Zero(std::pow(K,D),1))
		 , splitCoefs(propagators::splitting_parameters::coefLT)
		 , shape_enum(enumerator.generate(wavepackets::shapes::HyperCubicShape<D>(K)))
		{
			coeffs[0] = 1.0;
			packet.eps() = eps;
			packet.shape() = std::make_shared<wavepackets::shapes::ShapeEnum<D,MultiIndex>>(shape_enum);
			packet.parameters() = param_set;
			packet.coefficients() = coeffs;
		}

};

int main() {

	using P = Parameters_Harmonic_1D;

	const real_t Dt_gold = 1e-4;

	// define a grid for evaluation
    const int G = 1000;
    RMatrix<P::D,G> grid(P::D,G);
    for (int d = 0; d < P::D; d++) {
        for (int i = 0; i < G; i++) {
            grid(d,i) = 5*(-1.0 + 2.0*(i+1)/real_t(G));
        }
    }

	// define the high-precision propagator and evolve
	Parameters_Harmonic_1D param_gold;
	propagators::McL84Propagator<P::N,P::D,P::MultiIndex,P::QR,P::Potential,P::Packet_t,propagators::SplitCoefs<34,34>> pGold(param_gold.packet,param_gold.V,propagators::splitting_parameters::coefKL10);
	pGold.evolve(param_gold.T,Dt_gold);

	// evaluate gold solution on grid
	const CArray<1,G> eval_gold = param_gold.packet.evaluate(grid);
	real_t norm_gold = 0;
	for(unsigned i=0; i<G; i++) {
		real_t x = std::abs(eval_gold[i]);
		norm_gold += x*x;
	}
	norm_gold = std::sqrt(norm_gold);

	auto errorL2 = [&] (P::Packet_t& pack) {
		const CArray<1,G> eval_approx = pack.evaluate(grid);
		real_t error = 0;
		for(unsigned i=0; i<G; i++) {
			real_t diff = std::abs(eval_approx[i]) - std::abs(eval_gold[i]);
			error += diff*diff;
		}
		return std::sqrt(error)/norm_gold;
	};


	real_t Dt = .25;

	while(Dt>=2*Dt_gold) {

		{ // Hagedorn
			Parameters_Harmonic_1D param_Hagedorn;
			propagators::HagedornPropagator<P::N,P::D,P::MultiIndex,P::QR,P::Potential,P::Packet_t> propagator_Hagedorn(param_Hagedorn.packet,param_Hagedorn.V);
			propagator_Hagedorn.evolve(param_Hagedorn.T,Dt);
			std::cout << "\nHagedorn - Dt: " << std::fixed << Dt << "\tError: " << std::scientific << errorL2(param_Hagedorn.packet) << "\n";
		}
		
		{ // Semiclassical
			Parameters_Harmonic_1D param_Semiclassical;
			propagators::SemiclassicalPropagator<P::N,P::D,P::MultiIndex,P::QR,P::Potential,P::Packet_t,P::SplitCoefs_t> propagator_Semiclassical(param_Semiclassical.packet,param_Semiclassical.V,param_Semiclassical.splitCoefs);
			propagator_Semiclassical.evolve(param_Semiclassical.T,Dt);
			std::cout << "\nSemiclassical - Dt: " << std::fixed << Dt << "\tError: " << std::scientific << errorL2(param_Semiclassical.packet) << "\n";
		}

		{ // MG4
			Parameters_Harmonic_1D param_MG4;
			propagators::MG4Propagator<P::N,P::D,P::MultiIndex,P::QR,P::Potential,P::Packet_t,P::SplitCoefs_t> propagator_MG4(param_MG4.packet,param_MG4.V,param_MG4.splitCoefs);
			propagator_MG4.evolve(param_MG4.T,Dt);
			std::cout << "\nMG4 - Dt: " << std::fixed << Dt << "\tError: " << std::scientific << errorL2(param_MG4.packet) << "\n";
		}

		{ // McL42
			Parameters_Harmonic_1D param_McL42;
			propagators::McL42Propagator<P::N,P::D,P::MultiIndex,P::QR,P::Potential,P::Packet_t,P::SplitCoefs_t> propagator_McL42(param_McL42.packet,param_McL42.V,param_McL42.splitCoefs);
			propagator_McL42.evolve(param_McL42.T,Dt);
			std::cout << "\nMcL42 - Dt: " << std::fixed << Dt << "\tError: " << std::scientific << errorL2(param_McL42.packet) << "\n";
		}

		{ // McL84
			Parameters_Harmonic_1D param_McL84;
			propagators::McL84Propagator<P::N,P::D,P::MultiIndex,P::QR,P::Potential,P::Packet_t,P::SplitCoefs_t> propagator_McL84(param_McL84.packet,param_McL84.V,param_McL84.splitCoefs);
			propagator_McL84.evolve(param_McL84.T,Dt);
			std::cout << "\nMcL84 - Dt: " << std::fixed << Dt << "\tError: " << std::scientific << errorL2(param_McL84.packet) << "\n";
		}

		{ // Pre764
			Parameters_Harmonic_1D param_Pre764;
			propagators::Pre764Propagator<P::N,P::D,P::MultiIndex,P::QR,P::Potential,P::Packet_t,P::SplitCoefs_t> propagator_Pre764(param_Pre764.packet,param_Pre764.V,param_Pre764.splitCoefs);
			propagator_Pre764.evolve(param_Pre764.T,Dt);
			std::cout << "\nPre764 - Dt: " << std::fixed << Dt << "\tError: " << std::scientific << errorL2(param_Pre764.packet) << "\n";
		}

		Dt /= 2;

	}

    return 0;
}
