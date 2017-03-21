#include <fstream>
#include <functional>
#include <iostream>

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
#include "waveblocks/propagators/McL84Propagator.hpp"
#include "waveblocks/utilities/prettyprint.hpp"


using namespace waveblocks;
using namespace std::placeholders;
namespace print = utilities::prettyprint;


class Parameters_Morse_1D {

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
					return 0.004164 * (std::exp(-2.*0.896696*(x-5.542567)) - 2.*std::exp(-0.896696*(x-5.542567)));
				}
				
				inline Taylor::jacobian_evaluation_type evalJ(const Taylor::argument_type& x) const {
					return 0.004164 * (-2.*0.896696*std::exp(-2.*0.896696*(x-5.542567)) + 2.*0.896696*std::exp(-0.896696*(x-5.542567)));
				}
				
				inline Taylor::hessian_evaluation_type evalH(const Taylor::argument_type& x) const {
					return 0.004164 * (4.*0.896696*0.896696*std::exp(-2.*0.896696*(x-5.542567)) - 2.*0.896696*0.896696*std::exp(-0.896696*(x-5.542567)));
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
		using SplitCoefs_t = propagators::SplitCoefs<4,4>;

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
	
		Parameters_Morse_1D()
		 : sigma_x(1.)
		 , T(50)
		 , eps(0.0484)
		 , Q(complex_t(3.4957,0.) * CMatrix<D,D>::Identity())
		 , P(complex_t(0.,0.2861) * CMatrix<D,D>::Identity())
		 , q(6.0426 * RVector<D>::Ones())
		 , p(-0.1100 * RVector<D>::Ones())
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

int main() {

	using P = Parameters_Morse_1D;

	const real_t Dt_gold = 5e-4;

	// define a grid for evaluation
    const int G = 1000;
    RMatrix<P::D,G> grid(P::D,G);
    for (int d = 0; d < P::D; d++) {
        for (int i = 0; i < G; i++) {
            grid(d,i) = 5+5*(-1.0 + 2.0*(i+1)/real_t(G));
        }
    }

	// define the high-precision propagator and evolve
	P param_gold;
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

	std::ofstream file_Hagedorn, file_Semiclassical, file_MG4, file_McL42, file_McL84, file_Pre764;
	file_Hagedorn.open("morse_1D_error_Hagedorn.csv");
	file_Semiclassical.open("morse_1D_error_Semiclassical.csv");
	file_MG4.open("morse_1D_error_MG4.csv");
	file_McL42.open("morse_1D_error_McL42.csv");
	file_McL84.open("morse_1D_error_McL84.csv");
	file_Pre764.open("morse_1D_error_Pre764.csv");
	file_Hagedorn << "# Propagator, split coefs, T\nHagedorn" << ", " << "Y4" << ", " << param_gold.T << "\n\n# Step size, error\n";
	file_Semiclassical << "# Propagator, split coefs, T\nSemiclassical" << ", " << "Y4" << ", " << param_gold.T << "\n\n# Step size, error\n";
	file_MG4 << "# Propagator, split coefs, T\nMG4" << ", " << "Y4" << ", " << param_gold.T << "\n\n# Step size, error\n";
	file_McL42 << "# Propagator, split coefs, T\nMcL42" << ", " << "Y4" << ", " << param_gold.T << "\n\n# Step size, error\n";
	file_McL84 << "# Propagator, split coefs, T\nMcL84" << ", " << "Y4" << ", " << param_gold.T << "\n\n# Step size, error\n";
	file_Pre764 << "# Propagator, split coefs, T\nPre764" << ", " << "Y4" << ", " << param_gold.T << "\n\n# Step size, error\n";

	while(Dt>=2*Dt_gold) {

		// file for error measurements

		{ // Hagedorn
			P param_Hagedorn;
			propagators::HagedornPropagator<P::N,P::D,P::MultiIndex,P::QR,P::Potential,P::Packet_t> propagator_Hagedorn(param_Hagedorn.packet,param_Hagedorn.V);
			propagator_Hagedorn.evolve(param_Hagedorn.T,Dt);
			real_t err = errorL2(param_Hagedorn.packet);
			file_Hagedorn << Dt << ", " << err << "\n";
			std::cout << "\n\tError: " << std::scientific << err << "\n";
		}
		
		{ // Semiclassical
			P param_Semiclassical;
			propagators::SemiclassicalPropagator<P::N,P::D,P::MultiIndex,P::QR,P::Potential,P::Packet_t,P::SplitCoefs_t> propagator_Semiclassical(param_Semiclassical.packet,param_Semiclassical.V,param_Semiclassical.splitCoefs);
			propagator_Semiclassical.evolve(param_Semiclassical.T,Dt);
			real_t err = errorL2(param_Semiclassical.packet);
			file_Semiclassical << Dt << ", " << err << "\n";
			std::cout << "\n\tError: " << std::scientific << err << "\n";
		}

		{ // MG4
			P param_MG4;
			propagators::MG4Propagator<P::N,P::D,P::MultiIndex,P::QR,P::Potential,P::Packet_t,P::SplitCoefs_t> propagator_MG4(param_MG4.packet,param_MG4.V,param_MG4.splitCoefs);
			propagator_MG4.evolve(param_MG4.T,Dt);
			real_t err = errorL2(param_MG4.packet);
			file_MG4 << Dt << ", " << err << "\n";
			std::cout << "\n\tError: " << std::scientific << err << "\n";
		}

		{ // McL42
			P param_McL42;
			propagators::McL42Propagator<P::N,P::D,P::MultiIndex,P::QR,P::Potential,P::Packet_t,P::SplitCoefs_t> propagator_McL42(param_McL42.packet,param_McL42.V,param_McL42.splitCoefs);
			propagator_McL42.evolve(param_McL42.T,Dt);
			real_t err = errorL2(param_McL42.packet);
			file_McL42 << Dt << ", " << err << "\n";
			std::cout << "\n\tError: " << std::scientific << err << "\n";
		}

		{ // McL84
			P param_McL84;
			propagators::McL84Propagator<P::N,P::D,P::MultiIndex,P::QR,P::Potential,P::Packet_t,P::SplitCoefs_t> propagator_McL84(param_McL84.packet,param_McL84.V,param_McL84.splitCoefs);
			propagator_McL84.evolve(param_McL84.T,Dt);
			real_t err = errorL2(param_McL84.packet);
			file_McL84 << Dt << ", " << err << "\n";
			std::cout << "\n\tError: " << std::scientific << err << "\n";
		}


		{ // Pre764
			P param_Pre764;
			propagators::Pre764Propagator<P::N,P::D,P::MultiIndex,P::QR,P::Potential,P::Packet_t,P::SplitCoefs_t> propagator_Pre764(param_Pre764.packet,param_Pre764.V,param_Pre764.splitCoefs);
			propagator_Pre764.evolve(param_Pre764.T,Dt);
			real_t err = errorL2(param_Pre764.packet);
			file_Pre764 << Dt << ", " << err << "\n";
			std::cout << "\n\tError: " << std::scientific << err << "\n";
		}
	
		Dt /= 2;

	}
		
	file_Hagedorn.close();
	file_Semiclassical.close();
	file_MG4.close();
	file_McL42.close();
	file_McL84.close();
	file_Pre764.close();

    return 0;
}
