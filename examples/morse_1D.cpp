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
		const real_t Dt;
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

		io::hdf5writer<D> writer;

		Parameters_Morse_1D(std::string name)
		 : sigma_x(1.)
		 , T(50)
		 , Dt(0.005)
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
		 , writer("morse_1D_" + name + ".hdf5")
		{
			coeffs[0] = 1.0;
			packet.eps() = eps;
			packet.shape() = std::make_shared<wavepackets::shapes::ShapeEnum<D,MultiIndex>>(shape_enum);
			packet.parameters() = param_set;
			packet.coefficients() = coeffs;

			writer.set_write_norm(true);
			writer.set_write_energies(true);
			writer.prestructuring<MultiIndex>(packet,Dt);
		}

		~Parameters_Morse_1D() {
			writer.poststructuring();
		}

		void callback(unsigned /*m*/, real_t /*t*/) {
			real_t ekin = observables::kinetic_energy<D,MultiIndex>(packet);
			real_t epot = observables::potential_energy<Potential,D,MultiIndex,QR>(packet,V);
			writer.store_packet(packet);
			writer.store_norm(packet);
			writer.store_energies(epot,ekin);
		}

};

int main() {

	using P = Parameters_Morse_1D; // all the parameters for the simulation are contained in this class

	// { // Hagedorn
	// 	Parameters_Morse_1D param_Hagedorn("Hagedorn");
	// 	propagators::HagedornPropagator<P::N,P::D,P::MultiIndex,P::QR,P::Potential,P::Packet_t> propagator_Hagedorn(param_Hagedorn.packet,param_Hagedorn.V);
	// 	propagator_Hagedorn.evolve(param_Hagedorn.T,param_Hagedorn.Dt,std::bind(&P::callback,std::ref(param_Hagedorn),_1,_2));
	// }

	// { // Semiclassical
	// 	Parameters_Morse_1D param_Semiclassical("Semiclassical");
	// 	propagators::SemiclassicalPropagator<P::N,P::D,P::MultiIndex,P::QR,P::Potential,P::Packet_t,P::SplitCoefs_t> propagator_Semiclassical(param_Semiclassical.packet,param_Semiclassical.V,param_Semiclassical.splitCoefs);
	// 	propagator_Semiclassical.evolve(param_Semiclassical.T,param_Semiclassical.Dt,std::bind(&P::callback,std::ref(param_Semiclassical),_1,_2));
	// }
	
	// { // MG4
	// 	Parameters_Morse_1D param_MG4("MG4");
	// 	propagators::MG4Propagator<P::N,P::D,P::MultiIndex,P::QR,P::Potential,P::Packet_t,P::SplitCoefs_t> propagator_MG4(param_MG4.packet,param_MG4.V,param_MG4.splitCoefs);
	// 	propagator_MG4.evolve(param_MG4.T,param_MG4.Dt,std::bind(&P::callback,std::ref(param_MG4),_1,_2));
	// }

	// { // McL42
	// 	Parameters_Morse_1D param_McL42("McL42");
	// 	propagators::McL42Propagator<P::N,P::D,P::MultiIndex,P::QR,P::Potential,P::Packet_t,P::SplitCoefs_t> propagator_McL42(param_McL42.packet,param_McL42.V,param_McL42.splitCoefs);
	// 	propagator_McL42.evolve(param_McL42.T,param_McL42.Dt,std::bind(&P::callback,std::ref(param_McL42),_1,_2));
	// }

	// { // McL84
	// 	Parameters_Morse_1D param_McL84("McL84");
	// 	propagators::McL84Propagator<P::N,P::D,P::MultiIndex,P::QR,P::Potential,P::Packet_t,P::SplitCoefs_t> propagator_McL84(param_McL84.packet,param_McL84.V,param_McL84.splitCoefs);
	// 	propagator_McL84.evolve(param_McL84.T,param_McL84.Dt,std::bind(&P::callback,std::ref(param_McL84),_1,_2));
	// }

	{ // Pre764
		Parameters_Morse_1D param_Pre764("Pre764");
		propagators::Pre764Propagator<P::N,P::D,P::MultiIndex,P::QR,P::Potential,P::Packet_t,P::SplitCoefs_t> propagator_Pre764(param_Pre764.packet,param_Pre764.V,param_Pre764.splitCoefs);
		propagator_Pre764.evolve(param_Pre764.T,param_Pre764.Dt,std::bind(&P::callback,std::ref(param_Pre764),_1,_2));
	}


    return 0;
}
