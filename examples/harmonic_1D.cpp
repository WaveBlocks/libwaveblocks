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

		Parameters_Harmonic_1D(std::string name)
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
		 , writer("harmonic_1D_" + name + ".hdf5")
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

		~Parameters_Harmonic_1D() {
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

// TODO: move more stuff into parameter class for convenience!
// TODO: keep testing on the go

int main() {

	using P = Parameters_Harmonic_1D; // all the parameters for the simulation are contained in this class

	{ // Semiclassical
		Parameters_Harmonic_1D param("Semiclassical");
		propagators::SemiclassicalPropagator<P::N,P::D,P::MultiIndex,P::QR,P::Potential,P::Packet_t,P::SplitCoefs_t> pSemiclassical(param.packet,param.V,param.splitCoefs);
		pSemiclassical.evolve(param.T,param.Dt,std::bind(&P::callback,std::ref(param),_1,_2));
	}

    return 0;
}
