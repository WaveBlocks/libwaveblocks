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


class Parameters_Torsional_1D {

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
					return 1. - std::cos(x.real());
				}
				
				inline Taylor::jacobian_evaluation_type evalJ(const Taylor::argument_type& x) const {
					return std::sin(x.real());
				}
				
				inline Taylor::hessian_evaluation_type evalH(const Taylor::argument_type& x) const {
					return std::cos(x.real());
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

		// basis shapes
		wavepackets::shapes::ShapeEnumerator<D, MultiIndex> enumerator;
		wavepackets::shapes::ShapeEnum<D, MultiIndex> shape_enum;

		Parameters_Torsional_1D()
		 : sigma_x(1.)
		 , T(10)
		 , Dt(0.001)
		 , eps(0.01)
		 , Q(CMatrix<D,D>::Identity())
		 , P(complex_t(0,1) * CMatrix<D,D>::Identity())
		 , q(RVector<D>::Ones())
		 , p(RVector<D>::Zero())
		 , S(0.)
		 , param_set(q,p,Q,P,S)
		 , coeffs(Coefficients::Zero(std::pow(K,D),1))
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

	using P = Parameters_Torsional_1D; // all the parameters for the simulation are contained in this class

	{ // Semiclassical
		Parameters_Torsional_1D param_Semiclassical;
		propagators::SemiclassicalPropagator<P::N,P::D,P::MultiIndex,P::QR,P::Potential,P::Packet_t,propagators::SplitCoefs<1,1>> propagator_Semiclassical(param_Semiclassical.packet,param_Semiclassical.V,propagators::splitting_parameters::coefLT);
		propagator_Semiclassical.evolve(param_Semiclassical.T,param_Semiclassical.Dt);
	}
	
	{ // Semiclassical
		Parameters_Torsional_1D param_Semiclassical;
		propagators::SemiclassicalPropagator<P::N,P::D,P::MultiIndex,P::QR,P::Potential,P::Packet_t,propagators::SplitCoefs<2,2>> propagator_Semiclassical(param_Semiclassical.packet,param_Semiclassical.V,propagators::splitting_parameters::coefS2);
		propagator_Semiclassical.evolve(param_Semiclassical.T,param_Semiclassical.Dt);
	}
	
	{ // Semiclassical
		Parameters_Torsional_1D param_Semiclassical;
		propagators::SemiclassicalPropagator<P::N,P::D,P::MultiIndex,P::QR,P::Potential,P::Packet_t,propagators::SplitCoefs<4,4>> propagator_Semiclassical(param_Semiclassical.packet,param_Semiclassical.V,propagators::splitting_parameters::coefY4);
		propagator_Semiclassical.evolve(param_Semiclassical.T,param_Semiclassical.Dt);
	}
	
	{ // Semiclassical
		Parameters_Torsional_1D param_Semiclassical;
		propagators::SemiclassicalPropagator<P::N,P::D,P::MultiIndex,P::QR,P::Potential,P::Packet_t,propagators::SplitCoefs<7,7>> propagator_Semiclassical(param_Semiclassical.packet,param_Semiclassical.V,propagators::splitting_parameters::coefPRKS6);
		propagator_Semiclassical.evolve(param_Semiclassical.T,param_Semiclassical.Dt);
	}
	
	{ // Semiclassical
		Parameters_Torsional_1D param_Semiclassical;
		propagators::SemiclassicalPropagator<P::N,P::D,P::MultiIndex,P::QR,P::Potential,P::Packet_t,propagators::SplitCoefs<8,8>> propagator_Semiclassical(param_Semiclassical.packet,param_Semiclassical.V,propagators::splitting_parameters::coefY61);
		propagator_Semiclassical.evolve(param_Semiclassical.T,param_Semiclassical.Dt);
	}
	
	{ // Semiclassical
		Parameters_Torsional_1D param_Semiclassical;
		propagators::SemiclassicalPropagator<P::N,P::D,P::MultiIndex,P::QR,P::Potential,P::Packet_t,propagators::SplitCoefs<10,10>> propagator_Semiclassical(param_Semiclassical.packet,param_Semiclassical.V,propagators::splitting_parameters::coefKL6);
		propagator_Semiclassical.evolve(param_Semiclassical.T,param_Semiclassical.Dt);
	}
	
	{ // Semiclassical
		Parameters_Torsional_1D param_Semiclassical;
		propagators::SemiclassicalPropagator<P::N,P::D,P::MultiIndex,P::QR,P::Potential,P::Packet_t,propagators::SplitCoefs<15,15>> propagator_Semiclassical(param_Semiclassical.packet,param_Semiclassical.V,propagators::splitting_parameters::coefBM63);
		propagator_Semiclassical.evolve(param_Semiclassical.T,param_Semiclassical.Dt);
	}
	
	{ // Semiclassical
		Parameters_Torsional_1D param_Semiclassical;
		propagators::SemiclassicalPropagator<P::N,P::D,P::MultiIndex,P::QR,P::Potential,P::Packet_t,propagators::SplitCoefs<18,18>> propagator_Semiclassical(param_Semiclassical.packet,param_Semiclassical.V,propagators::splitting_parameters::coefKL8);
		propagator_Semiclassical.evolve(param_Semiclassical.T,param_Semiclassical.Dt);
	}
	
	{ // Semiclassical
		Parameters_Torsional_1D param_Semiclassical;
		propagators::SemiclassicalPropagator<P::N,P::D,P::MultiIndex,P::QR,P::Potential,P::Packet_t,propagators::SplitCoefs<34,34>> propagator_Semiclassical(param_Semiclassical.packet,param_Semiclassical.V,propagators::splitting_parameters::coefKL10);
		propagator_Semiclassical.evolve(param_Semiclassical.T,param_Semiclassical.Dt);
	}
	
    return 0;
}
