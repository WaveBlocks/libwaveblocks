#include <iostream>
#include <iomanip>

#include "../innerproducts/homogeneous_inner_product.hpp"
#include "../innerproducts/vector_inner_product.hpp"
#include "../potentials/bases.hpp"
#include "../potentials/potentials.hpp"
#include "../types.hpp"
#include "../utilities/adaptors.hpp"
#include "../utilities/squeeze.hpp"
#include "waveblocks/observables/energy.hpp"
#include "waveblocks/io/hdf5writer.hpp"

/** \file
 * \brief Abstract Propagator Class
 * 
 * \tparam Type of wave packet
 * 
 */

namespace waveblocks {
namespace propagators {


using utilities::Squeeze;
using utilities::PacketToCoefficients;
using utilities::Unsqueeze;
using wavepackets::HaWpParamSet;

template <typename Potential_t>
class Propagator {

	public:


		// TODO: genrealize to more dimensions, more multiindices, Quadrature Rules
		static const int D = 1;
		using MultiIndex = wavepackets::shapes::TinyMultiIndex<unsigned short, D>;
		using QR = innerproducts::GaussHermiteQR</*K+4*/132>;
		using Packet_t = wavepackets::ScalarHaWp<D,MultiIndex>; // What's up with ScalarHaWp

		Propagator(Packet_t& pack, Potential_t& V)
		 : t_(0)
		 , wpacket_(pack)
		 , V_(V)
		{
			unsigned size = 0;
			size = 1; // Scalar Case
			// for(auto& c : wpacket_.components()) {
			// 	size += c.coefficients.size();
			// }
			// F_.resize(size,size);
		}

		// TODO: make CRTP instead
		void simulate(const real_t T, const real_t Dt, const std::string outfilename="data") {

			// TODO: introduce an extra callback function that is called in every iteration
			// TODO: provide public functions setupfile, savefile

			//////////////////////////////////////////////////////////////////////////////
			// TODO: wrap in its own function
			
			// Preparing the file and I/O writer
			io::hdf5writer<D> mywriter(outfilename+".hdf5");
			mywriter.set_write_norm(true);
			mywriter.set_write_energies(true);
			mywriter.prestructuring<MultiIndex>(wpacket_,Dt);

			writeData(mywriter);
			//////////////////////////////////////////////////////////////////////////////

			unsigned M = std::round(T/Dt);
			pre_propagate(Dt);
			std::cout << "\n\n";
			for(unsigned m=0; m<M; ++m) {
				propagate(Dt);
				t_ += Dt;
				writeData(mywriter); // TODO: if condition
				std::cout << "\rProgress: " << std::setw(6) << std::fixed << std::setprecision(1) << std::right << (100.*(m+1))/M << "%"
				          << "\t\tTime: " << std::setw(10) << std::right << t_ << std::flush;
			}
			std::cout << "\n\n";
			post_propagate(Dt);

			mywriter.poststructuring();

		}

		virtual void writeData(io::hdf5writer<D>& writer) const {
			// Compute energies
			real_t ekin = observables::kinetic_energy<D,MultiIndex>(wpacket_);
			real_t epot = 0; // observables::potential_energy<Packet_t,D,MultiIndex,QR>(wpacket_,V_);

			writer.store_packet(wpacket_);
			writer.store_norm(wpacket_);
			writer.store_energies(epot,ekin);
		}

		virtual void propagate(real_t) = 0;
		virtual void pre_propagate(real_t) {}
		virtual void post_propagate(real_t) {}

	protected:
		
		void buildF() {
			buildF(F_);
		}

		void buildF(CMatrix<Eigen::Dynamic,Eigen::Dynamic>& M) {

			// TODO: generalize to multidimensional case (efficiently)

			auto op = [this] (const CMatrix<D,Eigen::Dynamic>& x, const RMatrix<D,1>& q)
				{
					///> R: order of the quadrature rule
					///> x: nodal points of quadrature rule (dimension DxR)
					///> q: position (dimension Dx1)
					const dim_t R = x.cols(); ///> Order of the quadrature rule
					CMatrix<1,Eigen::Dynamic> f(1,R); ///> Result f(x,[q_1,...,q_R])

					#pragma omp parallel for schedule(guided)
					for(int r=0; r<R; ++r) {
						f(0,r) = V_.evaluate_local_remainder_at(
						            Squeeze<D,CMatrix<D,Eigen::Dynamic>>::apply(x,r),
						     	    Squeeze<D,CVector<D>>::apply(complex_t(1,0)*q)
						     	   );
					}
					return f;
				};

			M = HomogeneousInnerProduct<D,MultiIndex,QR>::build_matrix(wpacket_,op);

		}

		void stepU(real_t h) {
			///> taylorV[0,1,2] = [V,DV,DDV] = [evaluation,jacobian,hessian]
			auto& params = wpacket_.parameters();
			// TODO: what does get_leading_level do???
			const auto& taylorV = V_.get_leading_level().taylor_at(/* q = */ complex_t(1,0) * Squeeze<D,RVector<D>>::apply(params.q()));
			params.updatep( -h * Unsqueeze<D,RVector<D>>::apply(std::get<1>(taylorV).real()) );
			params.updateP( -h * std::get<2>(taylorV)*params.Q() );
			params.updateS( -h * std::get<0>(taylorV) );
		}
		void stepT(real_t h) {
			// TODO: add inverse mass Minv
			real_t Minv = 1.;
			auto& params = wpacket_.parameters();
			params.updateq(h*Minv*params.p()); ///> q = q + h * M^{-1} * p
			params.updateQ(h*Minv*params.P()); ///> Q = Q + h * M^{-1} * P
			params.updateS(.5*h*params.p().dot(Minv*params.p())); ///> S = S + h/2 * p^T M p
		}

		void stepW(real_t h) {
			/*
			// IMPROVEMENT SUGGESTION: change signature of PacketToCoefficients
			//  --> this would be much more readable + self explanatory with the syntax
			buildF();
			complex_t factor(0,-h/(wpacket_.eps()*wpacket_.eps());
			CVector<Eigen::Dynamic> coefs;
			PacketConverter<Packet_t>::PackToCoef(packet,coefs)
			coefs = (factor*F_).exp() * coefs; ///> c = exp(-i*h/eps^2 * F) * c
			PacketConverter<Packet_t>::CoefToPack(packet,coefs)
			*/

			buildF();
			CVector<Eigen::Dynamic> coefs = PacketToCoefficients<Packet_t>::to(wpacket_); // get coefficients from packet
			complex_t factor(0,-h/(wpacket_.eps()*wpacket_.eps()));
			coefs = (factor*F_).exp() * coefs; ///> c = exp(-i*h/eps^2 * F) * c
			PacketToCoefficients<Packet_t>::from(coefs,wpacket_); // update packet from coefficients
		}

		void intSplit(real_t Dt, unsigned N) {
			real_t dt = Dt/N;
			for(unsigned n=0; n<N; ++n) {
				// TODO: do this properly!! With weights!!
				// alternating templates!!
				stepU(dt);
				stepT(dt);
			}
		}

		real_t t_;
		Packet_t& wpacket_;
		Potential_t& V_;
		CMatrix<Eigen::Dynamic,Eigen::Dynamic> F_; // TODO: consider renaming from F_ to something more "correct"

};


/** \file
 * \brief Implements the Hagedorn Propagator
 */
template <typename Potential_t>
class HagedornPropagator : public Propagator<Potential_t> /* TODO: implement SplittingParameters , public SplittingParameters */ /* TODO: make CRTP */ {

	public:

		using Propagator<Potential_t>::Propagator; // inherit constructor

		void propagate(real_t Dt) override {

			// Hagedorn
			this->stepT(Dt/2);
			this->stepU(Dt);
			this->stepW(Dt);
			this->stepT(Dt/2);


			// Semiclassical
			// TODO: compute N
			// TODO: make sure that N is even, otherwise N/2 + N/2 != N
			int N = 2;
			this->intSplit(Dt/2,N/2);
			this->stepW(Dt);
			this->intSplit(Dt/2,N/2);


		// 	// Magnus
		// 	// TODO: introduce additional temporary variables A1_, A2_ for derived class
		// 	// TODO: consider saving sqrt(3) as a private class constant
		// 	real_t h1 = (3-std::sqrt(3))/6;
		// 	real_t h2 = std::sqrt(3)/6;
		// 	int N1 = 2; // TODO: put right formula
		// 	int N2 = 2; // TODO: put right formula
		// 	intSplit(h1,N1);
		// 	buildF(A1_);
		// 	A1_ *= complex_t(0,-Dt/(wpacket_.eps()*wpacket_.eps()));
		// 	intSplit(h2,N2);
		// 	buildF(A2_);
		// 	A2_ *= complex_t(0,-Dt/(wpacket_.eps()*wpacket_.eps()));
		// 	F_ = .5*(A1_+A2_) + std::sqrt(3)/12*(A1_*A2_-A2_*A1_);
		// 	// TODO: matrix exponential c = exp(F_) * c
		// 	intSplit(h1,N1);

		// 	// Pre764
		// 	// TODO: introduce additional coefficient vectors Z_, Y_, alpha_, beta_ for derived class
		// 	const int N = 1 + std::floor(std::sqrt(Dt*std::pow(wpacket_.eps(),-.75)));
		// 	// TODO: v, k could even be class constants v_, k_
		// 	const int v = 6;
		// 	const int k = 4;

		// 	// PRE: TODO: move to separate function, not in propagate
		// 	for(unsigned j=0; j<v; ++j) {
		// 		intSplit(-Z_[j]*Dt,N);
		// 		propW(-Y_[j]*Dt);
		// 	}
		// 	// PROPAGATE:
		// 	for(unsigned j=0; j<k; ++j) {
		// 		propW(alpha_[j]*Dt);
		// 		intSplit(beta_[j]*Dt);
		// 	}
		// 	// POST: TODO: move to separate function, not in propagate
		// 	for(unsigned j=0; j<v; ++j) {
		// 		propW(Y_[j]*Dt);
		// 		intSplit(Z_[j]*Dt,N);
		// 	}


		// 	// McL42
		// 	intSplit(A_[0]*Dt,N);
		// 	propW(B_[0]*Dt);
		// 	intSplit(A_[1]*Dt,N);
		// 	propW(B_[1]*Dt);
		// 	intSplit(A_[2]*Dt,N);

		}

};


} // namespace propagators
} // namespace waveblocks
