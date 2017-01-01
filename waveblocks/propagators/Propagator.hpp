#pragma once

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



/**
 * \brief Implements the Hagedorn propagator for
 * vector valued wavepackets. Offers a method for
 * time propagation.
 *
 * \tparam N
 * Number of energy levels
 * \tparam D
 * Dimension of space
 * \tparam MultiIndex
 * Type of multi index used in the basis shape
 * \tparam MDQR
 * Multi-dimensional quadrature rule
 */

template <int N, typename MultiIndex, typename MDQR, typename Potential_t>
class Propagator {

	public:

		static const int D = 1;
		// TODO: genrealize to more dimensions, more multiindices, Quadrature Rules
		using Packet_t = wavepackets::ScalarHaWp<D,MultiIndex>; // What's up with ScalarHaWp

		Propagator(Packet_t& pack, Potential_t& V)
		 : t_(0)
		 , wpacket_(pack)
		 , V_(V)
		{
			// unsigned size = 0;
			// size = 1; // Scalar Case
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
				          << "\t\tTime: " << std::setw(10) << std::fixed << std::setprecision(4) << std::right << t_ << std::flush;
			}
			std::cout << "\n\n";
			post_propagate(Dt);

			mywriter.poststructuring();

		}

		virtual void writeData(io::hdf5writer<D>& writer) const {
			// Compute energies
			real_t ekin = observables::kinetic_energy<D,MultiIndex>(wpacket_);
			real_t epot = 0; // observables::potential_energy<Packet_t,D,MultiIndex,MDQR>(wpacket_,V_);

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

			M = HomogeneousInnerProduct<D,MultiIndex,MDQR>::build_matrix(wpacket_,op);

		}

		void stepU(real_t h) {
			///> taylorV[0,1,2] = [V,DV,DDV] = [evaluation,jacobian,hessian]
			auto& params = wpacket_.parameters();
			// TODO: what does get_leading_level do???
			const auto& taylorV = V_.get_leading_level().taylor_at(/* q = */ complex_t(1,0) * Squeeze<D,RVector<D>>::apply(params.q()));
			params.updatep( -h * Unsqueeze<D,RVector<D>>::apply(std::get<1>(taylorV).real()) ); ///> p = p - h * jac(V(q))
			params.updateP( -h * std::get<2>(taylorV)*params.Q() ); ///> P = P - h * hess(V(q)) * Q
			params.updateS( -h * std::get<0>(taylorV) ); ///> S = S - h * V(q)
		}
		void stepT(real_t h) {
			// TODO: add inverse mass Minv
			real_t Minv = 1.;
			auto& params = wpacket_.parameters();
			params.updateq( +h * Minv*params.p() ); ///> q = q + h * M^{-1} * p
			params.updateQ( +h * Minv*params.P() ); ///> Q = Q + h * M^{-1} * P
			params.updateS( +.5*h * params.p().dot(Minv*params.p()) ); ///> S = S + h/2 * p^T M p
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

		void intSplit(real_t Dt, unsigned M) {
			real_t dt = Dt/M;
			for(unsigned m=0; m<M; ++m) {
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


} // namespace propagators
} // namespace waveblocks
