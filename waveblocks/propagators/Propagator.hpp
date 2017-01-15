#pragma once

#include <iostream>
#include <iomanip>
#include <functional>
#include <type_traits>

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

using innerproducts::HomogeneousInnerProduct;
using innerproducts::VectorInnerProduct;


namespace print {
	static unsigned WIDTH = 60;

	template <typename T>
	void pair(const std::string&& s, const T v, const std::string&& s0="\n") {
		std::streamsize default_precision = std::cout.precision();
		std::cout << s0 << "\t"
			<< std::setw(WIDTH/2) << std::left << s
			<< std::setw(WIDTH/2) << std::right << std::fixed << std::setprecision(4) << v
			<< std::scientific << std::setprecision(default_precision) << std::flush;
	}
	
	inline void title(const std::string&& s, const char c='-') {
		int str_length = s.length();
		int w1 = (WIDTH-str_length-2)/2;
		int w2 = (WIDTH-str_length-2) - w1;
		std::cout
			<< "\n\t" << std::string(WIDTH,c)
			<< "\n\t" << c << std::string(w1,' ') << s << std::string(w2,' ') << c
			<< "\n\t" << std::string(WIDTH,c);
	}

	inline void separator(const char c='-') {
		std::cout << "\n\t" << std::string(WIDTH,c) << std::flush;
	}

}



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

template <int N, int D, typename MultiIndex, typename MDQR, typename Potential_t, typename Packet_t>
class Propagator {

	public:

		Propagator(Packet_t& pack, Potential_t& V)
		 : wpacket_(pack)
		 , V_(V)
		{
			// unsigned size = 0;
			// size = 1; // Scalar Case
			// for(auto& c : wpacket_.components()) {
			// 	size += c.coefficients.size();
			// }
			// F_.resize(size,size);
			const bool scalar = std::is_same<Packet_t,wavepackets::ScalarHaWp<D,MultiIndex>>::value;
			static_assert(scalar == (N==1),"Scalar wave packets must have N==1");
		}

		// TODO: make CRTP instead
		/**
		 * \param callback An optional callback function that is called before doing any time step and at the end of the propagation
		 *  The callback function must take two arguments: the index of the current iteration (unsigned integer) and the current time (real_t)
		 */
		void evolve(const real_t T, const real_t Dt,
				const std::function<void(unsigned, real_t)> callback = [](unsigned i, real_t t) { (void)i; (void)t; }) {

			if(T<0) return;

			real_t t = 0;

			// TODO: ifdef verbose
			{
				const bool scalar = std::is_same<Packet_t,wavepackets::ScalarHaWp<D,MultiIndex>>::value;
				std::cout << "\n\n";
				print::title(getName() + " Propagator");
				print::pair((scalar ? "Scalar" : "Vectorial") + std::string(" Wave Packet"),"");
				// TODO: homogeneous / inhomogeneous?
				print::pair("D (number of dimensions)",D);
				print::pair("N (number of energy levels)",N);
				print::pair("Quadrature Rule","TODO: which quadrature rule");
				print::pair("Order of the scheme","TODO: order of method");
				print::separator();
				print::pair("Time Span T",T);
				print::pair("Stepsize Dt",Dt);
				print::separator();
				std::cout << "\n";
			}

			unsigned M = std::round(T/Dt);
			pre_propagate(Dt);

			for(unsigned m=0; m<M; ++m) {
				t += Dt;
				callback(m,t);
				print::pair("Time t",t,"\r");
				propagate(Dt);
			}
			callback(M,t);
			print::pair("","COMPLETE","\r");
			print::separator();

			post_propagate(Dt);

		}

		virtual std::string getName() const = 0;
		virtual void propagate(const real_t) = 0;
		virtual void pre_propagate(const real_t) {}
		virtual void post_propagate(const real_t) {}

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

					// #pragma omp parallel for schedule(guided)
					for(int r=0; r<R; ++r) {
						f(0,r) = V_.evaluate_local_remainder_at(
						            Squeeze<D,CMatrix<D,Eigen::Dynamic>>::apply(x,r),
						     	    Squeeze<D,CVector<D>>::apply(complex_t(1,0)*q)
						     	   );
					}
					return f;
				};

			using Innerproduct_t = typename std::conditional<
			                           std::is_same<Packet_t,wavepackets::ScalarHaWp<D,MultiIndex>>::value, // scalar wavepacket?
			                           HomogeneousInnerProduct<D,MultiIndex,MDQR>, // 1 dimensional
			                           VectorInnerProduct<D,MultiIndex,MDQR>       // N dimensional
			                        >::type;
			M = Innerproduct_t::build_matrix(wpacket_,op);

		}

		void stepU(const real_t h) {
			///> taylorV[0,1,2] = [V,DV,DDV] = [evaluation,jacobian,hessian]
			auto& params = wpacket_.parameters();
			// TODO: what does get_leading_level do???
			// TODO: add inhomogeneous implementation: loop over components
			//  --> get leading level for homogeneous packets only
			const auto& taylorV = V_.get_leading_level().taylor_at(/* q = */ complex_t(1,0) * Squeeze<D,RVector<D>>::apply(params.q()));
			params.updatep( -h * Unsqueeze<D,RVector<D>>::apply(std::get<1>(taylorV).real()) ); ///> p = p - h * jac(V(q))
			params.updateP( -h * std::get<2>(taylorV)*params.Q() ); ///> P = P - h * hess(V(q)) * Q
			params.updateS( -h * std::get<0>(taylorV) ); ///> S = S - h * V(q)
		}
		void stepT(const real_t h) {
			// TODO: add inhomogeneous implementation: loop over components
			real_t Minv = 1.f; // inverse mass
			auto& params = wpacket_.parameters();
			params.updateq( +h * Minv*params.p() ); ///> q = q + h * M^{-1} * p
			params.updateQ( +h * Minv*params.P() ); ///> Q = Q + h * M^{-1} * P
			params.updateS( +.5*h * params.p().dot(Minv*params.p()) ); ///> S = S + h/2 * p^T M p
		}

		void stepW(const real_t h) {
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

			// TODO: Pi' = ?
			// TODO: Pi_i' = ? for all i=0,..,N-1
			buildF(); // TODO: hom/inhom
			CVector<Eigen::Dynamic> coefs = PacketToCoefficients<Packet_t>::to(wpacket_); // get coefficients from packet
			complex_t factor(0,-h/(wpacket_.eps()*wpacket_.eps()));
			coefs = (factor*F_).exp() * coefs; ///> c = exp(-i*h/eps^2 * F) * c
			PacketToCoefficients<Packet_t>::from(coefs,wpacket_); // update packet from coefficients
		}

		void intSplit(const real_t Dt, const unsigned M) {
			real_t dt = Dt/M;
			const std::vector<real_t> a = { 1,2,3 };
			const std::vector<real_t> b = { 4,5,6 };
			for(unsigned m=0; m<M; ++m) {
				// TODO: do this properly!! With weights!!
				// alternating templates!!
				splitTU(a,b,dt);
			}
		}
		
		void splitTU(const std::vector<real_t>& w_T, const std::vector<real_t>& w_U, const real_t dt) {
			assert(w_T.size() == w_U.size() || w_T.size() == w_U.size()+1);
			if(w_T.size()>0) {
				stepT(w_T.at(0)*dt); // do a step of size w_T[0]
				std::vector<real_t> w_T_new(w_T.begin()+1,w_T.end()); // pass on all element but first w_T[1..N-1]
				splitUT(dt,w_U,w_T_new);
			}
		}
		
		void splitUT(const std::vector<real_t>& w_U, const std::vector<real_t>& w_T, const real_t dt) {
			assert(w_U.size() == w_T.size() || w_U.size() == w_T.size()+1);
			if(w_U.size()>0) {
				stepU(w_U.at(0)*dt); // do a step of size w_U[0]
				std::vector<real_t> w_U_new(w_U.begin()+1,w_U.end()); // pass on all element but first w_U[1..N-1]
				splitTU(dt,w_T,w_U_new);
			}
		}


		Packet_t& wpacket_;
		Potential_t& V_;
		CMatrix<Eigen::Dynamic,Eigen::Dynamic> F_; // TODO: consider renaming from F_ to something more "correct"

};


} // namespace propagators
} // namespace waveblocks
