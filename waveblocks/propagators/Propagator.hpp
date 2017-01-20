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
#include "../utilities/prettyprint.hpp"
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

namespace print = utilities::prettyprint;
namespace utils = utilities;

/**
 * \brief generic propagator class for Hagedorn Wave Packets
 *
 * \tparam N Number of energy levels
 * \tparam D Dimension of space
 * \tparam MultiIndex Type of multi index used in the basis shape
 * \tparam MDQR Multi-dimensional quadrature rule
 * \tparam Potential_t Type of the Potential to be used
 * \tparam Packet_t Type of the Wavepacket to be propagated
 */

// TODO: consider making Potential_t a parameter of propagate only (but then needs to be passed around a lot, right?)
template <int N, int D, typename MultiIndex, typename MDQR, typename Potential_t, typename Packet_t>
class Propagator {
	
	private:

		Packet_t& wpacket_;
		Potential_t& V_;
		CMatrix<Eigen::Dynamic,Eigen::Dynamic> F_; // TODO: consider renaming from F_ to something more "correct"


	public:

		Propagator(Packet_t& pack, Potential_t& V)
		 : wpacket_(pack)
		 , V_(V)
		{

			/*
                    int size = 0;
                    for (auto& component : packet.components()) {
                        size += component.coefficients().size();
                    }
                    F.resize(size,size);
            */

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

		virtual std::string getName() const = 0; ///< get the name of the current propagator
		virtual void propagate(const real_t) = 0; ///< do the main propagation loop
		virtual void pre_propagate(const real_t) {} ///< pre-propagation work
		virtual void post_propagate(const real_t) {} ///< post-propagation work

	protected:
		
		void buildF() {
			buildF(F_);
		}

		// N>1
        template <int N_LEVELS=N>
        typename std::enable_if<(N_LEVELS>1),void>::type
		buildF(CMatrix<Eigen::Dynamic,Eigen::Dynamic>& M) {

			auto op = [&] (const CMatrix<D,Eigen::Dynamic>& x, const RMatrix<D,1>& q, const dim_t i, const dim_t j)
				{
					///> R: order of the quadrature rule
					///> x: nodal points of quadrature rule (dimension DxR)
					///> q: position (dimension Dx1)
					const dim_t R = x.cols(); ///> Order of the quadrature rule
					CMatrix<1,Eigen::Dynamic> f(1,R); ///> Result f(x,[q_1,...,q_R])

					// #pragma omp parallel for schedule(guided)
					for(int r=0; r<R; ++r) {
						f(0,r) = V_.evaluate_local_remainder_at(
						            utils::Squeeze<D,CMatrix<D,Eigen::Dynamic>>::apply(x,r),
						     	    utils::Squeeze<D,CVector<D>>::apply(complex_t(1,0)*q)
						     	   ) (i,j);
					}
					return f;
				};

			M = innerproducts::VectorInnerProduct<D,MultiIndex,MDQR>::build_matrix(wpacket_,op);

		}

		// N=1
        template <int N_LEVELS=N>
        typename std::enable_if<(N_LEVELS==1),void>::type
		buildF(CMatrix<Eigen::Dynamic,Eigen::Dynamic>& M) {

			auto op = [&] (const CMatrix<D,Eigen::Dynamic>& x, const RMatrix<D,1>& q)
				{
					///> R: order of the quadrature rule
					///> x: nodal points of quadrature rule (dimension DxR)
					///> q: position (dimension Dx1)
					const dim_t R = x.cols(); ///> Order of the quadrature rule
					CMatrix<1,Eigen::Dynamic> f(1,R); ///> Result f(x,[q_1,...,q_R])

					// #pragma omp parallel for schedule(guided)
					for(int r=0; r<R; ++r) {
						f(0,r) = V_.evaluate_local_remainder_at(
						            utils::Squeeze<D,CMatrix<D,Eigen::Dynamic>>::apply(x,r),
						     	    utils::Squeeze<D,CVector<D>>::apply(complex_t(1,0)*q)
						     	   );
					}
					return f;
				};

			M = innerproducts::HomogeneousInnerProduct<D,MultiIndex,MDQR>::build_matrix(wpacket_,op);
		}

		/**
		 * \brief single step with quadratic potential energy U 
		 */
		// Homogeneous
        template <typename P=Packet_t>
        typename std::enable_if<!std::is_same<P,InhomogeneousHaWp<D,MultiIndex>>::value,void>::type
		stepU(const real_t h) {
			// Homogeneous
			// TODO: implement
			///> taylorV[0,1,2] = [V,DV,DDV] = [evaluation,jacobian,hessian]
			// TODO: is this only valid for homogeneous? 
			// const auto& taylorV = V_.get_leading_level().taylor_at(/* q = */ complex_t(1,0) * utils::Squeeze<D,RVector<D>>::apply(params.q()));
			//  --> get leading level for homogeneous packets only
			// TODO: stepU params pass leading level
			stepU_params(h,wpacket_.parameters());
		}
		// Inhomogeneous
        template <typename P=Packet_t>
        typename std::enable_if<std::is_same<P,InhomogeneousHaWp<D,MultiIndex>>::value,void>::type
		stepU(const real_t h) {
			// Inhomogeneous
			// TODO: implement
			///> taylorV[0,1,2] = [V,DV,DDV] = [evaluation,jacobian,hessian]
			// TODO: add inhomogeneous implementation: loop over components
			// TODO: is this only valid for homogeneous? 
			// const auto& taylorV = V_.get_leading_level().taylor_at(/* q = */ complex_t(1,0) * utils::Squeeze<D,RVector<D>>::apply(params.q()));
			//  --> get leading level for homogeneous packets only
			// TODO: stepU params pass leading level
			stepU_params(h,wpacket_.parameters());
		}

		// TODO: is this inlined as a template?? cause it is a normal function!!
		// TODO: make this inline // template <typename T=void>

		/**
		 * \brief single step with kinetic energy operator T
		 */
		// Homogeneous
        template <typename P=Packet_t>
        typename std::enable_if<!std::is_same<P,InhomogeneousHaWp<D,MultiIndex>>::value,void>::type
		stepT(const real_t h) {
			// Homogeneous
			stepT_params(h,wpacket_.parameters());
		}
		// Inhomogeneous
        template <typename P=Packet_t>
        typename std::enable_if<std::is_same<P,InhomogeneousHaWp<D,MultiIndex>>::value,void>::type
		stepT(const real_t h) {
			// Inhomogeneous
			for(auto& comp : wpacket_.components()) {
				stepT_params(h,comp.parameters());
			}
		}
		
		/**
		 * \brief single step with potential energy remainder W
		 */
		void stepW(const real_t h) {
			/*
			// IMPROVEMENT SUGGESTION: change signature of PacketToCoefficients
			//  --> this would be much more readable + self explanatory with the syntax
			buildF();
			complex_t factor(0,-h/(wpacket_.eps()*wpacket_.eps());
			CVector<Eigen::Dynamic> coefs;
			PacketConverter<Packet_t>::PackToCoef(packet,coefs)
			coefs = (factor*F_).exp() * coefs; ///> c = exp(-i*h/eps^2 * F) * c
			PacketConverter<Packet_t>::CoefToPack(coefs,packet)
			*/

			// TODO: Pi' = ?
			// TODO: Pi_i' = ? for all i=0,..,N-1
			buildF(); // TODO: hom/inhom
			CVector<Eigen::Dynamic> coefs = utils::PacketToCoefficients<Packet_t>::to(wpacket_); // get coefficients from packet
			complex_t factor(0,-h/(wpacket_.eps()*wpacket_.eps()));
			coefs = (factor*F_).exp() * coefs; ///> c = exp(-i*h/eps^2 * F) * c
			utils::PacketToCoefficients<Packet_t>::from(coefs,wpacket_); // update packet from coefficients
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
		
		/**
		 * \brief alternately apply T and U, starting with T
		 */
		void splitTU(const std::vector<real_t>& w_T, const std::vector<real_t>& w_U, const real_t dt) {
			assert(w_T.size() == w_U.size() || w_T.size() == w_U.size()+1);
			if(w_T.size()>0) {
				stepT(w_T.at(0)*dt); // do a step of size w_T[0]
				std::vector<real_t> w_T_new(w_T.begin()+1,w_T.end()); // pass on all element but first w_T[1..N-1]
				splitUT(dt,w_U,w_T_new);
			}
		}
		
		/**
		 * \brief alternately apply T and U, starting with U
		 */
		void splitUT(const std::vector<real_t>& w_U, const std::vector<real_t>& w_T, const real_t dt) {
			assert(w_U.size() == w_T.size() || w_U.size() == w_T.size()+1);
			if(w_U.size()>0) {
				stepU(w_U.at(0)*dt); // do a step of size w_U[0]
				std::vector<real_t> w_U_new(w_U.begin()+1,w_U.end()); // pass on all element but first w_U[1..N-1]
				splitTU(dt,w_T,w_U_new);
			}
		}

	
	private:

		/**
		 * \brief update q,Q,S according to operator T
		 */
		void stepT_params(const real_t h, wavepackets::HaWpParamSet<D>& params) {
			// NB: remove mass if not required
			RMatrix<Eigen::Dynamic,Eigen::Dynamic> Minv = RMatrix<D,D>::Identity();
			params.updateq( +h * Minv*params.p() ); ///> q = q + h * M^{-1} * p
			params.updateQ( +h * Minv*params.P() ); ///> Q = Q + h * M^{-1} * P
			params.updateS( +.5f*h * params.p().dot(Minv*params.p()) ); ///> S = S + h/2 * p^T M p
		}

		/**
		 * \brief update p,P,S according to operator U
		 */
		void stepU_params(const real_t h, wavepackets::HaWpParamSet<D>& params) {
			///> taylorV[0,1,2] = [V,DV,DDV] = [evaluation,jacobian,hessian]
			// TODO: what does get_leading_level do???
			//  --> get leading level for homogeneous packets only // TODO: check this for inhom
			//  TODO: consider passing Level& as an additional parameter
			const auto& taylorV = V_.get_leading_level().taylor_at(/* q = */ complex_t(1,0) * utils::Squeeze<D,RVector<D>>::apply(params.q()));
			params.updatep( -h * utils::Unsqueeze<D,RVector<D>>::apply(std::get<1>(taylorV).real()) ); ///> p = p - h * jac(V(q))
			params.updateP( -h * std::get<2>(taylorV)*params.Q() ); ///> P = P - h * hess(V(q)) * Q
			params.updateS( -h * std::get<0>(taylorV) ); ///> S = S - h * V(q)
		}

		// TODO: introduce private functions update_qQS(params,dt) and update_pPS(params,dt)


};


} // namespace propagators
} // namespace waveblocks
