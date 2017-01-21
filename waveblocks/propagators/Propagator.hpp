#pragma once

#include <array>
#include <functional>
#include <iomanip>
#include <iostream>
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

		/**
		 * \brief Propagator constructor for scalar wave packets
		 *
		 * Constructs a propagator for N=1 energy levels
		 *
		 * \param pack wave packet to be propagated
		 * \param V potential to use for propagation
		 */
		template <typename U=Packet_t, typename std::enable_if<std::is_same<U,ScalarHaWp<D,MultiIndex>>::value,int>::type = 0>
		Propagator(Packet_t& pack, Potential_t& V)
		 : wpacket_(pack)
		 , V_(V)
		{
			static_assert(N==1,"Scalar wave packets must have N==1");
		}

		/**
		 * \brief Propagator constructor for multi-level wave packets
		 *
		 * Constructs a propagator for N>1 energy levels, resizes the matrix F_ to the correct size for future computations
		 *
		 * \param pack wave packet to be propagated
		 * \param V potential to use for propagation
		 */
		template <typename U=Packet_t, typename std::enable_if<!std::is_same<U,ScalarHaWp<D,MultiIndex>>::value,int>::type = 0>
		Propagator(Packet_t& pack, Potential_t& V)
		 : wpacket_(pack)
		 , V_(V)
		{
			static_assert(N>1,"Multi-Level wave packets must have N>1");
			unsigned size = 0;
			for(auto& comp : wpacket_.components()) {
				size += comp.coefficients.size();
			}
			F_.resize(size,size);
		}



/////////////////////////////////////////////////////////////////////////////////
// Quadratic Potential Energy Operator U
/////////////////////////////////////////////////////////////////////////////////

		/**
		 * \brief function for time evolution
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



/////////////////////////////////////////////////////////////////////////////////
// Functions to be provided by derived classes (CRTP)
/////////////////////////////////////////////////////////////////////////////////

	// TODO: make CRTP instead
	protected:
		virtual std::string getName() const = 0; ///< get the name of the current propagator
		virtual void propagate(const real_t) = 0; ///< do the main propagation loop
		virtual void pre_propagate(const real_t) {} ///< pre-propagation work
		virtual void post_propagate(const real_t) {} ///< post-propagation work



		
		// TODO: is this inlined as a template?? cause it is a normal function!!
		// TODO: make this inline // template <typename T=void>

/////////////////////////////////////////////////////////////////////////////////
// Kinetic Energy Operator T
/////////////////////////////////////////////////////////////////////////////////

	protected:
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
	
	private:
		/**
		 * \brief update q,Q,S according to kinetic operator T
		 */
		void stepT_params(const real_t h, wavepackets::HaWpParamSet<D>& params) {
			// NB: remove mass if not required
			RMatrix<Eigen::Dynamic,Eigen::Dynamic> Minv = RMatrix<D,D>::Identity();
			params.updateq( +h * Minv*params.p() ); ///< q = q + h * M^{-1} * p
			params.updateQ( +h * Minv*params.P() ); ///< Q = Q + h * M^{-1} * P
			params.updateS( +.5f*h * params.p().dot(Minv*params.p()) ); ///< S = S + h/2 * p^T M p
		}



/////////////////////////////////////////////////////////////////////////////////
// Quadratic Potential Energy Operator U
/////////////////////////////////////////////////////////////////////////////////
		
	protected:
		/**
		 * \brief single step with quadratic potential energy U 
		 */
		// Homogeneous
        template <typename P=Packet_t>
        typename std::enable_if<!std::is_same<P,InhomogeneousHaWp<D,MultiIndex>>::value,void>::type
		stepU(const real_t h) {
			// Homogeneous
			// taylorV[0,1,2] = [V,DV,DDV] = [evaluation,jacobian,hessian]
			auto& params = wpacket_.parameters();
			const auto& taylorV = V_.get_leading_level().taylor_at(/* q = */ complex_t(1,0) * utils::Squeeze<D,RVector<D>>::apply(params.q()));
			stepU_params(h,params,std::get<0>(taylorV),std::get<1>(taylorV),std::get<2>(taylorV));
		}
		// Inhomogeneous
        template <typename P=Packet_t>
        typename std::enable_if<std::is_same<P,InhomogeneousHaWp<D,MultiIndex>>::value,void>::type
		stepU(const real_t h) {
			// Inhomogeneous
			// taylorV[0,1,2] = [V,DV,DDV] = [evaluation,jacobian,hessian]
			// NB: inefficient, since all levels 0 <= i < N are evaluated for each single level
			int i=0;
			for(auto& comp : wpacket_.components()) {
				auto& params = comp.parameters();
				const auto& taylorV = V_.get_leading_level().taylor_at(/* q = */ complex_t(1,0) * utils::Squeeze<D,RVector<D>>::apply(params.q()));
				stepU_params(h,params,std::get<0>(taylorV)[i],std::get<1>(taylorV)[i],std::get<2>(taylorV)[i]);
				++i;
			}
		}
	
	private:
		/**
		 * \brief update p,P,S according to quadratic potential operator U
		 */
		template <typename V_T, typename DV_T, typename DDV_T>
		void stepU_params(const real_t h, wavepackets::HaWpParamSet<D>& params, const V_T& potential, const DV_T& jacobian, const DDV_T& hessian) {
			params.updatep( -h * utils::Unsqueeze<D,RVector<D>>::apply(jacobian.real()) ); ///< p = p - h * jac(V(q))
			params.updateP( -h * hessian*params.Q() ); ///< P = P - h * hess(V(q)) * Q
			params.updateS( -h * potential ); ///< S = S - h * V(q)
		}



/////////////////////////////////////////////////////////////////////////////////
// Non-Quadratic Potential Energy Operator W
/////////////////////////////////////////////////////////////////////////////////

	protected:
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
			coefs = (factor*F_).exp() * coefs; ///< c = exp(-i*h/eps^2 * F) * c
			PacketConverter<Packet_t>::CoefToPack(coefs,packet)
			*/

			// TODO: Pi' = ?
			// TODO: Pi_i' = ? for all i=0,..,N-1
			buildF(); // TODO: hom/inhom
			CVector<Eigen::Dynamic> coefs = utils::PacketToCoefficients<Packet_t>::to(wpacket_); // get coefficients from packet
			complex_t factor(0,-h/(wpacket_.eps()*wpacket_.eps()));
			coefs = (factor*F_).exp() * coefs; ///< c = exp(-i*h/eps^2 * F) * c
			utils::PacketToCoefficients<Packet_t>::from(coefs,wpacket_); // update packet from coefficients
		}

		void buildF() {
			buildF(F_);
		}

		// N>1
        template <int N_LEVELS=N>
        typename std::enable_if<(N_LEVELS>1),void>::type
		buildF(CMatrix<Eigen::Dynamic,Eigen::Dynamic>& M) {

			auto op = [&] (const CMatrix<D,Eigen::Dynamic>& x, const RMatrix<D,1>& q, const dim_t i, const dim_t j)
				{
					///< R: order of the quadrature rule
					///< x: nodal points of quadrature rule (dimension DxR)
					///< q: position (dimension Dx1)
					const dim_t R = x.cols(); ///< Order of the quadrature rule
					CMatrix<1,Eigen::Dynamic> f(1,R); ///< Result f(x,[q_1,...,q_R])

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
					///< R: order of the quadrature rule
					///< x: nodal points of quadrature rule (dimension DxR)
					///< q: position (dimension Dx1)
					const dim_t R = x.cols(); ///< Order of the quadrature rule
					CMatrix<1,Eigen::Dynamic> f(1,R); ///< Result f(x,[q_1,...,q_R])

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



/////////////////////////////////////////////////////////////////////////////////
// Int-Split
/////////////////////////////////////////////////////////////////////////////////

	protected:
		/**
		 * \brief split the timestep of size Dt into M smaller timesteps
		 */
		template <typename PROPAGATOR, const std::size_t S_A, const std::size_t S_B>
		void intSplit(PROPAGATOR& prop, const real_t Dt, const unsigned M, const std::pair<std::array<real_t,S_A>,std::array<real_t,S_B>>& coefs) {
			// TODO: consider providing intSplitTU, intSplitUT
			static_assert(S_A==S_B || S_A==S_B+1,"invalid combination of coefficient array sizes");
			const real_t dt = Dt/M;
			for(unsigned m=0; m<M; ++m) {
				TU<PROPAGATOR,0,0,S_A,S_B>::split(prop,coefs.first,coefs.second,dt);
			}
		}
		
		template <typename PROPAGATOR, std::size_t N_T, std::size_t N_U, std::size_t S_T, std::size_t S_U> struct TU;
		template <typename PROPAGATOR, std::size_t N_U, std::size_t N_T, std::size_t S_U, std::size_t S_T> struct UT;

		// termination condition
		template <typename PROPAGATOR, std::size_t S_T, std::size_t S_U>
		struct TU<PROPAGATOR,S_T,S_U,S_T,S_U> {
			static void split(PROPAGATOR& prop, const std::array<real_t,S_T>& coefT, const std::array<real_t,S_U>& coefU, const real_t dt) {
				(void) prop; (void) coefT; (void) coefU; (void) dt;
			}
		};

		// termination condition
		template <typename PROPAGATOR, std::size_t S_U, std::size_t S_T>
		struct UT<PROPAGATOR,S_U,S_T,S_U,S_T> {
			static void split(PROPAGATOR& prop, const std::array<real_t,S_U>& coefU, const std::array<real_t,S_T>& coefT, const real_t dt) {
				(void) prop; (void) coefT; (void) coefU; (void) dt;
			}
		};


		/**
		 * \brief alternately apply T and U, starting with T
		 */
		template <typename PROPAGATOR, std::size_t N_T, std::size_t N_U, std::size_t S_T, std::size_t S_U>
		struct TU {
			static void split(PROPAGATOR& prop, const std::array<real_t,S_T>& coefT, const std::array<real_t,S_U>& coefU, const real_t dt) {

				static_assert(N_T==N_U || N_T==N_U-1,"ERROR: invalid alternating index pair for TU");
				std::cout << "\nTU - (" << N_T << "," << N_U << ")";

				// Apply Steps
				/* T     */ prop.stepT(coefT.at(N_T));
				/* UT... */ UT<PROPAGATOR,N_U,N_T+1,S_U,S_T>::split(prop,coefU,coefT,dt);

			}
		};

		/**
		 * \brief alternately apply T and U, starting with U
		 */
		template <typename PROPAGATOR, std::size_t N_U, std::size_t N_T, std::size_t S_U, std::size_t S_T>
		struct UT {
			static void split(PROPAGATOR& prop, const std::array<real_t,S_U>& coefU, const std::array<real_t,S_T>& coefT, const real_t dt) {

				static_assert(N_U==N_T || N_U==N_T-1,"ERROR: invalid alternating index pair for UT");
				std::cout << "\nUT - (" << N_U << "," << N_T << ")";

				// Apply Steps
				/* U     */ prop.stepU(coefU.at(N_U));
				/* TU... */ TU<PROPAGATOR,N_T,N_U+1,S_T,S_U>::split(prop,coefT,coefU,dt);

			}
		};

};
	

} // namespace propagators
} // namespace waveblocks
