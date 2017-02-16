#pragma once

#include <array>
#include <functional>
#include <iomanip>
#include <iostream>
#include <string>
#include <type_traits>

#include "../innerproducts/homogeneous_inner_product.hpp"
#include "../innerproducts/vector_inner_product.hpp"
#include "../potentials/bases.hpp"
#include "../potentials/potentials.hpp"
#include "../types.hpp"
#include "../utilities/adaptors.hpp"
#include "../utilities/prettyprint.hpp"
#include "../utilities/squeeze.hpp"
#include "../utilities/timer.hpp"
#include "waveblocks/propagators/Parameters.hpp"

/** \file */

namespace waveblocks {
namespace propagators {
        
namespace print = utilities::prettyprint;
namespace utils = utilities;

/**
 * \brief generic abstract propagator class for Hagedorn Wave Packets
 *
 * \tparam N Number of energy levels
 * \tparam D Dimension of space
 * \tparam MultiIndex_t Type of multi index used in the basis shape
 * \tparam MDQR_t Multi-dimensional quadrature rule
 * \tparam Potential_t Type of the Potential to be used
 * \tparam Packet_t Type of the Wavepacket to be propagated
 */
template <typename Propagator_t, int N, int D, typename MultiIndex_t, typename MDQR_t, typename Potential_t, typename Packet_t, typename Coef_t>
class Propagator {
	
	protected:

		Packet_t& wpacket_; ///< wave packet to be propagated
		Potential_t& V_; ///< potential energy
		Coef_t splitCoef_;
		CMatrix<Eigen::Dynamic,Eigen::Dynamic> F_;


	public:

		/**
		 * \brief Propagator constructor
		 *
		 * \param pack wave packet to be propagated
		 * \param V potential to use for propagation
		 */
		Propagator(Packet_t& pack, Potential_t& V, Coef_t coef = SplitCoefs<0,0>({},{}))
		 : wpacket_(pack)
		 , V_(V)
		 , splitCoef_(coef)
		{
			this->onConstruct();
		}

		/**
		 * \brief Helper for constructing propagator for scalar wave packets
		 *
		 * This function is needed because enable_if does not work on constructors
		 * Asserts that number of levels is compatible with the packet.
		 *
		 * \tparam U Dummy template parameter, neccessary for enable_if
		 */
		template <typename U=Packet_t>
		typename std::enable_if<std::is_same<U,wavepackets::ScalarHaWp<D,MultiIndex_t>>::value,void>::type
		onConstruct() {
			static_assert(N==1,"Scalar wave packets must have N==1");
		}


		/**
		 * \brief Helper for constructing propagator for multi-level wave packets
		 *
		 * This function is needed because enable_if does not work on constructors.
		 *
		 * Asserts that number of levels is compatible with the packet and
		 * resizes the matrix F_ to the correct size for future computations.
		 *
		 * \tparam U Dummy template parameter, neccessary for enable_if
		 */
		template <typename U=Packet_t>
		typename std::enable_if<!std::is_same<U,wavepackets::ScalarHaWp<D,MultiIndex_t>>::value,void>::type
		onConstruct() {
			static_assert(N>1,"Multi-Level wave packets must have N>1");
			unsigned size = 0;
			for(auto& comp : wpacket_.components()) {
				size += comp.coefficients().size();
			}
			F_.resize(size,size);
		}


/////////////////////////////////////////////////////////////////////////////////
// Quadratic Potential Energy Operator U
/////////////////////////////////////////////////////////////////////////////////

		/**
		 * \brief function for time evolution
		 *
		 * \param T Size of total time interval
		 * \param Dt Size of each time step
		 * \param callback An optional callback function that is called before doing any time step and at the end of the propagation
		 *  The callback function must take two arguments: the index of the current iteration (unsigned integer) and the current time (real_t)
		 */
		inline void evolve(const real_t T, const real_t Dt, const std::function<void(unsigned, real_t)> callback = [](unsigned,real_t) {}) {

			if(T<0) return;

			real_t t = 0;
			unsigned M = std::round(T/Dt);

			{
				const bool scalar = std::is_same<Packet_t,wavepackets::ScalarHaWp<D,MultiIndex_t>>::value;
				const bool hom = std::is_same<Packet_t,wavepackets::InhomogeneousHaWp<D,MultiIndex_t>>::value;
				std::cout << "\n\n";
				print::title(getName() + " Propagator");
				print::pair("Wave Packet",(scalar ? "Single-Level" : (std::string("Multi-Level") + (hom ? ", Homogeneous" : ", Inhomogeneous"))));
				print::pair("D (number of dimensions)",D);
				print::pair("N (number of energy levels)",N);
				// print::pair("Quadrature Rule","MDQR");
				// print::pair("Order of the scheme","ORDER");
				print::separator();
				print::pair("Time Span T",T);
				print::pair("Stepsize Dt",Dt);
				print::pair("Number of Timesteps",M);
				print::separator();
				std::cout << "\n";
			}

			utilities::Timer timer;

			timer.start(); // ------------------------------- TIMER START

			pre_propagate(Dt);

			callback(0,t);
			for(unsigned m=1; m<=M; ++m) {
				t += Dt;
				propagate(Dt);
				callback(m,t);
				print::pair("Time t",t,"\r");
			}

			timer.stop(); // -------------------------------- TIMER STOP

			print::pair("","COMPLETE","\r");
			print::separator();
			print::pair("Computation time [s]",timer.seconds());
			print::separator();

			post_propagate(Dt);

		}



/////////////////////////////////////////////////////////////////////////////////
// Functions to be provided by derived classes (CRTP)
/////////////////////////////////////////////////////////////////////////////////

	protected:

		/** \brief get the name of the current propagator */
		std::string getName() { return static_cast<Propagator_t*>(this)->getName(); }

		/**
		 * \brief do the main propagation loop
		 * 
		 * \param Dt Size of timestep to be performed
		 */
		inline void propagate(const real_t Dt) { static_cast<Propagator_t*>(this)->propagate(Dt); }

		/**
		 * \brief pre-propagation work
		 *
		 * \param Dt Size of timestep to be performed
		 */
		inline void pre_propagate(const real_t Dt) { static_cast<Propagator_t*>(this)->pre_propagate(Dt); }

		/**
		 * \brief post-propagation work
		 *
		 * \param Dt Size of timestep to be performed
		 */
		inline void post_propagate(const real_t Dt) { static_cast<Propagator_t*>(this)->post_propagate(Dt); }



/////////////////////////////////////////////////////////////////////////////////
// Kinetic Energy Operator T
/////////////////////////////////////////////////////////////////////////////////

	protected:
		/**
		 * \brief single step with kinetic energy operator T, homogeneous wave packet
		 *
		 * Carry out a single timestep of size h with the kinetic energy operator T for a homogeneous wavepacket.
		 * Calls stepT_params on the complete wavepacket
		 *
		 * \param h Size of timestep
		 */
		// Homogeneous
        template <typename P=Packet_t>
        typename std::enable_if<!std::is_same<P,wavepackets::InhomogeneousHaWp<D,MultiIndex_t>>::value,void>::type
		stepT(const real_t h) {
			// Homogeneous
			stepT_params(h,wpacket_.parameters());
		}
		/**
		 * \brief single step with kinetic energy operator T, inhomogeneous wave packet
		 *
		 * Carry out a single timestep of size h with the kinetic energy operator T for an inhomogeneous wavepacket.
		 * Calls stepT_params on each component (energy level) of the wavepacket
		 *
		 * \param h Size of timestep
		 */
        template <typename P=Packet_t>
        typename std::enable_if<std::is_same<P,wavepackets::InhomogeneousHaWp<D,MultiIndex_t>>::value,void>::type
		stepT(const real_t h) {
			// Inhomogeneous
			for(auto& comp : wpacket_.components()) {
				stepT_params(h,comp.parameters());
			}
		}
	
	private:
		/**
		 * \brief update q,Q,S according to kinetic operator T
		 *
		 * \param h Size of timestep
		 * \param params Wave packet parameters to be updated
		 */
		inline void stepT_params(const real_t h, wavepackets::HaWpParamSet<D>& params) {
			// NB: remove mass if not required
			// RMatrix<Eigen::Dynamic,Eigen::Dynamic> Minv = RMatrix<D,D>::Identity();
			params.updateq( +h * /*Minv*/params.p() ); ///< q = q + h * M^{-1} * p
			params.updateQ( +h * /*Minv*/params.P() ); ///< Q = Q + h * M^{-1} * P
			params.updateS( +.5f*h * params.p().dot(/*Minv*/params.p()) ); ///< S = S + h/2 * p^T M p
		}



/////////////////////////////////////////////////////////////////////////////////
// Quadratic Potential Energy Operator U
/////////////////////////////////////////////////////////////////////////////////
		
	protected:
		/**
		 * \brief single step with quadratic potential energy operator U, homogeneous wave packet
		 *
		 * Carry out a single timestep of size h with the quadratic potential energy operator U for a homogeneous wavepacket.
		 * Calls stepU_params on the complete wavepacket
		 *
		 * \param h Size of timestep
		 */
        template <typename P=Packet_t>
        typename std::enable_if<!std::is_same<P,wavepackets::InhomogeneousHaWp<D,MultiIndex_t>>::value,void>::type
		stepU(const real_t h) {
			// Homogeneous
			// taylorV[0,1,2] = [V,DV,DDV] = [evaluation,jacobian,hessian]
			auto& params = wpacket_.parameters();
			const auto& taylorV = V_.get_leading_level().taylor_at(/* q = */ complex_t(1,0) * utils::Squeeze<D,RVector<D>>::apply(params.q()));
			stepU_params(h,params,std::get<0>(taylorV),std::get<1>(taylorV),std::get<2>(taylorV));
		}
		/**
		 * \brief single step with quadratic potential energy operator U, inhomogeneous wave packet
		 *
		 * Carry out a single timestep of size h with the quadratic potential energy operator U for an inhomogeneous wavepacket.
		 * Calls stepU_params on each component (energy level) of the wavepacket
		 *
		 * \param h Size of timestep
		 */
        template <typename P=Packet_t>
        typename std::enable_if<std::is_same<P,wavepackets::InhomogeneousHaWp<D,MultiIndex_t>>::value,void>::type
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
		 *
		 * \tparam V_T Evaluation type for the potential
		 * \tparam DV_T Evaluation type for the jacobian
		 * \tparam DDV_T Evaluation type for the hessian
		 *
		 * \param h Size of timestep
		 * \param params Wave packet parameters to be updated
		 * \param potential Evaluation of the Potential
		 * \param jacobian Evaluation of the Jacobian of the potential
		 * \param hessian Evaluation of the Hessian of the potential
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
		 * Carry out a single timestep of size h with the potential energy remainder W=V-U
		 * 
		 * \param h Size of timestep
		 */
		inline void stepW(const real_t h) {
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

			buildF();
			CVector<Eigen::Dynamic> coefs = utils::PacketToCoefficients<Packet_t>::to(wpacket_); // get coefficients from packet
			complex_t factor(0,-h/(wpacket_.eps()*wpacket_.eps()));
			coefs = (factor*F_).exp() * coefs; ///< c = exp(-i*h/eps^2 * F) * c
			utils::PacketToCoefficients<Packet_t>::from(coefs,wpacket_); // update packet from coefficients
		}

		/** \brief convenience function that calls buildF(F_) */
		inline void buildF() { buildF(F_); }

		/**
		 * \brief build the interaction matrix F for a multi-level system
		 *
		 * \param M Result matrix
		 */
        template <int N_LEVELS=N>
        typename std::enable_if<(N_LEVELS>1),void>::type // N>1
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

			M = innerproducts::VectorInnerProduct<D,MultiIndex_t,MDQR_t>::build_matrix(wpacket_,op);

		}

		/**
		 * \brief build the interaction matrix F for a 1-level system
		 *
		 * \param M Result matrix
		 */
        template <int N_LEVELS=N>
        typename std::enable_if<(N_LEVELS==1),void>::type // N=1
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

			M = innerproducts::HomogeneousInnerProduct<D,MultiIndex_t,MDQR_t>::build_matrix(wpacket_,op);
		}



/////////////////////////////////////////////////////////////////////////////////
// Int-Split
/////////////////////////////////////////////////////////////////////////////////

	protected:
		/**
		 * \brief split the timestep of size Dt into M smaller timesteps
		 *
		 * Splits the timestep Dt into M smaller intervals of size dt=Dt/M and propagates
		 * the wavepacket in an alternating way with the operators T and U (starting with T).
		 * The coefficients of the splitting are given by a splitting scheme.
		 *
		 * \tparam S_T Size of the array of coefficients for operator T
		 * \tparam S_U Size of the array of coefficients for operator U
		 * \param Dt Size of time interval
		 * \param M Number of sub-intervals
		 * \param coefs A pair of arrays of size S_T and S_U respectively.
		 *     Since the operators T and U are applied alternately (starting with T),
		 *     it must hold that S_T==S_U OR S_T == S_U+1
		 */
		template <const std::size_t S_T, const std::size_t S_U>
		void intSplit(const real_t Dt, const unsigned M, const SplitCoefs<S_T,S_U> coefs) {
			static_assert(S_T==S_U || S_T==S_U+1,"invalid combination of coefficient array sizes");
			const real_t dt = Dt/M;
			for(unsigned m=0; m<M; ++m) {
				TU<Propagator_t,0,0,S_T,S_U>::split(*static_cast<Propagator_t*>(this),coefs.a,coefs.b,dt);
			}
		}
		

	private:

		// forward declarations
		/** \brief templated struct for splitting into T and remainder */
		template <typename PROPAGATOR, std::size_t N_T, std::size_t N_U, std::size_t S_T, std::size_t S_U> struct TU;
		/** \brief templated struct for splitting into U and remainder */
		template <typename PROPAGATOR, std::size_t N_U, std::size_t N_T, std::size_t S_U, std::size_t S_T> struct UT;

		// termination conditions
		/** \brief termination condition N_T == S_T and N_U == S_U */
		template <typename PROPAGATOR, std::size_t S_T, std::size_t S_U>
		struct TU<PROPAGATOR,S_T,S_U,S_T,S_U> {
			static void split(PROPAGATOR&, const std::array<real_t,S_T>&, const std::array<real_t,S_U>&, const real_t) {}
		};
		/** \brief termination condition N_T == S_T and N_U == S_U */
		template <typename PROPAGATOR, std::size_t S_U, std::size_t S_T>
		struct UT<PROPAGATOR,S_U,S_T,S_U,S_T> {
			static void split(PROPAGATOR&, const std::array<real_t,S_U>&, const std::array<real_t,S_T>&, const real_t) {}
		};


		/**
		 * A templated struct with one static function that alternately applies the operators T and U, starting with T
		 * Neccessary since partial specialization of function templates is not allowed in cpp
		 * (only partial specialization of types/classes)
		 *
		 * \tparam PROPAGATOR type of the propagator
		 * \tparam N_T current index of coefficient for operator T
		 * \tparam N_U current index of coefficient for operator U
		 * \tparam S_T size of coefficient array coefT
		 * \tparam S_U size of coefficient array coefU
		 */
		template <typename PROPAGATOR, std::size_t N_T, std::size_t N_U, std::size_t S_T, std::size_t S_U>
		struct TU {
			/**
			 * \brief alternately apply T and U, starting with U
			 *
			 * \param prop Reference to the current propagator (this pointer is not available from within nested classes)
			 * \param coefT Array of coefficients for propagation with operator T
			 * \param coefU Array of coefficients for propagation with operator U
			 * \param dt Size of the current timestep sub-interval
			 */
			static void split(PROPAGATOR& prop, const std::array<real_t,S_T>& coefT, const std::array<real_t,S_U>& coefU, const real_t dt) {
				static_assert(N_T==N_U || N_T==N_U-1,"ERROR: invalid alternating index pair for TU");
				// print::pair("TU",std::to_string(N_T)+","+std::to_string(N_U));
				/* T     */ prop.stepT(coefT.at(N_T)*dt);
				/* UT... */ UT<PROPAGATOR,N_U,N_T+1,S_U,S_T>::split(prop,coefU,coefT,dt);
			}
		};
		/**
		 * A templated struct with one static function that alternately applies the operators U and T, starting with U
		 * Neccessary since partial specialization of function templates is not allowed in cpp
		 * (only partial specialization of types/classes)
		 *
		 * \tparam PROPAGAUOR type of the propagator
		 * \tparam N_U current index of coefficient for operator U
		 * \tparam N_T current index of coefficient for operator T
		 * \tparam S_U size of coefficient array coefU
		 * \tparam S_T size of coefficient array coefT
		 */
		template <typename PROPAGATOR, std::size_t N_U, std::size_t N_T, std::size_t S_U, std::size_t S_T>
		struct UT {
			/**
			 * \brief alternately apply U and T, starting with U
			 *
			 * \param prop Reference to the current propagator (this pointer is not available from within nested classes)
			 * \param coefU Array of coefficients for propagation with operator U
			 * \param coefT Array of coefficients for propagation with operator T
			 * \param dt Size of the current timestep sub-interval
			 */
			static void split(PROPAGATOR& prop, const std::array<real_t,S_U>& coefU, const std::array<real_t,S_T>& coefT, const real_t dt) {
				static_assert(N_U==N_T || N_U==N_T-1,"ERROR: invalid alternating index pair for UT");
				// print::pair("UT",std::to_string(N_U)+","+std::to_string(N_T));
				/* U     */ prop.stepU(coefU.at(N_U)*dt);
				/* TU... */ TU<PROPAGATOR,N_T,N_U+1,S_T,S_U>::split(prop,coefT,coefU,dt);
			}
		};

};
	
} // namespace propagators
} // namespace waveblocks
