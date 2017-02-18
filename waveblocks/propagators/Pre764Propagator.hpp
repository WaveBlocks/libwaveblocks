#pragma once

#include "waveblocks/propagators/Propagator.hpp"

/** \file */

namespace waveblocks {
namespace propagators {

namespace utils = utilities;

/**
 * \brief implements the Pre764 Propagator
 *
 * This class implements the (pre-/post-)propagate functions for the Pre764 Propagator
 */
template <int N, int D, typename MultiIndex_t, typename MDQR_t, typename Potential_t, typename Packet_t, typename Coef_t = SplitCoefs<0,0>>
class Pre764Propagator : public Propagator<Pre764Propagator<N,D,MultiIndex_t,MDQR_t,Potential_t,Packet_t,Coef_t>,N,D,MultiIndex_t,MDQR_t,Potential_t,Packet_t,Coef_t> {

	public:
	
		// inherit constructor
		using Propagator<Pre764Propagator<N,D,MultiIndex_t,MDQR_t,Potential_t,Packet_t,Coef_t>,N,D,MultiIndex_t,MDQR_t,Potential_t,Packet_t,Coef_t>::Propagator;

		std::string getName() { return "Pre764"; }

		void propagate(const real_t Dt) {
			
			constexpr auto& alpha = processing_splitting_parameters::coefBCR764.a;
			constexpr auto& beta = processing_splitting_parameters::coefBCR764.b;

			static_assert(alpha.size() == k_, "Incorrect size of coefficient array alpha");
			static_assert(beta.size() == k_,"Incorrect size of coefficient array beta");

			// auto& packet = this->wpacket_;
			const int M = 4; // 1 + std::floor(std::sqrt(Dt*std::pow(packet.eps(),-.75)));

			for(int j=0; j<k_; ++j) {
				this->stepW(alpha[j]*Dt);
				this->intSplit(beta[j]*Dt,M,this->splitCoef_);
			}

		}

		void pre_propagate(const real_t Dt) {

			constexpr auto& y = processing_splitting_parameters::prepostBCR764.a;
			constexpr auto& z = processing_splitting_parameters::prepostBCR764.b;

			static_assert(y.size() == v_,"Incorrect size of coefficient array Y");
			static_assert(z.size() == v_,"Incorrect size of coefficient array Z");

			// auto& packet = this->wpacket_;
			const int M = 4; // 1 + std::floor(std::sqrt(Dt*std::pow(packet.eps(),-.75)));

			for(int j=0; j<v_; ++j) {
				this->intSplit(-z[j]*Dt,M,this->splitCoef_);
				this->stepW(-y[j]*Dt);
			}

		}

		void post_propagate(const real_t Dt) {
			
			constexpr auto& y = processing_splitting_parameters::prepostBCR764.a;
			constexpr auto& z = processing_splitting_parameters::prepostBCR764.b;

			static_assert(y.size() == v_,"Incorrect size of coefficient array Y");
			static_assert(z.size() == v_,"Incorrect size of coefficient array Z");

			// auto& packet = this->wpacket_;
			const int M = 4; // 1 + std::floor(std::sqrt(Dt*std::pow(packet.eps(),-.75)));

			for(int j=(v_-1); j>=0; --j) {
				this->stepW(y[j]*Dt);
				this->intSplit(z[j]*Dt,M,this->splitCoef_);
			}

		}
			
		static const int k_ = 4;
		static const int v_ = 6;

};

} // namespace propagators
} // namespace waveblocks
