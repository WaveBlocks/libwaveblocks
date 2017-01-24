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

			auto& packet = this->wpacket_;
			const int M = 1 + std::floor(std::sqrt(Dt*std::pow(packet.eps(),-.75)));

			for(unsigned j=0; j<k_; ++j) {
				this->stepW(alpha_[j]*Dt);
				this->intSplit(beta_[j]*Dt,M,this->splitCoef_);
			}

		}

		void pre_propagate(const real_t Dt) {

			auto& packet = this->wpacket_;
			const int M = 1 + std::floor(std::sqrt(Dt*std::pow(packet.eps(),-.75)));

			for(unsigned j=0; j<v_; ++j) {
				this->intSplit(-Z_[j]*Dt,M,this->splitCoef_);
				this->stepW(-Y_[j]*Dt);
			}

		}

		void post_propagate(const real_t Dt) {
			
			auto& packet = this->wpacket_;
			const int M = 1 + std::floor(std::sqrt(Dt*std::pow(packet.eps(),-.75)));

			for(unsigned j=0; j<v_; ++j) {
				this->stepW(Y_[j]*Dt);
				this->intSplit(Z_[j]*Dt,M,this->splitCoef_);
			}

		}
			
		static const int k_ = 4;
		static const int v_ = 6;
		static constexpr std::array<real_t,k_> alpha_ = { 0.0, 1.5171479707207228, -2.0342959414414454, 1.5171479707207228 };
		static constexpr std::array<real_t,k_> beta_ = { 0.5600879810924619, -0.06008798109246194, -0.06008798109246194, 0.5600879810924619 };
		static constexpr std::array<real_t,v_> Y_ = { -1.621810118086801, 0.0061709468110142, 0.8348493592472594, -0.0511253369989315, 0.5633782670698199, -0.5 };
		static constexpr std::array<real_t,v_> Z_ = { -0.3346222298730, 1.097567990732164, -1.038088746096783, 0.6234776317921379, -1.102753206303191, -0.0141183222088869 };

};

} // namespace propagators
} // namespace waveblocks
