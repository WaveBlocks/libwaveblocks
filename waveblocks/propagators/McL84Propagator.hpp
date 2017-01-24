#pragma once

#include <cmath>

#include "waveblocks/propagators/Propagator.hpp"
#include "waveblocks/propagators/SplittingParameters.hpp"

namespace waveblocks {
namespace propagators {

namespace utils = utilities;

/** \file
 * \brief implements the McL84 Propagator
 *
 * This class implements the (pre-/post-)propagate functions for the McL84 Propagator
 */
template <int N, int D, typename MultiIndex_t, typename MDQR_t, typename Potential_t, typename Packet_t, typename Coef_t = SplitCoefs<0,0>>
class McL84Propagator : public Propagator<McL84Propagator<N,D,MultiIndex_t,MDQR_t,Potential_t,Packet_t,Coef_t>,N,D,MultiIndex_t,MDQR_t,Potential_t,Packet_t,Coef_t> {

	public:

		// inherit constructor
		using Propagator<McL84Propagator<N,D,MultiIndex_t,MDQR_t,Potential_t,Packet_t,Coef_t>,N,D,MultiIndex_t,MDQR_t,Potential_t,Packet_t,Coef_t>::Propagator;

		std::string getName() { return "McL84"; }

		void propagate(const real_t Dt) {

			auto& packet = this->wpacket_;
			int r = 2; // TODO: r = innerorder
            real_t alpha = 2.;
            real_t beta = 2.;
            const int M = std::pow(Dt,r-beta) / std::pow(std::pow(packet.eps(),alpha+2.) , (1./r));

			this->intSplit(a_.at(0)*Dt,M,this->splitCoef_);
			this->stepW(b_.at(0)*Dt);
			this->intSplit(a_.at(1)*Dt,M,this->splitCoef_);
			this->stepW(b_.at(1)*Dt);
			this->intSplit(a_.at(2)*Dt,M,this->splitCoef_);
			this->stepW(b_.at(2)*Dt);

		}

		void pre_propagate(const real_t) { /* nothing to do */ }
		void post_propagate(const real_t) { /* nothing to do */ }

		static constexpr std::array<real_t,3> a_ = {0.07534696026989288842,0.51791685468825678230,-0.09326381495814967072};
		static constexpr std::array<real_t,3> b_ = {0.19022593937367661925,0.84652407044352625706,-1.07350001963440575260};

};

} // namespace propagators
} // namespace waveblocks
