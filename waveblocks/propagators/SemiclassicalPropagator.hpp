#pragma once

#include "waveblocks/propagators/Propagator.hpp"

/** \file */

namespace waveblocks {
namespace propagators {

namespace utils = utilities;

/**
 * \brief implements the Semiclassical Propagator
 *
 * This class implements the (pre-/post-)propagate functions for the Semiclassical Propagator
 */
template <int N, int D, typename MultiIndex_t, typename MDQR_t, typename Potential_t, typename Packet_t, typename Coef_t = SplitCoefs<0,0>>
class SemiclassicalPropagator : public Propagator<SemiclassicalPropagator<N,D,MultiIndex_t,MDQR_t,Potential_t,Packet_t,Coef_t>,N,D,MultiIndex_t,MDQR_t,Potential_t,Packet_t,Coef_t> {

	public:

		// inherit constructor
		using Propagator<SemiclassicalPropagator<N,D,MultiIndex_t,MDQR_t,Potential_t,Packet_t,Coef_t>,N,D,MultiIndex_t,MDQR_t,Potential_t,Packet_t,Coef_t>::Propagator;

		std::string getName() { return "Semiclassical"; }

		void propagate(const real_t Dt) {

			auto& packet = this->wpacket_;
			int M = std::ceil(1 + std::sqrt{Dt}/(std::pow(packet.eps(),.75)));
			if(M%2) M++; // make sure M is even such that M/2 + M/2 = M

			this->intSplit(Dt/2,M/2,this->splitCoef_);
			this->stepW(Dt);
			this->intSplit(Dt/2,M/2,this->splitCoef_);

		}

		void pre_propagate(const real_t) { /* nothing to do */ }
		void post_propagate(const real_t) { /* nothing to do */ }

};

} // namespace propagators
} // namespace waveblocks
