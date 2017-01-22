#pragma once

#include "waveblocks/propagators/Propagator.hpp"
#include "waveblocks/propagators/SplittingParameters.hpp"

namespace waveblocks {
namespace propagators {

namespace utils = utilities;

/** \file
 * \brief implements the Semiclassical Propagator
 *
 * This class implements the (pre-/post-)propagate functions for the Hagedorn Propagator
 */
template <int N, int D, typename MultiIndex_t, typename MDQR_t, typename Potential_t, typename Packet_t>
class SemiclassicalPropagator : public Propagator<SemiclassicalPropagator<N,D,MultiIndex_t,MDQR_t,Potential_t,Packet_t>,N,D,MultiIndex_t,MDQR_t,Potential_t,Packet_t> {

	public:

		// inherit constructor
		using Propagator<SemiclassicalPropagator<N,D,MultiIndex_t,MDQR_t,Potential_t,Packet_t>,N,D,MultiIndex_t,MDQR_t,Potential_t,Packet_t>::Propagator;

		std::string getName() { return "Semiclassical"; }

		template <typename T=std::pair<std::array<real_t,1>,std::array<real_t,1>>>
		void propagate(const real_t Dt) {

			auto& packet = this->wpacket_;
			int M = std::ceil(1 + Dt/(std::pow(packet.eps(),.75)));
			if(M%2) M++; // make sure M is even such that M/2 + M/2 = M

			this->intSplit(Dt/2,M/2,splitting_parameters::coefLT);
			this->stepW(Dt);
			this->intSplit(Dt/2,M/2,splitting_parameters::coefLT);

		}

		void pre_propagate(const real_t) { /* nothing to do */ }
		void post_propagate(const real_t) { /* nothing to do */ }

};

} // namespace propagators
} // namespace waveblocks
