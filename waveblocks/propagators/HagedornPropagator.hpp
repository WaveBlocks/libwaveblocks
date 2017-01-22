#pragma once

#include "waveblocks/propagators/Propagator.hpp"
#include "waveblocks/propagators/SplittingParameters.hpp"

namespace waveblocks {
namespace propagators {

namespace utils = utilities;

/** \file
 * \brief implements the Hagedorn Propagator
 *
 * This class implements the (pre-/post-)propagate functions for the Hagedorn Propagator
 */
template <int N, int D, typename MultiIndex_t, typename MDQR_t, typename Potential_t, typename Packet_t>
class HagedornPropagator : public Propagator<HagedornPropagator<N,D,MultiIndex_t,MDQR_t,Potential_t,Packet_t>,N,D,MultiIndex_t,MDQR_t,Potential_t,Packet_t>, public SplittingParameters {

	public:

		// inherit constructor
		using Propagator<HagedornPropagator<N,D,MultiIndex_t,MDQR_t,Potential_t,Packet_t>,N,D,MultiIndex_t,MDQR_t,Potential_t,Packet_t>::Propagator;

		std::string getName() { return "Hagedorn"; }

		void propagate(const real_t Dt) {
			this->stepT(Dt/2);
			this->stepU(Dt);
			this->stepW(Dt);
			this->stepT(Dt/2);
		}

		void pre_propagate(const real_t) { /* nothing to do */ }
		void post_propagate(const real_t) { /* nothing to do */ }

};


} // namespace propagators
} // namespace waveblocks
