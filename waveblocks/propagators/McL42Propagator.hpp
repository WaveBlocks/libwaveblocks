#pragma once

#include <cmath>

#include "waveblocks/propagators/Propagator.hpp"
#include "waveblocks/propagators/SplittingParameters.hpp"

namespace waveblocks {
namespace propagators {

namespace utils = utilities;

/** \file
 * \brief implements the McL42 Propagator
 *
 * This class implements the (pre-/post-)propagate functions for the McL42 Propagator
 */
template <int N, int D, typename MultiIndex_t, typename MDQR_t, typename Potential_t, typename Packet_t, typename Coef_t = SplitCoefs<0,0>>
class McL42Propagator : public Propagator<McL42Propagator<N,D,MultiIndex_t,MDQR_t,Potential_t,Packet_t,Coef_t>,N,D,MultiIndex_t,MDQR_t,Potential_t,Packet_t,Coef_t> {

	public:

		// inherit constructor
		using Propagator<McL42Propagator<N,D,MultiIndex_t,MDQR_t,Potential_t,Packet_t,Coef_t>,N,D,MultiIndex_t,MDQR_t,Potential_t,Packet_t,Coef_t>::Propagator;

		std::string getName() { return "McL42"; }

		void propagate(const real_t Dt) {

			auto& packet = this->wpacket_;
			int r = 2; // TODO: r = innerorder
            real_t alpha = 2.;
            real_t beta = 2.;
            int M = std::pow(Dt,r-beta) / std::pow(std::pow(packet.eps(),alpha+2.) , (1./r));
			print::separator();
			print::pair("M",M);
			print::separator();
			this->intSplit(a_.at(0)*Dt,M,this->splitCoef_);
			this->stepW(b_.at(0)*Dt);
			this->intSplit(a_.at(1)*Dt,M,this->splitCoef_);
			this->stepW(b_.at(1)*Dt);
			this->intSplit(a_.at(2)*Dt,M,this->splitCoef_);

		}

		void pre_propagate(const real_t) { /* nothing to do */ }
		void post_propagate(const real_t) { /* nothing to do */ }

		static constexpr std::array<real_t,3> a_ = { (3.-std::sqrt(3))/6., 1./std::sqrt(3), (3.-std::sqrt(3))/6. };
		static constexpr std::array<real_t,2> b_ = { .5,.5 };

};

} // namespace propagators
} // namespace waveblocks
