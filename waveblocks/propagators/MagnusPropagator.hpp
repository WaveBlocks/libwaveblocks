#pragma once

#include "waveblocks/propagators/Propagator.hpp"

/** \file */

namespace waveblocks {
namespace propagators {

namespace utils = utilities;

/**
 * \brief implements the Magnus Propagator
 *
 * This class implements the (pre-/post-)propagate functions for the Magnus Propagator
 */
template <int N, int D, typename MultiIndex_t, typename MDQR_t, typename Potential_t, typename Packet_t, typename Coef_t = SplitCoefs<0,0>>
class MagnusPropagator : public Propagator<MagnusPropagator<N,D,MultiIndex_t,MDQR_t,Potential_t,Packet_t,Coef_t>,N,D,MultiIndex_t,MDQR_t,Potential_t,Packet_t,Coef_t> {

	public:

		// inherit constructor
		using Propagator<MagnusPropagator<N,D,MultiIndex_t,MDQR_t,Potential_t,Packet_t,Coef_t>,N,D,MultiIndex_t,MDQR_t,Potential_t,Packet_t,Coef_t>::Propagator;

		std::string getName() { return "Magnus"; }

		void propagate(const real_t Dt) {

			auto& packet = this->wpacket_;
			CMatrix<Eigen::Dynamic,Eigen::Dynamic> A1, A2;
			// Gauss Legendre coefficients on [0,Dt]:
			real_t h1 = (3.-std::sqrt(3.))/6. * Dt;
			real_t h2 = std::sqrt(3.)/6. * Dt;
			int M1 = std::ceil(1 + std::sqrt(h1*std::pow(packet.eps(),-3./8.)));
			int M2 = std::ceil(1 + std::sqrt(h2*std::pow(packet.eps(),-3./8.)));

			this->intSplit(h1,M1,this->splitCoef_);
			this->buildF(A1);
			A1 *= complex_t(0,-Dt/(packet.eps()*packet.eps()));
			this->intSplit(h2,M2,this->splitCoef_);
			this->buildF(A2);
			A2 *= complex_t(0,-Dt/(packet.eps()*packet.eps()));
			this->F_ = .5*(A1+A2) + std::sqrt(3)/12*(A1*A2-A2*A1);
			CVector<Eigen::Dynamic> coefs = utils::PacketToCoefficients<Packet_t>::to(packet); // get coefficients from packet
			coefs = (this->F_).exp() * coefs;
			utils::PacketToCoefficients<Packet_t>::from(coefs,packet); // update packet from coefficients
			this->intSplit(h1,M1,this->splitCoef_);

		}

		void pre_propagate(const real_t) { /* nothing to do */ }
		void post_propagate(const real_t) { /* nothing to do */ }

};

} // namespace propagators
} // namespace waveblocks
