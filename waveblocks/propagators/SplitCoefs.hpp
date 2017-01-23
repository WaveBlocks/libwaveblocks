#pragma once

#include <array>
#include "waveblocks/types.hpp"

namespace waveblocks {
namespace propagators {

template <std::size_t S_A=0, std::size_t S_B=S_A>
struct SplitCoefs {
	SplitCoefs(std::array<real_t,S_A> coefA, std::array<real_t,S_B> coefB)
	 : coefs_(coefA,coefB)
	 , a(coefs_.first)
	 , b(coefs_.second)
	{}
	const std::pair<std::array<real_t,S_A>,std::array<real_t,S_B>> coefs_;
	const std::array<real_t,S_A>& a;
	const std::array<real_t,S_B>& b;
};

} // namespace waveblocks
} // namespace propagators
