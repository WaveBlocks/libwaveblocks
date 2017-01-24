#pragma once

#include <array>
#include "waveblocks/types.hpp"

namespace waveblocks {
namespace propagators {

template <std::size_t S_A=0, std::size_t S_B=S_A>
struct SplitCoefs {
	SplitCoefs(std::array<real_t,S_A> coefA, std::array<real_t,S_B> coefB)
	 : a(coefA)
	 , b(coefB)
	{}
	const std::array<real_t,S_A> a;
	const std::array<real_t,S_B> b;
};

} // namespace waveblocks
} // namespace propagators
