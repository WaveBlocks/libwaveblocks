#pragma once

#include <array>

#include "basic_types.hpp"

template<size_t order>
struct QuadratureRule
{
  std::array<real_t, order> nodes;
  std::array<real_t, order> weights;
};
