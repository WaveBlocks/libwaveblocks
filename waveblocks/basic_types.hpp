#ifndef WAVEBLOCKS_BASIC_TYPES_HPP
#define WAVEBLOCKS_BASIC_TYPES_HPP

#include <complex>
#include <Eigen/Core>

namespace waveblocks {

typedef double real_t;
typedef std::complex<real_t> complex_t;
typedef int dim_t;

template<int R, int C>
using CMatrix = Eigen::Matrix<complex_t,R,C>;

template<int R, int C>
using RMatrix = Eigen::Matrix<real_t,R,C>;

template<int R, int C>
using CArray = Eigen::Array<complex_t,R,C>;

template<int R, int C>
using RArray = Eigen::Array<real_t,R,C>;

template<int N>
using HaWpBasisVector = CArray<Eigen::Dynamic, N>;

using Coefficients = Eigen::Matrix<complex_t, Eigen::Dynamic, 1>;
}

#endif
