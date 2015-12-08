#pragma once

#include <vector>
#include <tuple>
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <functional>
#include "basic_types.hpp"

namespace waveblocks
{

  // scalar types
  typedef double real_t;
  typedef std::complex<real_t> complex_t;
  typedef int dim_t;

  // matrix types
  template <class I, int R, int C>
  using GMatrix = Eigen::Matrix<I, R, C>;

  template <int R, int C>
  using CMatrix = GMatrix<complex_t, R, C>;

  template <int R, int C>
  using RMatrix = GMatrix<real_t, R, C>;

  // vector types
  template <class I, int R>
  using GVector = GMatrix<I, R, 1>;

  template <int R>
  using CVector = GVector<complex_t, R>;

  template <int R>
  using RVector = GVector<real_t, R>;

  // function type
  template <class P>
  using function_t = std::function<P>;

  template <int D>
  using rD_to_r = function_t<real_t( RVector<D> )>;

  template <int D>
  using rD_to_rD = function_t<RVector<D>( RVector<D> )>;

  template <int D>
  using rD_to_rDxD = function_t<RMatrix<D, D>( RVector<D> )>;

  template <int D>
  using rD_to_c = function_t<complex_t( RVector<D> )>;

  template <int D>
  using rD_to_cD = function_t<CVector<D>( RVector<D> )>;

  template <int D>
  using rD_to_cDxD = function_t<CMatrix<D, D>( RVector<D> )>;

  template <int D, int N>
  using rD_to_cNxN = function_t<CMatrix<N, N>( RVector<D> )>;


  template <int D>
  using cD_to_r = function_t<real_t( CVector<D> )>;

  template <int D>
  using cD_to_rD = function_t<RVector<D>( CVector<D> )>;

  template <int D>
  using cD_to_rDxD = function_t<RMatrix<D, D>( CVector<D> )>;

  using c_to_c = function_t<complex_t(complex_t)>;

  template <int D>
  using cD_to_c = function_t<complex_t( CVector<D> )>;

  template <int D>
  using cD_to_cD = function_t<CVector<D>( CVector<D> )>;

  template <int D>
  using cD_to_cDxD = function_t<CMatrix<D, D>( CVector<D> )>;

  template <int D, int N>
  using cD_to_cNxN = function_t<CMatrix<N, N>( CVector<D> )>;

  // function valued functions
  template <int D, class F>
  using rD_to_function = function_t<F( RVector<D> )>;

  template <int D, int N, class F>
  using rD_to_function_vector = function_t<GVector<F, N>( RVector<D> )>;

  template <int D, int N, class F>
  using rD_to_function_matrix = function_t<GMatrix<F, N, N>( RVector<D> )>;

  template <int D, class F>
  using cD_to_function = function_t<F( CVector<D> )>;

  template <int D, int N, class F>
  using cD_to_function_vector = function_t<GVector<F, N>( CVector<D> )>;

  template <int D, int N, class F>
  using cD_to_function_matrix = function_t<GMatrix<F, N, N>( CVector<D> )>;

  // eigenvalue
  template <int N>
  using eigenvalues_t =
    typename Eigen::EigenSolver<RMatrix<N, N> >::EigenvalueType;

  template <int N>
  using eigenvector_t =
    typename Eigen::EigenSolver<RMatrix<N, N> >::EigenvectorsType;
}
