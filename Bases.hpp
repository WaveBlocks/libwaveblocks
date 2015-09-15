#pragma once

template <int N, int D>
struct CanonicalBasis {
  using potential_type = GMatrix<rD_to_r<D>, N, N>;
  using jacobian_type = GMatrix<rD_to_rD<D>, N, N>;
  using hessian_type = GMatrix<rD_to_rDxD<D>, N, N>;
  using transformation_type = rD_to_cNxN<D, N>;
  using local_quadratic_type = rD_to_function_matrix<D, N, rD_to_r<D> >;
  
  using potential_evaluation_type = RMatrix<N, N>;
  using jacobian_evaluation_type = GMatrix<RVector<D>, N, N>;
  using hessian_evaluation_type = GMatrix<RMatrix<D, D>, N, N>;
  
  using potential_return_type = real_t;
  using jacobian_return_type = RVector<D>;
  using hessian_return_type = RMatrix<D, D>;
};

template <int N, int D>
struct EigenBasis {
  using potential_type = GVector<rD_to_c<D>, N>;
  using jacobian_type = GVector<rD_to_cD<D>, N>;
  using hessian_type = GVector<rD_to_cDxD<D>, N>;
  using transformation_type = rD_to_cNxN<D, N>;
  using local_quadratic_type = rD_to_function_vector<D, N, rD_to_c<D> >;
  using local_remainder_type = rD_to_function_vector<D, N, rD_to_c<D> >;
  
  using potential_evaluation_type = CVector<N>;
  using jacobian_evaluation_type = GVector<CVector<D>, N>;
  using hessian_evaluation_type = GVector<CMatrix<D, D>, N>;
  
  using potential_return_type = complex_t;
  using jacobian_return_type = CVector<D>;
  using hessian_return_type = CMatrix<D, D>;
};
