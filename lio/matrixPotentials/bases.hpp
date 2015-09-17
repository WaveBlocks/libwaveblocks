#pragma once
#define IMPORT_TYPES_FROM(B,N,D)  \
  using scalar_type = typename B<N,D>::scalar_type;                        \
  using potential_type = typename B<N,D>::potential_type;                       \
  using jacobian_type = typename B<N,D>::jacobian_type;                         \
  using hessian_type = typename B<N,D>::hessian_type;                           \
  using local_quadratic_type = typename B<N,D>::local_quadratic_type;           \
  using potential_evaluation_type = typename B<N,D>::potential_evaluation_type; \
  using jacobian_evaluation_type = typename B<N,D>::jacobian_evaluation_type;   \
  using hessian_evaluation_type = typename B<N,D>::hessian_evaluation_type;     \
  using potential_return_type = typename B<N,D>::potential_return_type;         \
  using jacobian_return_type = typename B<N,D>::jacobian_return_type;           \
  using hessian_return_type = typename B<N,D>::hessian_return_type              \
  
namespace lio { 
  namespace matrixPotentials { 
    namespace bases {
      template <int N, int D>
      struct Basis;

      template <int N, int D>
      using Canonical = typename Basis<N,D>::Canonical;

      template< int N, int D>
      using Eigen = typename Basis<N,D>::Eigen;


      template<int N, int D>
      struct Basis { 
        struct Canonical {
          using scalar_type = real_t;
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

        struct Eigen {
          using scalar_type = complex_t;
          using potential_type = GVector<rD_to_c<D>, N>;
          using jacobian_type = GVector<rD_to_cD<D>, N>;
          using hessian_type = GVector<rD_to_cDxD<D>, N>;
          using transformation_type = rD_to_cNxN<D, N>;
          using local_quadratic_type = rD_to_function_vector<D, N, rD_to_c<D> >;
          
          using potential_evaluation_type = CVector<N>;
          using jacobian_evaluation_type = GVector<CVector<D>, N>;
          using hessian_evaluation_type = GVector<CMatrix<D, D>, N>;
          
          using potential_return_type = complex_t;
          using jacobian_return_type = CVector<D>;
          using hessian_return_type = CMatrix<D, D>;
        };
      };
      template<int D>
      struct Basis<1,D> {
        struct Canonical {
          using scalar_type = real_t;
          using potential_type = rD_to_r<D>;
          using jacobian_type = rD_to_rD<D>;
          using hessian_type = rD_to_rDxD<D>;
          using local_quadratic_type = rD_to_function<D, rD_to_r<D> >;
          
          using potential_evaluation_type = real_t;
          using jacobian_evaluation_type = RVector<D>;
          using hessian_evaluation_type = RMatrix<D, D>;
          
          using potential_return_type = real_t;
          using jacobian_return_type = RVector<D>;
          using hessian_return_type = RMatrix<D, D>;
        };
        
        //~ using Eigen = Canonical;
        struct Eigen {
          using scalar_type = complex_t;
          using potential_type = rD_to_c<D>;
          using jacobian_type = rD_to_cD<D>;
          using hessian_type = rD_to_cDxD<D>;
          using local_quadratic_type = rD_to_function<D, rD_to_c<D> >;
          
          using potential_evaluation_type = complex_t;
          using jacobian_evaluation_type = CVector<D>;
          using hessian_evaluation_type = CMatrix<D, D>;
          
          using potential_return_type = complex_t;
          using jacobian_return_type = CVector<D>;
          using hessian_return_type = CMatrix<D, D>;
         
        };
      };
    }
  }
}
