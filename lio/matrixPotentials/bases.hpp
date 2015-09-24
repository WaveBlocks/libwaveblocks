#pragma once
#define IMPORT_TYPES_FROM(B, N, D)                                             \
  using argument_type = typename B<N, D>::argument_type \ 
  using potential_type = typename B<N, D>::potential_type;                     \
  using jacobian_type = typename B<N, D>::jacobian_type;                       \
  using hessian_type = typename B<N, D>::hessian_type;                         \
  using potential_evaluation_type =                                            \
      typename B<N, D>::potential_evaluation_type;                             \
  using jacobian_evaluation_type = typename B<N, D>::jacobian_evaluation_type; \
  using hessian_evaluation_type = typename B<N, D>::hessian_evaluation_type;   \
  using potential_return_type = typename B<N, D>::potential_return_type;       \
  using jacobian_return_type = typename B<N, D>::jacobian_return_type;         \
  using hessian_return_type = typename B<N, D>::hessian_return_type;           \
  using local_quadratic_evaluation_type = typename B<N, D>::local_quadratic_evaluation_type;           \
  using local_quadratic_return_type = typename B<N, D>::local_quadratic_return_type;  



namespace waveblocks
{
  namespace matrixPotentials
  {
    namespace bases
    {
      template <int N, int D>
      struct Basis;
      
      /**
       * \brief Collection of types associated with a matrix potential in canonical basis
       * 
       * \tparam N
       * Number of levels (dimension of square matrix when evaluated)
       * \tparam D
       * Dimension of argument space
       */
      template <int N, int D>
      using Canonical = typename Basis<N, D>::Canonical;
      
      /**
       * \brief Collection of types associated with a matrix potential in eigen basis
       * 
       * \tparam N
       * Number of levels (dimension of diagonal vector when evaluated)
       * \tparam D
       * Dimension of argument space
       */
      template <int N, int D>
      using Eigen = typename Basis<N, D>::Eigen;
      
      /**
       * \brief Helper class to ease template specialzations
       * 
       * \tparam N
       * Number of levels
       * \tparam D
       * Dimension of argument space
       */
      template <int N, int D>
      struct Basis {
        struct Canonical {
          using argument_type = CVector<D>;
          using potential_type = GMatrix<cD_to_c<D>, N, N>;
          using jacobian_type = GMatrix<cD_to_cD<D>, N, N>;
          using hessian_type = GMatrix<cD_to_cDxD<D>, N, N>;
          
          using potential_evaluation_type = CMatrix<N, N>;
          using jacobian_evaluation_type = GMatrix<CVector<D>, N, N>;
          using hessian_evaluation_type = GMatrix<CMatrix<D, D>, N, N>;
          using local_quadratic_evaluation_type = CMatrix<N,N>;
          
          using potential_return_type = real_t;
          using jacobian_return_type = RVector<D>;
          using hessian_return_type = RMatrix<D, D>;
          using local_quadratic_return_type = complex_t;
        };
        
        struct Eigen {
          using argument_type = CVector<D>;
          using potential_type = GVector<cD_to_c<D>, N>;
          using jacobian_type = GVector<cD_to_cD<D>, N>;
          using hessian_type = GVector<cD_to_cDxD<D>, N>;
          
          using potential_evaluation_type = CVector<N>;
          using jacobian_evaluation_type = GVector<CVector<D>, N>;
          using hessian_evaluation_type = GVector<CMatrix<D, D>, N>;
          using local_quadratic_evaluation_type = CVector<N>;

          
          using potential_return_type = complex_t;
          using jacobian_return_type = CVector<D>;
          using hessian_return_type = CMatrix<D, D>;
          using local_quadratic_return_type = complex_t;

        };
      };
      template <int D>
      struct Basis<1, D> {
        struct Canonical {
          using argument_type = CVector<D>;
          using potential_type = cD_to_c<D>;
          using jacobian_type = cD_to_cD<D>;
          using hessian_type = cD_to_cDxD<D>;
          
          using potential_evaluation_type = complex_t;
          using jacobian_evaluation_type = CVector<D>;
          using hessian_evaluation_type = CMatrix<D, D>;
          using local_quadratic_evaluation_type = complex_t;

          
          using potential_return_type = complex_t;
          using jacobian_return_type = RVector<D>;
          using hessian_return_type = RMatrix<D, D>;
          using local_quadratic_return_type = complex_t;

        };
        
        using Eigen = Canonical;

      };
    }
  }
}
