#pragma once

#define IMPORT_TYPES_FROM(B)                                            \
    using argument_type = typename B::argument_type;                    \
    using potential_type = typename B::potential_type;                  \
    using jacobian_type = typename B::jacobian_type;                    \
    using hessian_type = typename B::hessian_type;                      \
    using potential_evaluation_type = typename B::potential_evaluation_type; \
    using jacobian_evaluation_type = typename B::jacobian_evaluation_type; \
    using hessian_evaluation_type = typename B::hessian_evaluation_type; \
    using potential_return_type = typename B::potential_return_type;    \
    using jacobian_return_type = typename B::jacobian_return_type;      \
    using hessian_return_type = typename B::hessian_return_type;        \
    using local_quadratic_evaluation_type = typename B::local_quadratic_evaluation_type; \
    using local_quadratic_return_type = typename B::local_quadratic_return_type;


namespace waveblocks
{
    namespace potentials
    {
        namespace bases
        {
            template <int N, int D, int C>
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
            using Canonical = Basis<N, D, N>;

            /**
             * \brief Collection of types associated with a matrix potential in eigen basis
             *
             * \tparam N
             * Number of levels (dimension of diagonal vector when evaluated)
             * \tparam D
             * Dimension of argument space
             */
            template <int N, int D>
            using Eigen = Basis<N, D, 1>;

            /**
             * \brief Helper class to ease template specialzations
             *
             * \tparam N
             * Number of levels
             * \tparam D
             * Dimension of argument space
             * \tparam C
             * Number of columns (N for canonical basis, 1 for eigen basis)
             */
            template <int N, int D, int C>
            struct Basis {
                static const int argument_dimension = D;
                static const int number_of_levels = N;
                static const int number_of_columns = C;

                using argument_type = CVector<D>;
                using potential_type = GMatrix<cD_to_c<D>, N, C>;
                using jacobian_type = GMatrix<cD_to_cD<D>, N, C>;
                using hessian_type = GMatrix<cD_to_cDxD<D>, N, C>;

                using potential_evaluation_type = CMatrix<N, C>;
                using jacobian_evaluation_type = GMatrix<CVector<D>, N, C>;
                using hessian_evaluation_type = GMatrix<CMatrix<D, D>, N, C>;
                using local_quadratic_evaluation_type = CMatrix<N, C>;

                using potential_return_type = complex_t;
                using jacobian_return_type = CVector<D>;
                using hessian_return_type = CMatrix<D, D>;
                using local_quadratic_return_type = complex_t;
            };

            template <int N, int C>
            struct Basis<N, 1, C> {
                static const int argument_dimension = 1;
                static const int number_of_levels = N;
                static const int number_of_columns = C;

                using argument_type = complex_t;
                using potential_type = GMatrix<c_to_c, N, C>;
                using jacobian_type = GMatrix<c_to_c, N, C>;
                using hessian_type = GMatrix<c_to_c, N, C>;

                using potential_evaluation_type = CMatrix<N, C>;
                using jacobian_evaluation_type = CMatrix<N, C>;
                using hessian_evaluation_type = CMatrix<N, C>;
                using local_quadratic_evaluation_type = CMatrix<N, C>;

                using potential_return_type = complex_t;
                using jacobian_return_type = complex_t;
                using hessian_return_type = complex_t;
                using local_quadratic_return_type = complex_t;
            };

            template <int D, int C>
            struct Basis<1, D, C> {
                const int argument_dimension = D;
                static const int number_of_levels = 1;
                static const int number_of_columns = C;

                using argument_type = CVector<D>;
                using potential_type = cD_to_c<D>;
                using jacobian_type = cD_to_cD<D>;
                using hessian_type = cD_to_cDxD<D>;

                using potential_evaluation_type = complex_t;
                using jacobian_evaluation_type = CVector<D>;
                using hessian_evaluation_type = CMatrix<D, D>;
                using local_quadratic_evaluation_type = complex_t;


                using potential_return_type = complex_t;
                using jacobian_return_type = CVector<D>;
                using hessian_return_type = CMatrix<D, D>;
                using local_quadratic_return_type = complex_t;
            };

            template <int C>
            struct Basis<1, 1, C> {
                static const int argument_dimension = 1;
                static const int number_of_levels = 1;
                static const int number_of_columns = C;

                using argument_type = complex_t;
                using potential_type = c_to_c;
                using jacobian_type = c_to_c;
                using hessian_type = c_to_c;

                using potential_evaluation_type = complex_t;
                using jacobian_evaluation_type = complex_t;
                using hessian_evaluation_type = complex_t;
                using local_quadratic_evaluation_type = complex_t;


                using potential_return_type = complex_t;
                using jacobian_return_type = complex_t;
                using hessian_return_type = complex_t;
                using local_quadratic_return_type = complex_t;

            };
        }
    }

    template<int N, int D>
    using CanonicalBasis = potentials::bases::Canonical<N,D>;

    template<int N, int D>
    using EigenBasis = potentials::bases::Eigen<N,D>;
}
