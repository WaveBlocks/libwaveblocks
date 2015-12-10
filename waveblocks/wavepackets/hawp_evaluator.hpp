#pragma once

#include <functional>

#include <Eigen/Core>

#include "../math/pi.hpp"
#include "../math/kahan_sum.hpp"

#include "hawp_paramset.hpp"
#include "shapes/shape_enum.hpp"
#include "shapes/shape_enum_subset.hpp"


namespace waveblocks {
    namespace wavepackets {
        using shapes::ShapeEnum;
        using shapes::ShapeSlice;

        /**
         * \brief Evaluates a wavepacket slice by slice.
         *
         * This class is low-level. You should not use it directly.
         * You should use the high-level member functions
         * AbstractScalarHaWp::evaluate() and AbstractScalarHaWpBasis::evaluate_basis().
         *
         * The only reason you may want to use HaWpEvaluator directly is when
         * you gain an advantage by evaluating a wavepacket slice-by-slice.
         * The slice-by-slice evaluation reduces memory since you
         * don't have to store all basis function values, but only the 'active' ones.
         * Take a look at the implementation of the member functions all() or reduce()
         * to learn, how to evaluate a wave packet slice-by-slice.
         *
         * \tparam D dimensionality of wavepacket
         * \tparam MulitIndex The type used to represent multi-indices.
         * \tparam N
         * Number of quadrature points.
         * Don't use Eigen::Dynamic. It works, but performance is bad.
         */
        template<dim_t D, class MultiIndex, int N>
        class HaWpEvaluator
        {
        public:
            typedef Eigen::Matrix<complex_t,D,D> CMatrixDD;
            typedef Eigen::Array<complex_t,1,N> CArray1N;
            typedef Eigen::Matrix<complex_t,D,N> CMatrixDN;
            typedef Eigen::Matrix<real_t,D,N> RMatrixDN;

        private:
            real_t eps_;
            const HaWpParamSet<D>* parameters_;
            const ShapeEnum<D,MultiIndex>* enumeration_;

            /**
             * number of quadrature points
             */
            int npts_;

            /**
             * precomputed expression: x - q
             */
            CMatrixDN dx_;

            /**
             * precomputed expression: Q^{-1}
             */
            CMatrixDD Qinv_;

            /**
             * precomputed expression: Q^H * Q^{-T}
             */
            CMatrixDD Qh_Qinvt_;

            /**
             * precomputed expression: Q^{-1} * (x - q)
             */
            CMatrixDN Qinv_dx_;

            /**
             * lookup-table for sqrt
             */
            std::vector<real_t> sqrt_;

        public:
            /**
             * \param[in] eps The semi-classical scaling parameter \f$ \varepsilon \f$ of the wavepacket.
             * \param[in] parameters The Hagedorn parameter set \f$ \Pi \f$ of the wavepacket.
             * \param[in] enumeration
             * The basis shape \f$ \mathfrak{K} \f$ of the wavepacket.
             * Or the union of all component's basis shapes, if you want to evaluate a vectorial wavepacket.
             * \param[in] x Quadrature points: Complex matrix of shape \f$ (D \times N) \f$, where \f$ N \f$ number of quadrature points).
             */
            HaWpEvaluator(real_t eps,
                          const HaWpParamSet<D>* parameters,
                          const ShapeEnum<D,MultiIndex>* enumeration,
                          const CMatrixDN &x)
                : eps_(eps)
                , parameters_(parameters)
                , enumeration_(enumeration)
                , npts_(x.cols())
                , sqrt_()
            {
                RMatrix<D,1> const& q = parameters_->q();
                CMatrix<D,D> const& Q = parameters_->Q();

                // precompute ...
                dx_ = x.colwise() - q.template cast<complex_t>();
                Qinv_ = Q.inverse();
                Qh_Qinvt_ = Q.adjoint()*Qinv_.transpose();
                Qinv_dx_ = Qinv_*dx_;

                //precompute sqrt lookup table
                {
                    int limit = 0;
                    for (dim_t d = 0; d < D; d++)
                        limit = std::max(limit, enumeration_->limit(d) );

                    sqrt_.resize(limit+2);

                    for (int i = 0; i <= limit+1; i++)
                        sqrt_[i] = std::sqrt( real_t(i) );
                }
            }

            /**
             * \brief Evaluates basis function on node \f$ \underline{0} \f$ (ground-state).
             *
             * \return Complex 2D-Array of shape \f$ (D \times N) \f$, where \f$ N \f$ is the number of quadrature points.
             */
            CArray1N seed() const
            {
                RMatrix<D,1> const& p = parameters_->p();
                CMatrix<D,D> const& P = parameters_->P();

                CMatrixDN P_Qinv_dx = P*Qinv_dx_;

                CArray1N pr1 = ( dx_.array() * P_Qinv_dx.array() ).colwise().sum();
                CArray1N pr2 = ( p.transpose()*dx_ ).array();

                CArray1N e = complex_t(0.0, 1.0)/(eps_*eps_) * (0.5*pr1 + pr2);

                return e.exp() / std::pow(math::pi<real_t>()*eps_*eps_, D/4.0);
            }

            /**
             * \brief Having basis function values on previous and current slice,
             * this member function computes basis function values
             * on the next slice (using recursive evaluation formula).
             *
             * Hint: Use function seed() to bootstrap recursion.
             *
             * \param[in] islice Ordinal of current slice.
             * \param[in] prev_basis
             * Basis values on previous slice.
             * Type: Complex 2D-array of shape \f$ (S \times N) \f$,
             * where \f$ S \f$ is the number of nodes in the current slice
             * and \f$ N \f$ is the number of quadrature points.
             * \param[in] curr_basis B
             * Basis values on current slice.
             * Type: Complex 2D-array of shape \f$ (S \times N) \f$
             *
             * \return
             * Computed basis values on next slice.
             * Type: Complex 2D-Array of shape \f$ (S^+ \times N) \f$,
             * where \f$ S^+ \f$ is the number of nodes in the next slice.
             */
            HaWpBasisVector<N> step(std::size_t islice,
                                    const HaWpBasisVector<N>& prev_basis,
                                    const HaWpBasisVector<N>& curr_basis) const
            {
                auto & prev_enum = enumeration_->slice(islice-1);
                auto & curr_enum = enumeration_->slice(islice);
                auto & next_enum = enumeration_->slice(islice+1);

                assert ((int)prev_enum.size() == prev_basis.rows());
                assert ((int)curr_enum.size() == curr_basis.rows());

                HaWpBasisVector<N> next_basis(next_enum.size(), npts_);

                #pragma omp parallel
                {
                    // pre-allocate
                    CArray<1,N> pr1(1,npts_), pr2(1,npts_);

                    //loop over all multi-indices within next slice [j = position of multi-index within next slice]
                    #pragma omp for
                    for (std::size_t j = 0; j < next_enum.size(); j++) {
                        std::array<int,D> next_index = next_enum[j];
                        //find valid precursor: find first non-zero entry
                        dim_t axis = D;
                        for (dim_t d = 0; d < D; d++) {
                            if (next_index[d] != 0) {
                                axis = d;
                                break;
                            }
                        }

                        assert(axis != D); //assert that multi-index contains some non-zero entries


                        // compute contribution of current slice
                        std::array<int,D> curr_index = next_index;
                        curr_index[axis] -= 1; //get backward neighbour
                        std::size_t curr_ordinal = curr_enum.find(curr_index);

                        assert(curr_ordinal < curr_enum.size()); //assert that multi-index has been found within current slice

                        pr1 = curr_basis.row(curr_ordinal) * Qinv_dx_.row(axis).array() * std::sqrt(2.0)/eps_ ;


                        // compute contribution of previous slice
                        std::array< std::size_t,D > prev_ordinals = prev_enum.find_backward_neighbours(curr_index);

                        pr2.setZero();

                        for (dim_t d = 0; d < D; d++) {
                            if (curr_index[d] != 0) {
                                pr2 += prev_basis.row(prev_ordinals[d]) * Qh_Qinvt_(axis,d) * sqrt_[ curr_index[d] ];
                            }
                        }


                        // compute basis value within next slice
                        next_basis.row(j) = (pr1 - pr2) / sqrt_[ 1+curr_index[axis] ];
                    }

                }
                return next_basis;
            }

            /**
             * \brief Evaluates all basis functions.
             *
             * \return
             * Complex 2D-Array of shape \f$ (|\mathfrak{K}| \times N) \f$,
             * where \f$ N \f$ is the number of quadrature points.
             */
            HaWpBasisVector<N> all() const
            {
                HaWpBasisVector<N> complete_basis(enumeration_->n_entries(), npts_);
                //HaWpBasisVector<N> complete_basis = HaWpBasisVector<N>::Zero(enumeration_->n_entries(), npts_);
                //~ std::cout << "complete_basis (init):\n" << complete_basis << "\n";
                //~ std::cout << "complete_basis.shape: " << complete_basis.rows() << ", " << complete_basis.cols() << "\n";

                HaWpBasisVector<N> prev_basis(0,npts_);
                HaWpBasisVector<N> curr_basis(0,npts_);
                HaWpBasisVector<N> next_basis(1,npts_);

                complete_basis.block(0, 0, 1, npts_) = next_basis = seed();

                for (int islice = 0; islice < enumeration_->n_slices(); islice++) {
                    prev_basis = std::move(curr_basis);
                    curr_basis = std::move(next_basis);

                    next_basis = step(islice, prev_basis, curr_basis);

                    std::size_t offset = enumeration_->slice(islice+1).offset();

                    //~ std::cout << "offset: " << offset << ", next_basis.rows: " << next_basis.rows() << "\n";

                    complete_basis.block(offset, 0, next_basis.rows(), npts_) = next_basis;
                }

                return complete_basis;
            }

            /**
             * \brief Evaluates wavepacket in a memory efficient manner.
             *
             * This function computes the dot product (thus the name 'reduce') of the wavepacket coefficients
             * with the wavepacket basis. This is done by evaluating the wavepacket basis slice by slice
             * and multiplying the basis values with the coefficients on the fly. Computed
             * basis values are discarded, once they are not needed any more.
             *
             * \param[in] coefficients Vector of wavepacket coefficients, length is \f$ |\mathfrak{K}| \f$.
             * \return
             * Complex 2D-Array of shape \f$ (1 \times N) \f$,
             * where \f$ N \f$ is the number of quadrature points.
             */
            CArray<1,N> reduce(const Coefficients& coefficients) const
            {
                // use Kahan's algorithm to accumulate bases with O(1) numerical error instead of O(Sqrt(N))
                math::KahanSum< CArray<1,N> > psi( CArray<1,N>::Zero(1,npts_) );

                HaWpBasisVector<N> prev_basis(0,npts_);
                HaWpBasisVector<N> curr_basis(0,npts_);
                HaWpBasisVector<N> next_basis(1,npts_);

                next_basis = seed();

                psi += next_basis.row(0)*coefficients[0];

                for (int islice = 0; islice < enumeration_->n_slices(); islice++) {
                    prev_basis = std::move(curr_basis);
                    curr_basis = std::move(next_basis);

                    next_basis = step(islice, prev_basis, curr_basis);

                    std::size_t offset = enumeration_->slice(islice+1).offset();

                    for (long j = 0; j < next_basis.rows(); j++) {
                        complex_t cj = coefficients[offset + j];

                        //prints: multi-index -> basis -> coefficient
                        //std::cout << enumeration->slice(islice+1)[j] << " -> " << next_basis.row(j).matrix() << " * " << cj << std::endl;

                        psi += next_basis.row(j)*cj;
                    }
                }

                return psi();
            }

            /**
             * \brief Efficiently evaluates a vectorial wavepacket with shared Hagedorn parameter set,
             * but different basis shapes.
             *
             * The union of all component's basis shapes is passed to the constructor.
             *
             * This function is used to evaluate a homogeneous wavepacket (see HomogeneousHaWp::evaluate())
             *
             * \param[in] subset_enums The shape enumerations of each wavepacket component.
             * \param[in] subset_coeffs The coefficients of each wavepacket component.
             * \param[in] n_components The number of wavepacket components \f$ C \f$.
             * \return
             * Complex 2D-Array of shape \f$ (C \times N) \f$,
             * where \f$ N \f$ is the number of quadrature points.
             */
            CArray<Eigen::Dynamic,N> vector_reduce(ShapeEnum<D,MultiIndex> ** subset_enums,
                                                   complex_t const ** subset_coeffs,
                                                   std::size_t n_components) const
            {
                // use Kahan's algorithm to accumulate bases with O(1) numerical error instead of O(Sqrt(N))
                std::vector< math::KahanSum< CArray<1,N> > > psi(n_components);

                HaWpBasisVector<N> prev_basis(0,npts_);
                HaWpBasisVector<N> curr_basis(0,npts_);
                HaWpBasisVector<N> next_basis(1,npts_);

                next_basis = seed();

                for (std::size_t n = 0; n < n_components; n++) {
                    psi[n] = math::KahanSum< CArray<1,N> >( CArray<1,N>::Zero(1,npts_)); // zero initialize
                    psi[n] += subset_coeffs[n][0]*next_basis.row(0).matrix();
                }

                for (int islice = 0; islice < enumeration_->n_slices(); islice++) {
                    prev_basis = std::move(curr_basis);
                    curr_basis = std::move(next_basis);

                    next_basis = step(islice, prev_basis, curr_basis);

                    ShapeSlice<D,MultiIndex> const& superset_slice = enumeration_->slice(islice+1);

                    //             for (std::size_t n = 0; n < n_components; n++) {
                    //                 ShapeSlice<D,MultiIndex> const& subset_slice = subset_enums[n]->slice(islice+1);
                    //                 HaWpBasisVector<N> subset_basis = shape_enum::copy_subset(next_basis, superset_slice, subset_slice);
                    //
                    //                 for (long j = 0; j < subset_basis.rows(); j++) {
                    //                     complex_t cj = subset_coeffs[n][subset_slice.offset() + j];
                    //
                    //                     psi[n] += cj*subset_basis.row(j);
                    //                 }
                    //             }

                    std::vector<std::size_t> seek(n_components);
                    for (long j = 0; j < next_basis.rows(); j++) {
                        for (std::size_t n = 0; n < n_components; n++) {
                            ShapeSlice<D,MultiIndex> const& subset_slice = subset_enums[n]->slice(islice+1);

                            if (seek[n] < subset_slice.size() && subset_slice[seek[n]] == superset_slice[j]) {
                                complex_t cj = subset_coeffs[n][subset_slice.offset() + seek[n]];

                                psi[n] += next_basis.row(j)*cj;

                                seek[n]++;
                            }
                        }
                    }
                }

                CArray<Eigen::Dynamic,N> result(n_components, npts_);
                for (std::size_t n = 0; n < n_components; n++) {
                    result.row(n) = psi[n]();
                }

                return result;
            }

            /**
             * \brief Efficiently evaluates a vectorial wavepacket with shared Hagedorn parameter set and shared basis shapes.
             *
             * This function is used to evaluate a wavepacket gradient (see HaWpGradient::evaluate()).
             *
             * \param[in] coefficients The coefficients of each wavepacket component.
             * \param[in] n_components The number of wavepacket components.
             * \param[in] n_components The number of wavepacket components \f$ C \f$.
             * \return
             * Complex 2D-Array of shape \f$ (C \times N) \f$,
             * where \f$ N \f$ is the number of quadrature points.
             */
            CArray<Eigen::Dynamic,N> vector_reduce(complex_t const ** coefficients, std::size_t n_components) const
            {
                // use Kahan's algorithm to accumulate bases with O(1) numerical error instead of O(Sqrt(N))
                std::vector< math::KahanSum< CArray<1,N> > > psi(n_components);

                HaWpBasisVector<N> prev_basis(0,npts_);
                HaWpBasisVector<N> curr_basis(0,npts_);
                HaWpBasisVector<N> next_basis(1,npts_);

                next_basis = seed();

                for (std::size_t n = 0; n < n_components; n++) {
                    psi[n] = math::KahanSum< CArray<1,N> >( CArray<1,N>::Zero(1,npts_)); // zero initialize
                    psi[n] += coefficients[n][0]*next_basis.row(0);
                }

                for (int islice = 0; islice < enumeration_->n_slices(); islice++) {
                    prev_basis = std::move(curr_basis);
                    curr_basis = std::move(next_basis);

                    next_basis = step(islice, prev_basis, curr_basis);

                    ShapeSlice<D,MultiIndex> const& slice = enumeration_->slice(islice+1);


                    for (long j = 0; j < next_basis.rows(); j++) {
                        for (std::size_t n = 0; n < n_components; n++) {
                            complex_t cj = coefficients[n][slice.offset() + j];

                            psi[n] += next_basis.row(j)*cj;
                        }
                    }
                }

                CArray<Eigen::Dynamic,N> result(n_components, npts_);
                for (std::size_t n = 0; n < n_components; n++) {
                    result.row(n) = psi[n]();
                }

                return result;
            }
        };
    }
}
