#pragma once

#include <memory>
#include <array>

#include "../../basic_types.hpp"

#include "shape_enum.hpp"


namespace waveblocks {
    namespace wavepackets {
        namespace shapes {
            namespace shape_enum {
                template<dim_t D>
                bool _copy_subset__fast_equals(const std::array<int,D>& lhs, const std::array<int,D>& rhs)
                {
                    // multi-indices are lexically sorted therefore first entries will almost always be identical
                    // this function compares entries beginning on last entries

                    for (dim_t i = D; i > 0; i--) {
                        if (lhs[i-1] != rhs[i-1])
                            return false;
                    }
                    return true;
                }

                /**
                 * \attention This function shows \e undefined \e behaviour if \p subset_slice is not a subset of \p superset_slice.
                 *
                 * \tparam D dimension of multi-index
                 * \tparam T component data type of \p superset_data and return value
                 * \tparam N number of quadrature points or \e Eigen::Dynamic
                 * \param[in] superset_data (number of nodes in superset, number of quadrature points)-matrix
                 * \param[in] superset_slice nodes within superset slice
                 * \param[in] offset1 Offset within the basis shape
                 * \param[in] subset_data (number of nodes in subset, number of quadrature points)-matrix
                 * \param[in] subset_slice nodes within subset slice
                 * \param[in] offset2 Offset within the basis shape
                 * \return (number of nodes in subset, number of quadrature points)-matrix
                 */
                template<dim_t D, class MultiIndex, int N>
                void copy_subset(const HaWpBasisVector<N>& superset_data, std::size_t offset1,
                                 HaWpBasisVector<N>& subset_data, std::size_t offset2,
                                 const ShapeSlice<D,MultiIndex>& superset_slice,
                                 const ShapeSlice<D,MultiIndex>& subset_slice)
                {
                    auto superset_it = superset_slice.begin();
                    auto subset_it = subset_slice.begin();

                    while (superset_it != superset_slice.end() && subset_it != subset_slice.end()) {
                        if (*superset_it == *subset_it) {
                            std::size_t sub = offset2 + (subset_it - subset_slice.begin());
                            std::size_t super = offset1 + (superset_it - superset_slice.begin());

                            subset_data.row(sub) = superset_data.row(super);
                            ++subset_it;
                        }
                        ++superset_it;
                    }
                }

                template<dim_t D, class MultiIndex, int N>
                HaWpBasisVector<N> copy_subset(const HaWpBasisVector<N>& superset_data,
                                               const ShapeSlice<D,MultiIndex>& superset_slice,
                                               const ShapeSlice<D,MultiIndex>& subset_slice)
                {
                    HaWpBasisVector<N> subset_data(subset_slice.size(), superset_data.cols());
                    copy_subset(superset_data, 0, subset_data, 0, superset_slice, subset_slice);
                    return subset_data;
                }

                template<dim_t D, class MultiIndex, int N>
                HaWpBasisVector<N> copy_subset(const HaWpBasisVector<N>& superset_data,
                                               const ShapeEnum<D,MultiIndex>* superset_enum,
                                               const ShapeEnum<D,MultiIndex>* subset_enum)
                {
                    HaWpBasisVector<N> subset_data(subset_enum->n_entries(), superset_data.cols());

                    std::size_t i_off = 0;
                    std::size_t o_off = 0;
                    for (int islice = 0; islice < subset_enum->n_slices(); islice++) {
                        copy_subset(superset_data, i_off, subset_data, o_off, superset_enum->slice(islice), subset_enum->slice(islice));
                        i_off += superset_enum->slice(islice).size();
                        o_off += subset_enum->slice(islice).size();
                    }

                    return subset_data;
                }
            }
        }
    }
}
