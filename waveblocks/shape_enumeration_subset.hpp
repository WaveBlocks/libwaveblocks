#ifndef WAVEBLOCKS_SHAPE_ENUMERATION_SUBSET_HPP
#define WAVEBLOCKS_SHAPE_ENUMERATION_SUBSET_HPP

#include <memory>
#include <array>

#include <Eigen/Core>

#include "basic_types.hpp"
#include "shape_enum.hpp"

namespace waveblocks {

template<dim_t D, class MultiIndex>
bool _copy_subset__fast_equals(const MultiIndex& a, const MultiIndex& b)
{
    // since multi-indices are lexically sorted, first entries will almost always be equals
    // therefore this function compares entries beginning on last entries
    
    for (dim_t i = D; i > 0; i--) {
        if (a[i-1] != b[i-1])
            return false;
    }
    return true;
}

/**
 * \attention This function shows \e undefined \e behaviour if \p subset_enum is not a subset of \p superset_enum.
 * 
 * \tparam D dimension of multi-index
 * \tparam T component data type of \p superset_data and return value
 * \tparam N number of quadrature points or \e Eigen::Dynamic
 * \param[in] superset_data (#nodes in superset, #quadrature points)-matrix
 * \param[in] superset_enum nodes within superset slice
 * \param[in] subset_enum nodes within subset slice
 * \return (#nodes in subset, #quadrature points)-matrix
 */
template<dim_t D, class MultiIndex, class T, int N>
Eigen::Array<T, Eigen::Dynamic,N> copy_subset(const Eigen::Array<T, Eigen::Dynamic,N>& superset_data, 
                                               const ShapeSlice<D,MultiIndex>& superset_enum, 
                                               const ShapeSlice<D,MultiIndex>& subset_enum)
{
    if (&superset_enum == &subset_enum)
        return superset_data;
    
    Eigen::Array<T, Eigen::Dynamic,N> subset_data{subset_enum.size(), superset_data.cols()};
    
    auto superset_it = superset_enum.begin();
    auto subset_it = subset_enum.begin();
    
    std::size_t superset_ord = 0;
    std::size_t subset_ord = 0;
    
    while (superset_it != superset_enum.end() && subset_it != subset_enum.end()) {
        
        while ( !_copy_subset__fast_equals<D,MultiIndex>(*superset_it, *subset_it) ) {
            ++superset_it;
            ++superset_ord;
        }
        
        subset_data.row(subset_ord) = superset_data.row(superset_ord);
        ++subset_it;
        ++subset_ord;
    }
    
    return subset_data;
}

}

#endif