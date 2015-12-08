#pragma once

#include "shape_enum.hpp"
#include "shape_enum_union.hpp"


namespace waveblocks {

namespace shape_enum {

template<class MultiIndex>
std::vector<MultiIndex> _extend(const std::vector<MultiIndex>& source, dim_t start, dim_t len)
{
    assert (len != 0);

    if (len == 1) {
        std::vector<MultiIndex> next(source);
        for (std::size_t i = 0; i < source.size(); i++) {
            next[i][start] += 1;
        }
        return next;
    } else {
        // use divide and conquer approach
        std::vector<MultiIndex> lhs = _extend(source, start, len/2);
        std::vector<MultiIndex> rhs = _extend(source, start + len/2, len - len/2);

        std::vector<MultiIndex> sink(lhs.size() + rhs.size());

        auto seek = strict_union(lhs.begin(), lhs.end(), rhs.begin(), rhs.end(), sink.begin(), std::less<MultiIndex>{});

        sink.resize(seek - sink.begin());

        return sink;
    }
}

template<dim_t D, class MultiIndex>
ShapeSlice<D, MultiIndex> _extend(const ShapeSlice<D, MultiIndex>& slice, std::size_t offset)
{
    std::vector<MultiIndex> result = _extend(slice._table(), 0, D);

    return {std::move(result), offset};
}

/**
 * \ingroup ShapeExtension
 *
 * \brief For a given enumerated shape, enumerate its extension.
 *
 * \param source pointer to basic shape enumeration
 * \return enumerated shape extension
 */
template<dim_t D, class MultiIndex>
ShapeEnum<D, MultiIndex> extend(const ShapeEnum<D, MultiIndex>* source)
{
    std::vector< ShapeSlice<D, MultiIndex> > slices(source->n_slices()+1);
    std::size_t offset = 1;

    slices[0] = source->slice(0);
    for (int islice = 0; islice < source->n_slices(); islice++) {
        slices[islice+1] = _extend(source->slice(islice), offset);
        offset += slices[islice+1].size();
    }

    MultiIndex limits = source->limits();
    for (dim_t d = 0; d < D; d++) {
        limits[d] += 1;
    }

    return {std::move(slices), offset, limits};
}

} // namespace shape_enum

} // namespace waveblocks
