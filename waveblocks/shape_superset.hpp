#ifndef WAVEBLOCKS_SHAPE_SUPERSET
#define WAVEBLOCKS_SHAPE_SUPERSET

#include <algorithm>
#include <tuple>

#include "waveblocks/basic_types.hpp"

namespace waveblocks {

template<dim_t D, class S1, class... SS>
class SupersetShape
{
private:
    S1 first_;
    SupersetShape<D,SS...> rest_;
    
public:
    SupersetShape(const S1& first, const SS&... rest)
        : first_(first)
        , rest_(rest...)
    { }
    
    inline int bbox(dim_t axis) const
    {
        return std::max((int)first_.bbox(axis), rest_.bbox(axis));
    }
    
    template<class MultiIndex>
    inline int limit(const MultiIndex &coordinate, dim_t axis) const
    {
        return std::max( (int)first_.template limit<MultiIndex>(coordinate, axis),
                         rest_.template limit<MultiIndex>(coordinate, axis) );
    }
};

template<dim_t D, class S>
class SupersetShape<D,S>
{
private:
    S last_;
    
public:
    SupersetShape(const S& last)
        : last_(last)
    { }
    
    inline int bbox(dim_t axis) const
    {
        return (int)last_.bbox(axis);
    }
    
    template<class MultiIndex>
    inline int limit(const MultiIndex &coordinate, dim_t axis) const
    {
        return (int)last_.template limit<MultiIndex>(coordinate, axis);
    }
};

}

#endif