#ifndef WAVEBLOCKS_SHAPE_ENUMERATION_SUBSET_HPP
#define WAVEBLOCKS_SHAPE_ENUMERATION_SUBSET_HPP

#include <memory>

#include "basic_types.hpp"
#include "shape_enumeration_base.hpp"

namespace waveblocks {

template<dim_t D>
class ShapeEnumerationSubset : public ShapeEnumeration<D>
{
private:
    std::shared_ptr< ShapeEnumeration<D> > superset_;
    std::shared_ptr< ShapeEnumeration<D> > subset_;
    
public:
    
};

}

#endif