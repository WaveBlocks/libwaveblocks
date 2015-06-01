#ifndef HAGEDORN_BASIS_VECTOR_HPP
#define HAGEDORN_BASIS_VECTOR_HPP

#include <memory>

#include "basic_types.hpp"
#include "shape_enumeration_base.hpp"

namespace waveblocks {

/**
 * \tparam D dimension of shape
 * \tparam C numbers of components if wavepacket is vector-valued
 */
template<dim_t D, int C = 1>
class HagedornBasisVector
{
private:
    
    
public:
    std::shared_ptr< std::vector<  complex_t > > basis_;
    std::shared_ptr< ShapeEnumeration<D> > enumeration_;
    
    HagedornBasisVector &operator=()
    {
        
    }
    
    complex_t operator[] (std::size_t ordinal) const
    {
        return basis_[ordinal];
    }
};

}

#endif