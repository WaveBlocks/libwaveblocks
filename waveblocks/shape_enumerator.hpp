#ifndef WAVEBLOCKS_SHAPE_ENUMERATOR_HPP
#define WAVEBLOCKS_SHAPE_ENUMERATOR_HPP

#include <vector>
#include <memory>

#include "basic_types.hpp"
#include "shape_enum.hpp"
#include "shape_base.hpp"

namespace waveblocks {

/**
 * \brief Enumerates nodes of a basis shape.
 * 
 * This basis shape enumerator takes shape description of the form
 * AbstractShape and converts the information to a ShapeEnum.
 * 
 * \tparam D The basis shape dimensionality.
 * \tparam MultiIndex The type to represent a multi-index. TinyMultiIndex is a valid type.
 */
template<dim_t D, class MultiIndex>
class ShapeEnumerator
{
public:
    /**
     * \deprecated Use member function enumerate instead.
     */
    template<class Shape>
    ShapeEnum<D,MultiIndex> generate(const Shape& shape) const
    {
        std::vector< std::vector<MultiIndex> > mindices;
        
        // check multi-index type for compatibility
        {
            MultiIndex test;
            for (dim_t d = 0; d < D; d++) {
                test[d] = shape.bbox(d);
                if (test[d] != shape.bbox(d)) {
                    throw std::runtime_error("multi-index type is not suitable. reason: overflow");
                }
            }
        }
        
        // initialize slice vectors
        {
            std::size_t sum = 0;
            for (dim_t d = 0; d < D; d++) {
                sum += shape.bbox(d)+1;
            }
            
            mindices.resize(sum);
        }
        
        // enumerate shape & store all multi-indices
        {
            MultiIndex index{}; //zero initialize
            std::size_t islice = 0;
            
            while (true) {
                // iterate over last axis
                for (dim_t i = 0; i <= shape.limit(static_cast<std::array<int,D> >(index).data(),D-1); i++) {
                    index[D-1] = i;
                    
//                     if (use_dict)
//                         mindices[islice+i].dict_[index] = mindices[islice+i].table_.size();
                    
                    mindices[islice+i].push_back(index);
                }
                index[D-1] = 0;
                
                // iterate over other axes
                if (D > 1) {
                    dim_t j = D-2;
                    while ((int)index[j] == shape.limit(static_cast<std::array<int,D> >(index).data(),j)) {
                        islice -= index[j];
                        index[j] = 0;
                        if (j == 0)
                            goto enumeration_complete;
                        else
                            j = j-1;
                    }
                    islice += 1;
                    index[j] += 1;
                }
                else break;
            }
enumeration_complete:
            (void)0;
        }
        
        std::vector< ShapeSlice<D,MultiIndex> > slices(mindices.size());
        
        std::size_t offset = 0; // number of entries in all previous slices
        std::size_t islice = 0; // ordinal of current slice
        for (; islice < mindices.size() && mindices[islice].size() != 0; islice++) {
            ShapeSlice<D,MultiIndex> slice{std::move(mindices[islice]), offset};
            offset += slice.size();
            
            slices[islice] = std::move(slice);
        }
        
        slices.resize(islice);
        
        // total number of entries in all slices
        size_t size = offset;
        
        MultiIndex limits;
        for (dim_t d = 0; d < D; d++)
            limits[d] = shape.bbox(d);
        
        return {std::move(slices), size, limits};
    }
    
    /**
     * \brief Enumerates all nodes of a basis shape described by AbstractShape.
     * 
     * \param shape A reference to the basis shape description.
     * \return A shared pointer to the enumerated shape.
     */
    std::shared_ptr<ShapeEnum<D,MultiIndex> > enumerate(AbstractShape<D> const& shape) const
    {
        return std::make_shared<ShapeEnum<D,MultiIndex> >(generate<AbstractShape<D> >(shape));
    }
    
    /**
     * \brief Enumerates all nodes of basis shape described by AbstractShape.
     * 
     * \param shape A pointer to the basis shape description.
     * \return A shared pointer to the enumerated shape.
     */
    std::shared_ptr<ShapeEnum<D,MultiIndex> > enumerate(AbstractShape<D> const* shape) const
    {
        return std::make_shared<ShapeEnum<D,MultiIndex> >(generate<AbstractShape<D> >(*shape));
    }
};

}

#endif
