#ifndef WAVEBLOCKS_SHAPE_ENUMERATOR_HPP
#define WAVEBLOCKS_SHAPE_ENUMERATOR_HPP

#include <vector>

#include "basic_types.hpp"
#include "shape_enum.hpp"

namespace waveblocks {

template<dim_t D, class MultiIndex>
class ShapeEnumerator
{
public:
    /**
     * \brief Enumerates all nodes contained in a shape.
     * 
     * \param[in] shape description of a shape
     * \return shape enumeration
     * 
     * \see HyperCubicShape
     * \see HyperbolicCutShape
     * \see LimitedHyperbolicCutShape
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
                for (dim_t i = 0; i <= shape.template limit<MultiIndex>(index,D-1); i++) {
                    index[D-1] = i;
                    
//                     if (use_dict)
//                         mindices[islice+i].dict_[index] = mindices[islice+i].table_.size();
                    
                    mindices[islice+i].push_back(index);
                }
                index[D-1] = 0;
                
                // iterate over other axes
                if (D > 1) {
                    dim_t j = D-2;
                    while ((int)index[j] == shape.template limit<MultiIndex>(index,j)) {
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
    
    template<class Shape>
    ShapeEnum<D,MultiIndex> enumerate(const Shape& shape) const
    {
        return generate<Shape>(shape);
    }
};

}

#endif