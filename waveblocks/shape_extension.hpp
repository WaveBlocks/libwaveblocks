#ifndef WAVEBLOCKS_SHAPE_EXTENSION_HPP
#define WAVEBLOCKS_SHAPE_EXTENSION_HPP

#include <vector>
#include <memory>

#include "basic_types.hpp"
#include "multi_index.hpp"

namespace waveblocks {

template<dim_t D, class S>
class ShapeExtension
{
private:
    S shape_;
    
public:
    ShapeExtension(S shape)
        : shape_(shape)
    { }
    
    int getSurface(dim_t axis, const MultiIndex<D> &index) const
    {
        int value = shape_.getSurface(axis, index);
        
        //extend only when there actually are some nodes
        if (value > 0)
            value += 1;
        
        for (dim_t d = 0; d < D; d++) {
            if (d != axis && index[d] != 0) {
                //backward neighbour index
                MultiIndex<D> prev_index = index; prev_index[d] -= 1;
                
                value = std::max(value, shape_.getSurface(axis, prev_index));
            }
        }
        
        return value;
    }
};

template<dim_t D, class S>
class ShapeExtensionEnumeration
{
private:
    S shape_;
    ShapeExtension<D,S> extension_;
    std::shared_ptr<std::vector<MultiIndex<D>>> table_;
    
public:
    ShapeExtensionEnumeration(S shape)
        : shape_(shape)
        , extension_(shape)
        , table_(std::make_shared<std::vector<MultiIndex<D>>>())
    {
        MultiIndex<D> index = {{}};
        
        //get first extension node
        index[D-1] = 1 + shape_.getSurface(D-1,index);
        
        while (true) {
            for (dim_t i = 1 + shape_.getSurface(D-1,index); i <= extension_.getSurface(D-1,index); i++) {
                index[D-1] = i;
                table_->push_back(index);
            }
            index[D-1] = 0;
            
            if (D > 1) {
                dim_t j = D-2;
                while ((int)index[j] == 1 + shape_.getSurface(j, index)) {
                    index[j] = 0;
                    if (j == 0)
                        return;
                    else
                        j = j-1;
                }
                index[j] += 1;
            }
        }
    }
    
    ShapeExtensionEnumeration(const ShapeExtensionEnumeration &that)
        : shape_(that.shape_)
        , extension_(that.extension_)
        , table_(that.table_)
    { }
    
    ShapeExtensionEnumeration &operator=(const ShapeExtensionEnumeration &that)
    {
        shape_ = that.shape_;
        extension_ = that.extension_;
        table_ = that.table_;
        return *this;
    }
    
    
};

}

#endif