#ifndef WAVEBLOCKS_LEXICAL_SHAPE_ENUMERATOR_HPP
#define WAVEBLOCKS_LEXICAL_SHAPE_ENUMERATOR_HPP

#include <cassert>
#include <memory>
#include <string>
#include <sstream>
#include <iostream>
#include <vector>

namespace waveblocks {

template<dim_t D, class MultiIndex, class S>
class LexicalIndexGenerator
{
private:
    S shape_;
    
    MultiIndex index_;
    
public:
    LexicalIndexGenerator(S shape)
        : shape_(shape)
        , index_()
    { }
    
    LexicalIndexGenerator(S shape, MultiIndex index)
        : shape_(shape)
        , index_(index)
    { }
    
    LexicalIndexGenerator(const LexicalIndexGenerator &that)
        : shape_(that.shape_)
        , index_(that.index_)
    { }
    
    LexicalIndexGenerator &operator=(const LexicalIndexGenerator &that)
    {
        shape_ = that.shape_;
        index_ = that.index_;
        
        return *this;
    }
    
    MultiIndex index() const
    {
        return index_;
    }
    
    bool backward()
    {
        //find first non-zero entry
        std::size_t axis = D-1;
        while (index_[axis] == 0) {
            if (axis == 0)
                return false;
            else
                --axis;
        }
        
        index_[axis++] -= 1;
        for (; axis < D; axis++) {
            index_[axis] = shape_.template limit< MultiIndex >(index_, axis);
        }
        
        return true;
    }
    
    bool forward()
    {
        std::size_t axis = D-1;
        while ((int)index_[axis] == shape_.template limit< MultiIndex >(index_, axis)) {
            index_[axis] = 0;
            if (axis == 0)
                return false; //end of chain
            else
                --axis;
        }
        
        index_[axis] += 1;
        
        return true;
    }
};

}

#endif